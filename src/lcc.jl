# Least-Cost Corridor: two cost-accumulation surfaces (from source AND from
# ground), summed cell-by-cell. Threshold gives the corridor.
#
# All matrices are indexed [row, col] (Julia native); the IO boundary in
# io.jl handles the GDAL [col, row] transpose.

"""
    run_lcc(friction, ref, source_xy, ground_xy, run_dir; cfg, name,
            dem=nothing, buffer_m=0.0) -> Dict

Returns paths to the summed-cost raster, the binary corridor raster, and
optionally a vector polygon export (export_vector flag).

Threshold modes:
  percent  : value = % above global min (e.g. 5 means corridor = within +5%)
  absolute : value = absolute cost units above global min
"""
function run_lcc(friction::AbstractMatrix{<:Real},
                 ref::RasterData,
                 source_xy::Tuple{<:Real,<:Real},
                 ground_xy::Tuple{<:Real,<:Real},
                 run_dir::AbstractString;
                 cfg::AbstractDict,
                 name::AbstractString = "lcc",
                 dem::Union{Nothing,AbstractMatrix}=nothing,
                 buffer_m::Real = 0.0)

    lcc_cfg = cfg["lcc"]
    res = pixel_resolution(ref)
    nd = ref.nodata
    # LCC has its own knight flag, independent of LCP. Fall back to LCP's
    # only if the LCC-specific value is missing (back-compat with old configs).
    knight = Bool(get(lcc_cfg, "knight",
                       get(get(cfg, "lcp", Dict()), "knight", true)))
    rows, cols = size(friction)
    slope_cfg = get(get(cfg, "friction", Dict()), "slope", Dict())

    src_rc = world_to_pixel(ref, source_xy...)
    gnd_rc = world_to_pixel(ref, ground_xy...)

    rad_px = buffer_m > 0 ? max(Int(ceil(Float64(buffer_m) / res)), 1) : 0
    src_seed = rad_px > 0 ? _buffer_cells(src_rc, rows, cols, rad_px) : src_rc
    gnd_seed = rad_px > 0 ? _buffer_cells(gnd_rc, rows, cols, rad_px) : gnd_rc

    cum_src, _ = accumulate_costs(friction, src_seed;
                                  res=res, nodata=nd, knight=knight, dem=dem,
                                  slope_cfg=slope_cfg)
    cum_gnd, _ = accumulate_costs(friction, gnd_seed;
                                  res=res, nodata=nd, knight=knight, dem=dem,
                                  slope_cfg=slope_cfg)

    summed = cum_src .+ cum_gnd
    finite_mask = isfinite.(summed)
    !any(finite_mask) && throw(ErrorException("LCC: no finite cells between endpoints"))
    gmin = minimum(summed[finite_mask])

    thr = lcc_cfg["threshold"]
    method = get(thr, "method", "percent")
    value  = Float64(get(thr, "value", 5.0))
    cutoff = method == "percent" ? gmin * (1.0 + value / 100.0) : gmin + value

    summed_out = copy(summed)
    summed_out[.!finite_mask] .= Float32(nd)
    summed_path = joinpath(run_dir, "lcc", "lcc_sum_$(name).tif")
    mkpath(dirname(summed_path))
    write_raster(summed_path, Float32.(summed_out), ref)

    corridor = Float32.(finite_mask .& (summed .<= cutoff))
    corridor_path = joinpath(dirname(summed_path), "lcc_corridor_$(name).tif")
    write_raster(corridor_path, corridor, ref)

    out = Dict{String,Any}(
        "sum_raster"      => summed_path,
        "corridor_raster" => corridor_path,
        "global_min"      => gmin,
        "cutoff"          => cutoff,
#        "corridor_vector" => ""
    )

    # Honour lcc.export_vector: polygonise the binary corridor raster
    # into a GeoPackage polygon layer for QGIS use. Failure is non-fatal:
    # the rasters are already on disk.
    #disabled - will be added in future patch like LCP vector
    # if Bool(get(lcc_cfg, "export_vector", false))
    #        vec_path = joinpath(dirname(summed_path), "lcc_corridor_$(name).gpkg")
    #        try
    #            _polygonise_corridor(corridor_path, vec_path, ref; layer_name=name)
    #           out["corridor_vector"] = vec_path
    #       catch err
    #           @warn "  [lcc] vector export failed; raster is fine" err
    #        end
    # end

    return out
end

# Polygonise the binary corridor raster (cells with value > 0) into a
# single multi-polygon layer. Uses GDAL's Polygonize via ArchGDAL.
function _polygonise_corridor(in_tif::AbstractString,
                                out_gpkg::AbstractString,
                                ref::RasterData;
                                layer_name::AbstractString="corridor")
    AG = Base.invokelatest(getfield, Main, :ArchGDAL)
    isfile(out_gpkg) && rm(out_gpkg; force=true)
    AG.read(in_tif) do src_ds
        src_band = AG.getband(src_ds, 1)
        AG.create(out_gpkg; driver=AG.getdriver("GPKG")) do dst_ds
            srs = AG.importWKT(ref.projection)
            layer = AG.createlayer(name=layer_name, dataset=dst_ds,
                                    geom=AG.wkbPolygon, spatialref=srs)
            AG.addfielddefn!(layer, "value", AG.OFTInteger)
            AG.polygonize(src_band, layer, 0)
        end
    end
    return out_gpkg
end
