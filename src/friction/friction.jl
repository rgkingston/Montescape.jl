# Combine enabled layers into a weighted friction surface.

using ArchGDAL
const _AG = ArchGDAL

"""
    build_friction_surface(cfg, run_dir; dem_override=nothing,
                           weight_override=nothing, write=true) -> (path, array, ref_raster)

Build the weighted friction surface. Returns the output path (or empty
string if write=false), the friction array, and the reference RasterData
(for downstream geotransform/projection use).

- dem_override: alternative DEM Matrix{Float32}, used inside Monte Carlo
  to substitute a perturbed DEM without re-reading the file.
- weight_override: Dict{String,Float64} of layer_name => new_weight, used
  inside Monte Carlo to substitute sampled weights.
"""
function build_friction_surface(cfg::AbstractDict, run_dir::AbstractString;
                                 dem_override::Union{Nothing,AbstractMatrix}=nothing,
                                 weight_override::Union{Nothing,AbstractDict}=nothing,
                                 write::Bool=true)
    fric = cfg["friction"]
    # Prebuilt mode will be in future update
    # in Monte Carlo weight perturbation, so the framework always rebuilds

    dem_path = cfg["data"]["dem_raster"]
    dem_ref = read_raster(dem_path)
    dem = dem_override === nothing ? dem_ref.data : Float32.(dem_override)
    res = pixel_resolution(dem_ref)

    nodata_mask = dem .== dem_ref.nodata
    combined = zeros(Float32, size(dem))
    total_w = 0.0

    function w(name, default_section)
        if weight_override !== nothing && haskey(weight_override, name)
            return Float64(weight_override[name])
        end
        return Float64(get(default_section, "weight", 0.0))
    end

    # Slope
    sl = get(fric, "slope", nothing)
    if sl !== nothing && get(sl, "enabled", false)
        c = compute_slope_friction(dem, res, sl)
        wv = w("slope", sl)
        combined .+= c .* Float32(wv)
        total_w += wv
        @info "  [friction] slope weight=$wv"
    end

    # TRI
    tri_sec = get(fric, "tri", nothing)
    if tri_sec !== nothing && get(tri_sec, "enabled", false)
        c = compute_tri(dem)
        wv = w("tri", tri_sec)
        combined .+= c .* Float32(wv)
        total_w += wv
        @info "  [friction] tri weight=$wv"
    end

    # Aspect
    asp = get(fric, "aspect", nothing)
    if asp !== nothing && get(asp, "enabled", false)
        c = compute_aspect_penalty(dem, res, asp)
        wv = w("aspect", asp)
        combined .+= c .* Float32(wv)
        total_w += wv
        @info "  [friction] aspect weight=$wv"
    end

    # Vegetation
    veg = get(fric, "vegetation", nothing)
    if veg !== nothing && get(veg, "enabled", false)
        arr = _load_and_align(veg["raster"], dem_path, size(dem))
        c = compute_vegetation(arr, veg)
        wv = w("vegetation", veg)
        combined .+= c .* Float32(wv)
        total_w += wv
        @info "  [friction] vegetation weight=$wv"
    end

    # Hydrology
    hyd = get(fric, "hydrology", nothing)
    if hyd !== nothing && get(hyd, "enabled", false)
        arr = _load_and_align(hyd["raster"], dem_path, size(dem))
        c = compute_hydrology(arr, hyd)
        wv = w("hydrology", hyd)
        combined .+= c .* Float32(wv)
        total_w += wv
        @info "  [friction] hydrology weight=$wv"
    end

    # Custom (DEM-derived) layers, plugin-registered.
    # NOTE: custom layer files are include()'d at runtime, so the functions
    # they define live in a newer Julia world than this caller. Calling them
    # directly would raise a world-age MethodError. Base.invokelatest forces
    # a lookup against the latest world.
    for layer in get(fric, "custom_layers", Any[])
        get(layer, "enabled", false) || continue
        name = String(layer["name"])
        fn = lookup_layer(name)
        c = Base.invokelatest(fn, dem, res, layer)
        wv = w(name, layer)
        combined .+= Float32.(c) .* Float32(wv)
        total_w += wv
        @info "  [friction] $name weight=$wv"
    end

    if total_w > 0
        if !(abs(total_w - 1.0) < 0.01)
            throw(ArgumentError("Friction weights must sum to 1.0 (got $total_w)."))
        end
        combined ./= Float32(total_w)
    end

    # Ensure strictly positive friction so Dijkstra works
    combined .= max.(combined, 0.001f0)
    combined[nodata_mask] .= dem_ref.nodata

    out_path = ""
    if write
        outdir = joinpath(run_dir, "friction")
        mkpath(outdir)
        out_path = joinpath(outdir, get(fric, "output_name", "friction_surface.tif"))
        write_raster(out_path, combined, dem_ref)
        @info "  [friction] Saved: $out_path"
    end
    return out_path, combined, dem_ref
end

# Reproject + resample raster onto the DEM grid.
function _load_and_align(path::AbstractString, ref_path::AbstractString, shape::Tuple{Int,Int})
    out = zeros(Float32, shape)
    _AG.read(ref_path) do refds
        _AG.read(path) do srcds
            # Use gdalwarp via ArchGDAL to reproject onto ref grid
            warped = _AG.unsafe_gdalwarp([srcds],
                ["-t_srs", _AG.getproj(refds),
                 "-te", string.(_corner_extent(refds))...,
                 "-tr", string(abs(_AG.getgeotransform(refds)[2])),
                        string(abs(_AG.getgeotransform(refds)[6])),
                 "-r", "bilinear",
                 "-of", "MEM"])
            band = _AG.getband(warped, 1)
            arr = _AG.read(band)
            out .= Float32.(arr)
        end
    end
    return out
end

function _corner_extent(ds)
    gt = _AG.getgeotransform(ds)
    w = _AG.width(ds)
    h = _AG.height(ds)
    xmin = gt[1]
    xmax = gt[1] + w * gt[2]
    ymax = gt[4]
    ymin = gt[4] + h * gt[6]
    return (xmin, ymin, xmax, ymax)
end
