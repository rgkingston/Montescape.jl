# Pure-Julia anisotropic Dijkstra LCP on a friction raster.
#
# Edge cost = average friction of two cells * Euclidean distance between
# their centres * res * (optional Tobler-style anisotropy multiplier based
# on the DEM-derived slope along the edge).
#
# Indexing: every matrix here is [row, col] (Julia native). The IO layer
# in io.jl handles transposing to/from GDAL's [col, row] convention.
#
# This is the key fix vs. r.walk: every Monte Carlo iteration uses the
# friction surface as-passed and does NOT silently re-bake slope.

using DataStructures

const _NEIGH_8 = (
    (-1, -1, sqrt(2.0)), (-1, 0, 1.0), (-1, 1, sqrt(2.0)),
    ( 0, -1, 1.0),                     ( 0, 1, 1.0),
    ( 1, -1, sqrt(2.0)), ( 1, 0, 1.0), ( 1, 1, sqrt(2.0)),
)

# Knight-move (2,1) offsets, used when cfg.knight = true. Helps smooth
# staircase artefacts. Distance is sqrt(5).
const _NEIGH_KNIGHT = (
    (-2, -1, sqrt(5.0)), (-2, 1, sqrt(5.0)),
    ( 2, -1, sqrt(5.0)), ( 2, 1, sqrt(5.0)),
    (-1, -2, sqrt(5.0)), ( 1, -2, sqrt(5.0)),
    (-1,  2, sqrt(5.0)), ( 1,  2, sqrt(5.0)),
)

"""
    _buffer_cells(rc, rows, cols, radius_px) -> Vector{Tuple{Int,Int}}

Return every (row, col) within `radius_px` of `rc`, clipped to the raster.
Used to seed Dijkstra from a small disc rather than a single pixel, which
matches the GRASS r.walk buffer_m semantics from the original pipeline.
"""
function _buffer_cells(rc::Tuple{Int,Int}, rows::Int, cols::Int, radius_px::Int)
    out = Tuple{Int,Int}[]
    pr, pc = rc
    rad2 = radius_px * radius_px
    rmin = max(1, pr - radius_px); rmax = min(rows, pr + radius_px)
    cmin = max(1, pc - radius_px); cmax = min(cols, pc + radius_px)
    @inbounds for r in rmin:rmax, c in cmin:cmax
        if (r - pr)^2 + (c - pc)^2 <= rad2
            push!(out, (r, c))
        end
    end
    isempty(out) && push!(out, (pr, pc))
    return out
end

"""
    accumulate_costs(friction, start_rc; res, nodata, knight=true, dem=nothing)
        -> (cum_cost::Matrix{Float64}, came_from::Matrix{Int})

Dijkstra cost accumulation from a single source cell or a set of source
cells (start_rc may be `(row, col)` or `Vector{Tuple{Int,Int}}` for
buffered start regions).

`came_from[r, c]` stores the linear index of the parent cell, or 0 if root.
If `dem` is provided, edge costs are scaled by the configured slope-dependent
cost function using the signed elevation gradient along each edge, allowing
uphill and downhill movement to incur different costs
"""
function accumulate_costs(friction::AbstractMatrix{<:Real},
                          start_rc;
                          res::Real = 1.0,
                          nodata::Real = -9999.0,
                          knight::Bool = true,
                          dem::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
                          slope_cfg::AbstractDict = Dict())
    rows, cols = size(friction)
    cum = fill(Inf, rows, cols)
    came = zeros(Int, rows, cols)

    pq = BinaryMinHeap{Tuple{Float64,Int,Int}}()

    sources = start_rc isa Tuple{Int,Int} ? [start_rc] : collect(start_rc)
    for (r, c) in sources
        if 1 <= r <= rows && 1 <= c <= cols && friction[r, c] != nodata
            cum[r, c] = 0.0
            push!(pq, (0.0, r, c))
        end
    end

    offsets = knight ? Iterators.flatten((_NEIGH_8, _NEIGH_KNIGHT)) : _NEIGH_8

    while !isempty(pq)
        d, r, c = pop!(pq)
        d > cum[r, c] && continue

        f_here = Float64(friction[r, c])
        dem_here = dem === nothing ? 0.0 : Float64(dem[r, c])

        for (dr, dc, dist) in offsets
            nr = r + dr; nc = c + dc
            (1 <= nr <= rows && 1 <= nc <= cols) || continue
            f_n = Float64(friction[nr, nc])
            f_n == nodata && continue

            edge = 0.5 * (f_here + f_n) * dist * res
            if dem !== nothing
                Δh = Float64(dem[nr, nc]) - dem_here
                slope = Δh / (dist * res)
                # edge costs are scaled by the user-configured slope cost function
                # (see friction.slope.method in the YAML)
                edge *= edge_slope_cost(slope, slope_cfg)
            end
            nd = d + edge
            if nd < cum[nr, nc]
                cum[nr, nc] = nd
                came[nr, nc] = (r - 1) * cols + c
                push!(pq, (nd, nr, nc))
            end
        end
    end
    return cum, came
end

"""
    trace_path(came_from, target_rc) -> Vector{Tuple{Int,Int}}

Reconstruct path as a sequence of (row, col) from the target back to the
source. Returns an empty vector if target is unreachable.
"""
function trace_path(came_from::AbstractMatrix{Int}, target_rc::Tuple{Int,Int})
    rows, cols = size(came_from)
    r, c = target_rc
    path = Tuple{Int,Int}[]
    if !(1 <= r <= rows && 1 <= c <= cols) || came_from[r, c] == 0 && (r, c) != target_rc
        # unreachable: came_from of target is 0 and target itself is not the source
        # (the source has came_from==0 but was explicitly seeded)
    end
    while r >= 1 && c >= 1
        push!(path, (r, c))
        idx = came_from[r, c]
        idx == 0 && break
        r = div(idx - 1, cols) + 1
        c = mod(idx - 1, cols) + 1
    end
    reverse!(path)
    return path
end

"""
    run_lcp(friction, ref, source_xy, ground_xy, run_dir; cfg, name,
            waypoints=[], dem=nothing, buffer_m=0.0) -> Dict

Compute LCP between two world coordinates. Optional waypoints force the
path through intermediate points (computed as concatenated segments).
When `buffer_m > 0`, the Dijkstra wavefront is seeded from every cell
within `buffer_m` of each anchor, matching the GRASS r.walk buffer
semantics. Returns a dict with file paths, the path as cells, and the
total cost.
"""
function run_lcp(friction::AbstractMatrix{<:Real},
                 ref::RasterData,
                 source_xy::Tuple{<:Real,<:Real},
                 ground_xy::Tuple{<:Real,<:Real},
                 run_dir::AbstractString;
                 cfg::AbstractDict,
                 name::AbstractString = "lcp",
                 waypoints::Vector{<:Tuple{<:Real,<:Real}} = Tuple{Float64,Float64}[],
                 dem::Union{Nothing,AbstractMatrix}=nothing,
                 buffer_m::Real = 0.0)

    knight = get(get(cfg, "lcp", Dict()), "knight", true)
    res = pixel_resolution(ref)
    nd = ref.nodata
    rows, cols = size(friction)
    slope_cfg = get(get(cfg, "friction", Dict()), "slope", Dict())

    nodes = [source_xy, waypoints..., ground_xy]
    rcs = [world_to_pixel(ref, x, y) for (x, y) in nodes]

    rad_px = buffer_m > 0 ? max(Int(ceil(Float64(buffer_m) / res)), 1) : 0

    full_path = Tuple{Int,Int}[]
    total_cost = 0.0

    for k in 1:length(rcs)-1
        start = rad_px > 0 ? _buffer_cells(rcs[k], rows, cols, rad_px) : rcs[k]
        cum, came = accumulate_costs(friction, start;
                                     res=res, nodata=nd, knight=knight, dem=dem,
                                     slope_cfg=slope_cfg)
        # Pick the lowest-cost cell within the target buffer as the actual endpoint
        target = rcs[k+1]
        if rad_px > 0
            best = target
            best_cost = Inf
            for (tr, tc) in _buffer_cells(target, rows, cols, rad_px)
                if cum[tr, tc] < best_cost
                    best_cost = cum[tr, tc]
                    best = (tr, tc)
                end
            end
            target = best
        end
        seg = trace_path(came, target)
        isempty(seg) &&
            throw(ErrorException("LCP unreachable between waypoint $k and $(k+1)"))
        if k == 1
            append!(full_path, seg)
        else
            append!(full_path, seg[2:end])
        end
        total_cost += cum[target...]
    end

    out_dir = joinpath(run_dir, "lcp")
    mkpath(out_dir)
    raster = zeros(Float32, size(friction))
    for (r, c) in full_path
        raster[r, c] = 1.0f0
    end
    raster_path = joinpath(out_dir, "lcp_$(name).tif")
    write_raster(raster_path, raster, ref)

    # Also write a CSV of the path as world coords
    csv_path = joinpath(out_dir, "lcp_$(name).csv")
    open(csv_path, "w") do io
        println(io, "order,row,col,x,y")
        for (i, (r, c)) in enumerate(full_path)
            x, y = pixel_to_world(ref, r, c)
            println(io, "$i,$r,$c,$x,$y")
        end
    end

    # And a real GeoPackage LineString so the user can drop the path
    # straight into QGIS as a vector layer (rather than reconstructing
    # from the CSV). Failure here is non-fatal: the raster + CSV are
    # already on disk.
#    gpkg_path = joinpath(out_dir, "lcp_$(name).gpkg") #can be generated from the csv files - will be implemented in a future patch
#    try
#        _write_line_gpkg(gpkg_path, full_path, ref; layer_name=name)
#    catch err
#        @warn "  [lcp] vector export failed; raster + CSV are fine" err
#        gpkg_path = ""
#    end

    return Dict(
        "raster" => raster_path,
        "csv"    => csv_path,
#       "vector" => "", #can be generated from the csv files - will be implemented in a future patch
        "cells"  => full_path,
        "cost"   => total_cost,
    )
end

# Write a single LineString through the centre of each path cell to a
# GeoPackage in the reference raster's CRS. Used by run_lcp for the
# vector export. ArchGDAL is already loaded via src/io.jl.
function _write_line_gpkg(path::AbstractString,
                           cells::AbstractVector,
                           ref::RasterData;
                           layer_name::AbstractString="lcp")
    isempty(cells) && return path
    AG = Base.invokelatest(getfield, Main, :ArchGDAL)  # already in scope via io.jl
    isfile(path) && rm(path; force=true)
    AG.create(path; driver=AG.getdriver("GPKG")) do ds
        srs = AG.importWKT(ref.projection)
        layer = AG.createlayer(name=layer_name, dataset=ds, geom=AG.wkbLineString, spatialref=srs)
        AG.createfeature(layer) do feat
            line = AG.createlinestring()
            for (r, c) in cells
                x, y = pixel_to_world(ref, r, c)
                AG.addpoint!(line, x, y)
            end
            AG.setgeom!(feat, line)
        end
    end
    return path
end
