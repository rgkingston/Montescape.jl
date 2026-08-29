# Per-pass crossing-frequency aggregator (port of framework/pass_frequency.py).
#
# For one MC run, for each candidate pass point: buffer the point by
# pass_buffer_m, then count the fraction of iterations whose LCP raster
# has any "on path" cell inside that buffer.
#
# Writes <run_dir>/monte_carlo/aggregated/<endpoint>/pass_frequencies_<endpoint>.csv.

using Printf

const _DEFAULT_PASS_BUFFER_M = 2000.0
const _DEFAULT_NAME_FIELD = "name"
const _DEFAULT_TYPE_FIELD = "type"
const _DEFAULT_PASS_TYPE_VALUE = "pass"

"""
    compute_pass_frequencies_for_endpoint(cfg, run_dir, endpoint_name; buffer_m, ...)
        -> Union{Nothing,String}

Write pass_frequencies_<endpoint>.csv for one endpoint of one MC run.
Returns the CSV path, or `nothing` if nothing could be computed.
"""
function compute_pass_frequencies_for_endpoint(cfg::AbstractDict,
                                                run_dir::AbstractString,
                                                endpoint_name::AbstractString;
                                                buffer_m::Real=_DEFAULT_PASS_BUFFER_M,
                                                name_field::AbstractString=_DEFAULT_NAME_FIELD,
                                                type_field::AbstractString=_DEFAULT_TYPE_FIELD,
                                                pass_type_value::AbstractString=_DEFAULT_PASS_TYPE_VALUE)
    data = cfg["data"]
    passes_path = get(data, "passes_vector", nothing)
    passes_path = passes_path === nothing ? get(data, "sites_vector", nothing) : passes_path
    if passes_path === nothing || !isfile(passes_path)
        @info "  [pass_freq] No passes_vector / sites_vector found; skipping $endpoint_name"
        return nothing
    end

    ep_dir = joinpath(run_dir, "monte_carlo", endpoint_name)
    isdir(ep_dir) || (@info "  [pass_freq] No MC iter dir for $endpoint_name; skipping"; return nothing)

    iter_dirs = sort(filter(p -> isdir(p) && startswith(basename(p), "iter_"),
                            readdir(ep_dir; join=true)))
    isempty(iter_dirs) && (@info "  [pass_freq] No iter_* dirs under $ep_dir"; return nothing)

    # Pick a reference raster (first iter's LCP) for geotransform
    ref = nothing
    iter_lcp_paths = String[]
    for d in iter_dirs
        lcp_dir = joinpath(d, "lcp")
        isdir(lcp_dir) || continue
        candidates = filter(f -> endswith(f, ".tif"), readdir(lcp_dir; join=true))
        isempty(candidates) && continue
        chosen = first(sort(candidates))
        ref === nothing && (ref = read_raster(chosen))
        push!(iter_lcp_paths, chosen)
    end
    ref === nothing && (@info "  [pass_freq] No LCP rasters found for $endpoint_name"; return nothing)

    # Load passes (filter by type if available)
    passes_all = try
        read_points_full(passes_path; name_field=name_field, type_field=type_field)
    catch
        # Fall back to plain read_points if the file has no type field
        [(name=p.name, x=p.x, y=p.y, type="") for p in read_points(passes_path; name_field=name_field)]
    end
    passes = filter(p -> isempty(p.type) || p.type == pass_type_value, passes_all)
    isempty(passes) && (@info "  [pass_freq] No pass rows in $passes_path"; return nothing)

    res = pixel_resolution(ref)
    rad_px = max(Int(ceil(Float64(buffer_m) / res)), 1)
    n_iter = length(iter_lcp_paths)

    out_dir = joinpath(run_dir, "monte_carlo", "aggregated", endpoint_name)
    mkpath(out_dir)
    out_csv = joinpath(out_dir, "pass_frequencies_$(endpoint_name).csv")

    open(out_csv, "w") do io
        println(io, "pass_name,n_crossings,n_total_runs,fraction,buffer_m,endpoint,x,y")
        rows = NamedTuple[]
        for p in passes
            pr, pc = world_to_pixel(ref, p.x, p.y)
            hits = 0
            for lp in iter_lcp_paths
                arr = read_raster(lp).data
                if _circle_hits_lcp(arr, pr, pc, rad_px)
                    hits += 1
                end
            end
            frac = n_iter > 0 ? hits / n_iter : 0.0
            push!(rows, (name=p.name, hits=hits, n=n_iter,
                         frac=round(frac; digits=4), x=p.x, y=p.y))
        end
        sort!(rows; by=r->r.frac, rev=true)
        for r in rows
            @printf(io, "%s,%d,%d,%.4f,%.1f,%s,%.6f,%.6f\n",
                    r.name, r.hits, r.n, r.frac, Float64(buffer_m),
                    endpoint_name, r.x, r.y)
        end
    end
    @info "  [pass_freq] [$endpoint_name] $(length(passes)) passes x $n_iter iter -> $out_csv"
    return out_csv
end

"""
    compute_pass_frequencies_all_endpoints(cfg, run_dir; buffer_m) -> Dict

Loop over every endpoint in the config. Returns a Dict mapping endpoint
name to the written CSV path (or `nothing`).
"""
function compute_pass_frequencies_all_endpoints(cfg::AbstractDict,
                                                 run_dir::AbstractString;
                                                 buffer_m::Real=_DEFAULT_PASS_BUFFER_M)
    out = Dict{String,Any}()
    for ep in get(cfg, "endpoints", Any[])
        epname = String(ep["name"])
        try
            out[epname] = compute_pass_frequencies_for_endpoint(
                cfg, run_dir, epname; buffer_m=buffer_m)
        catch err
            @warn "  [pass_freq] endpoint '$epname' failed" err
            out[epname] = nothing
        end
    end
    return out
end

function _circle_hits_lcp(arr::AbstractMatrix, pr::Int, pc::Int, rad::Int)
    rows, cols = size(arr)
    rmin = max(1, pr - rad); rmax = min(rows, pr + rad)
    cmin = max(1, pc - rad); cmax = min(cols, pc + rad)
    rad2 = rad * rad
    @inbounds for r in rmin:rmax, c in cmin:cmax
        ((r - pr)^2 + (c - pc)^2 <= rad2) || continue
        arr[r, c] > 0 && return true
    end
    return false
end
