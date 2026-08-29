# Aggregate Monte Carlo iteration outputs into probability surfaces,
# per-pixel confidence-interval rasters, mean current maps, and pass-
# crossing frequencies.

using Statistics
using DataFrames

"""
    aggregate_monte_carlo(cfg, run_dir, endpoint, results) -> Dict

Reads each iteration's LCP/LCC rasters back from disk, stacks them, and
writes:
  - lcp_probability_<endpoint>.tif
  - lcc_probability_<endpoint>.tif
  - lcc_ci_low_<endpoint>.tif      (5th percentile by default)
  - lcc_ci_high_<endpoint>.tif     (95th percentile by default)
  - cs_mean_current_<endpoint>.tif
  - cs_ci_low_<endpoint>.tif       (Circuitscape lower CI)
  - cs_ci_high_<endpoint>.tif      (Circuitscape upper CI)
  - pass_frequencies_<endpoint>.csv

The percentiles come from monte_carlo.aggregate.confidence_interval
(default 0.95).
"""
function aggregate_monte_carlo(cfg::AbstractDict, run_dir::AbstractString,
                                endpoint::AbstractDict, results::AbstractDict)
    isempty(results) && return Dict()
    ana = String.(get(cfg["monte_carlo"], "analyses", String[]))
    agg_cfg = get(cfg["monte_carlo"], "aggregate", Dict{String,Any}())
    write_prob = Bool(get(agg_cfg, "probability_surface", true))
    ci = Float64(get(agg_cfg, "confidence_interval", 0.95))
    ci = clamp(ci, 0.0, 1.0)
    lo_q = (1.0 - ci) / 2.0
    hi_q = 1.0 - lo_q

    epname = endpoint["name"]
    out_dir = joinpath(run_dir, "monte_carlo", "aggregated", epname)
    mkpath(out_dir)

    summary = Dict{String,Any}()
    ref_raster = nothing

    # ---- LCP probability ---------------------------------------------------
    if "lcp" in ana && write_prob
        first_r = results[minimum(keys(results))]["lcp"]["raster"]
        ref_raster = read_raster(first_r)
        stack = _stack_rasters(results, "lcp", "raster")
        prob = vec_mean(stack)
        path = joinpath(out_dir, "lcp_probability_$(epname).tif")
        write_raster(path, prob, ref_raster)
        summary["lcp_probability"] = path
    end

    # ---- LCC probability + CI bands ---------------------------------------
    if "lcc" in ana
        first_r = results[minimum(keys(results))]["lcc"]["corridor_raster"]
        ref_raster = ref_raster === nothing ? read_raster(first_r) : ref_raster
        stack = _stack_rasters(results, "lcc", "corridor_raster")
        if write_prob
            prob = vec_mean(stack)
            path = joinpath(out_dir, "lcc_probability_$(epname).tif")
            write_raster(path, prob, ref_raster)
            summary["lcc_probability"] = path
        end

        lo = vec_quantile(stack, lo_q)
        hi = vec_quantile(stack, hi_q)
        lo_path = joinpath(out_dir, "lcc_ci_low_$(epname).tif")
        hi_path = joinpath(out_dir, "lcc_ci_high_$(epname).tif")
        write_raster(lo_path, lo, ref_raster)
        write_raster(hi_path, hi, ref_raster)
        summary["lcc_ci_low"]  = lo_path
        summary["lcc_ci_high"] = hi_path
    end

    # ---- Circuitscape mean + CI ------------------------------------------
    if "circuitscape" in ana && get(agg_cfg, "mean_current", true)
        arrays = Matrix{Float32}[]
        for (_, r) in results
            cmap = r["circuitscape"]["curmap"]
            isempty(cmap) && continue
            push!(arrays, Float32.(_read_ascii_grid(cmap)))
        end
        if !isempty(arrays) && ref_raster !== nothing
            mean_cur = mean(arrays)
            mpath = joinpath(out_dir, "cs_mean_current_$(epname).tif")
            write_raster(mpath, mean_cur, ref_raster)
            summary["cs_mean_current"] = mpath

            lo = vec_quantile_list(arrays, lo_q)
            hi = vec_quantile_list(arrays, hi_q)
            lopath = joinpath(out_dir, "cs_ci_low_$(epname).tif")
            hipath = joinpath(out_dir, "cs_ci_high_$(epname).tif")
            write_raster(lopath, lo, ref_raster)
            write_raster(hipath, hi, ref_raster)
            summary["cs_ci_low"]  = lopath
            summary["cs_ci_high"] = hipath
        end
    end

    # ---- Pass-crossing frequencies ----------------------------------------
    passes_path = get(cfg["data"], "passes_vector", nothing)
    if passes_path !== nothing && ref_raster !== nothing && "lcp" in ana
        buffer = Float64(get(get(cfg, "analysis", Dict()), "pass_buffer_m", 2000.0))
        passes = read_points(passes_path; name_field="name")
        rows = NamedTuple[]
        for p in passes
            pr, pc = world_to_pixel(ref_raster, p.x, p.y)
            resu = pixel_resolution(ref_raster)
            rad_px = max(Int(ceil(buffer / resu)), 1)
            hits = 0
            for (_, r) in results
                lcp_arr = read_raster(r["lcp"]["raster"]).data
                if _circle_hits(lcp_arr, pr, pc, rad_px)
                    hits += 1
                end
            end
            push!(rows, (pass=p.name, frequency=hits / length(results), n_hits=hits, n_iter=length(results)))
        end
        df = DataFrame(rows)
        path = joinpath(out_dir, "pass_frequencies_$(epname).csv")
        open(path, "w") do io
            println(io, "pass,frequency,n_hits,n_iter")
            for r in eachrow(df)
                println(io, "$(r.pass),$(r.frequency),$(r.n_hits),$(r.n_iter)")
            end
        end
        summary["pass_frequencies"] = path
    end

    return summary
end

# Stack a list of per-iteration rasters into Vector{Matrix} by reading them
# back from disk one at a time. Used for LCP probability + LCC CI rasters.
function _stack_rasters(results::AbstractDict, key1::AbstractString, key2::AbstractString)
    arrays = Matrix{Float32}[]
    for (_, r) in results
        path = r[key1][key2]
        push!(arrays, read_raster(path).data)
    end
    return arrays
end

vec_mean(stack::Vector{<:AbstractMatrix}) = mean(stack)

function vec_quantile(stack::Vector{<:AbstractMatrix}, q::Real)
    isempty(stack) && return Matrix{Float32}(undef, 0, 0)
    rows, cols = size(stack[1])
    out = zeros(Float32, rows, cols)
    @inbounds for i in 1:rows, j in 1:cols
        vals = [Float32(s[i, j]) for s in stack]
        out[i, j] = Float32(quantile(vals, q))
    end
    return out
end

vec_quantile_list(stack, q) = vec_quantile(stack, q)

function _circle_hits(arr::AbstractMatrix, pr::Int, pc::Int, rad::Int)
    rows, cols = size(arr)
    rmin = max(1, pr - rad); rmax = min(rows, pr + rad)
    cmin = max(1, pc - rad); cmax = min(cols, pc + rad)
    @inbounds for r in rmin:rmax, c in cmin:cmax
        ((r - pr)^2 + (c - pc)^2 <= rad^2) || continue
        arr[r, c] > 0 && return true
    end
    return false
end

function _read_ascii_grid(path::AbstractString)
    lines = readlines(path)
    data_lines = lines[7:end]
    rows = length(data_lines)
    cols = length(split(strip(data_lines[1])))
    arr = zeros(Float32, rows, cols)
    for (i, ln) in enumerate(data_lines)
        toks = split(strip(ln))
        for (j, t) in enumerate(toks)
            arr[i, j] = parse(Float32, t)
        end
    end
    return arr
end
