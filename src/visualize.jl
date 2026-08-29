# Quantitative run report generator.
#
# For each output raster found under a run directory, computes summary
# statistics and writes them to a sibling .txt file. Also writes a
# single `reports/summary.txt` indexing everything and a per-endpoint
# overview. No images are produced; the goal is numbers that can be
# diffed, tested, grep-ed and read on a headless server.

using Statistics
using Printf

"""
    generate_report(run_dir) -> Vector{String}

Scan a run directory for the standard raster outputs and produce
quantitative `.txt` reports under `<run_dir>/reports/`. Returns the
list of report files written.
"""
function generate_report(run_dir::AbstractString)
    reports_dir = joinpath(run_dir, "reports")
    mkpath(reports_dir)
    written = String[]
    index = Tuple{String,String,Dict{String,Any}}[]  # (category, label, stats)

    function _record!(category::AbstractString, label::AbstractString,
                      src_raster::AbstractString, dst_txt::AbstractString)
        try
            stats = _raster_stats(src_raster)
            _write_stats_file(dst_txt; category=category, label=label,
                              source=src_raster, stats=stats)
            push!(written, dst_txt)
            push!(index, (String(category), String(label), stats))
        catch err
            @warn "report failed" src_raster err
        end
    end

    friction_dir = joinpath(run_dir, "friction")
    if isdir(friction_dir)
        for f in readdir(friction_dir; join=true)
            endswith(f, ".tif") || continue
            _record!("friction", "friction", f,
                     joinpath(reports_dir, "friction_stats.txt"))
        end
    end

    lcp_dir = joinpath(run_dir, "lcp")
    if isdir(lcp_dir)
        for f in filter(p -> endswith(p, ".tif"), readdir(lcp_dir; join=true))
            ep = _extract_endpoint(basename(f), "lcp_")
            _record!("lcp", ep, f, joinpath(reports_dir, "$(ep)_lcp_stats.txt"))
        end
    end

    lcc_dir = joinpath(run_dir, "lcc")
    if isdir(lcc_dir)
        for f in filter(p -> endswith(p, ".tif") && startswith(basename(p), "lcc_corridor_"),
                        readdir(lcc_dir; join=true))
            ep = _extract_endpoint(basename(f), "lcc_corridor_")
            _record!("lcc", ep, f, joinpath(reports_dir, "$(ep)_lcc_stats.txt"))
        end
    end

    agg_root = joinpath(run_dir, "monte_carlo", "aggregated")
    if isdir(agg_root)
        for ep in readdir(agg_root; join=false)
            ep_dir = joinpath(agg_root, ep)
            isdir(ep_dir) || continue
            for f in readdir(ep_dir; join=true)
                endswith(f, ".tif") || continue
                bn = basename(f)
                if startswith(bn, "lcp_probability_")
                    _record!("lcp_probability", ep, f,
                             joinpath(reports_dir, "$(ep)_lcp_prob_stats.txt"))
                elseif startswith(bn, "lcc_probability_")
                    _record!("lcc_probability", ep, f,
                             joinpath(reports_dir, "$(ep)_lcc_prob_stats.txt"))
                elseif startswith(bn, "cs_mean_current_")
                    _record!("cs_mean_current", ep, f,
                             joinpath(reports_dir, "$(ep)_cs_mean_current_stats.txt"))
                end
            end
        end
    end

    summary_path = joinpath(reports_dir, "summary.txt")
    _write_summary_index(summary_path, run_dir, index)
    push!(written, summary_path)
    return written
end

# ---------------------------------------------------------------------------

"""
    _raster_stats(path) -> Dict{String,Any}

Load a raster and compute summary statistics over its finite (non-nodata)
cells. Returned keys: n_total, n_finite, n_nonzero, min, max, mean, median,
p90, p99, std, sum, fraction_nonzero.
"""
function _raster_stats(path::AbstractString)
    r = read_raster(path)
    arr = copy(r.data)
    arr[arr .== r.nodata] .= NaN32
    finite = filter(isfinite, arr)
    n_total = length(arr)
    n_finite = length(finite)
    if n_finite == 0
        return Dict{String,Any}(
            "n_total"=>n_total, "n_finite"=>0, "n_nonzero"=>0,
            "min"=>NaN, "max"=>NaN, "mean"=>NaN, "median"=>NaN,
            "p90"=>NaN, "p99"=>NaN, "std"=>NaN, "sum"=>0.0,
            "fraction_nonzero"=>0.0,
        )
    end
    nz = count(!iszero, finite)
    return Dict{String,Any}(
        "n_total"          => n_total,
        "n_finite"         => n_finite,
        "n_nonzero"        => nz,
        "min"              => Float64(minimum(finite)),
        "max"              => Float64(maximum(finite)),
        "mean"             => Float64(mean(finite)),
        "median"           => Float64(median(finite)),
        "p90"              => Float64(quantile(finite, 0.90)),
        "p99"              => Float64(quantile(finite, 0.99)),
        "std"              => Float64(std(finite)),
        "sum"              => Float64(sum(finite)),
        "fraction_nonzero" => Float64(nz / n_finite),
    )
end

function _write_stats_file(path::AbstractString;
                            category::AbstractString,
                            label::AbstractString,
                            source::AbstractString,
                            stats::AbstractDict)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# Raster statistics")
        println(io, "category: $category")
        println(io, "label:    $label")
        println(io, "source:   $source")
        println(io)
        println(io, "key\tvalue")
        for k in ("n_total", "n_finite", "n_nonzero", "fraction_nonzero",
                   "min", "max", "mean", "median", "p90", "p99", "std", "sum")
            v = stats[k]
            if v isa Integer
                @printf(io, "%s\t%d\n", k, v)
            else
                @printf(io, "%s\t%.6g\n", k, v)
            end
        end
    end
    return path
end

function _write_summary_index(path::AbstractString, run_dir::AbstractString,
                               index::AbstractVector)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# MovementFramework run summary")
        println(io, "run_dir: $run_dir")
        println(io, "n_reports: $(length(index))")
        println(io)
        println(io, "category\tlabel\tn_finite\tnonzero_frac\tmean\tmax")
        for (cat, label, stats) in index
            @printf(io, "%s\t%s\t%d\t%.4f\t%.4g\t%.4g\n",
                    cat, label,
                    Int(stats["n_finite"]),
                    Float64(stats["fraction_nonzero"]),
                    Float64(stats["mean"]),
                    Float64(stats["max"]))
        end
    end
    return path
end

function _extract_endpoint(filename::AbstractString, prefix::AbstractString)
    stem = replace(filename, ".tif" => "")
    startswith(stem, prefix) ? stem[length(prefix)+1:end] : stem
end
