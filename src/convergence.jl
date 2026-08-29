# Monte Carlo convergence diagnostic.
#
# Stacks per-iteration rasters at increasing subsets (10%, 20%, ..., 100%),
# computes a probability-of-use surface at each subset, and tracks the RMSE
# between consecutive surfaces. When RMSE drops below CONVERGENCE_THRESHOLD
# the MC ensemble is considered converged (Lewis 2021).
#
# Writes a plain-text quantitative report and returns a structured
# Dict so callers can rely on the output file existing and be tested against
# concrete numeric values.

using Statistics
using Printf

const MIN_RUNS_FOR_CONVERGENCE = 3
const CONVERGENCE_THRESHOLD = 0.01

"""
    check_convergence(mc_dir, analysis="lcc") -> Dict

Walks `<mc_dir>/iter_*` (or `run_*`), picks the representative output
raster for the chosen analysis, and writes a convergence report TXT
+ structured result dict into `<mc_dir>/aggregated/`.

`analysis` ∈ ("lcc", "circuitscape"). LCP is vector-like (a thin line)
and uses pass-frequency tables instead.

The returned dict always contains `convergence_report` pointing at a
file that exists. The legacy key `convergence_plot` is aliased to the
same path so older callers keep working.
"""
function check_convergence(mc_dir::AbstractString;
                            analysis::AbstractString="lcc",
                            out_dir::Union{Nothing,AbstractString}=nothing)
    mc_dir = abspath(mc_dir)
    agg_dir = out_dir !== nothing ? out_dir : joinpath(mc_dir, "aggregated")
    mkpath(agg_dir)
    out_txt = joinpath(agg_dir, "$(analysis)_convergence.txt")

    if analysis == "lcp"
        _write_convergence_stub(out_txt, analysis,
            "LCP is line-like; convergence not applicable. Use pass-frequency table instead.")
        return Dict(
            "converged" => nothing,
            "reason"    => "LCP is line-like; use pass-frequency instead",
            "n_runs"    => 0,
            "convergence_report" => out_txt,
            "convergence_plot"   => out_txt,  # legacy alias
            "recommendation"     => "See pass_frequency_*.csv",
            "final_rmse"         => nothing,
            "rmse_trajectory"    => Float64[],
            "step_sizes"         => Int[],
        )
    end

    iter_dirs = sort(filter(p -> isdir(p), readdir(mc_dir; join=true)))
    iter_dirs = filter(d -> startswith(basename(d), "iter_") || startswith(basename(d), "run_"), iter_dirs)

    raster_paths = String[]
    for d in iter_dirs
        candidates = String[]
        if analysis == "lcc"
            for (root, _, files) in walkdir(d)
                for f in files
                    (startswith(f, "lcc_corridor_") && endswith(f, ".tif")) ||
                    (f == "lcc.tif") ||
                    (startswith(f, "lcc_") && endswith(f, ".tif")) || continue
                    push!(candidates, joinpath(root, f))
                end
            end
        elseif analysis == "circuitscape"
            for (root, _, files) in walkdir(d)
                for f in files
                    occursin("curmap", f) && (endswith(f, ".asc") || endswith(f, ".tif")) || continue
                    push!(candidates, joinpath(root, f))
                end
            end
        end
        isempty(candidates) || push!(raster_paths, first(sort(candidates)))
    end

    n_total = length(raster_paths)
    if n_total < MIN_RUNS_FOR_CONVERGENCE
        _write_convergence_stub(out_txt, analysis,
            "Too few iterations to test convergence. " *
            "Found $n_total $(uppercase(analysis)) output(s); need >= $MIN_RUNS_FOR_CONVERGENCE.")
        return Dict(
            "converged" => nothing,
            "reason"    => "only $n_total runs (< $MIN_RUNS_FOR_CONVERGENCE)",
            "n_runs"    => n_total,
            "convergence_report" => out_txt,
            "convergence_plot"   => out_txt,
            "recommendation"     => "Increase monte_carlo.n_runs and re-run.",
            "final_rmse"         => nothing,
            "rmse_trajectory"    => Float64[],
            "step_sizes"         => Int[],
        )
    end

    step = max(1, n_total ÷ 10)
    steps = collect(max(MIN_RUNS_FOR_CONVERGENCE, step):step:n_total)
    if isempty(steps) || last(steps) != n_total
        push!(steps, n_total)
    end

    mode = analysis == "circuitscape" ? :continuous : :binary

    prev = nothing
    rmses = Float64[]
    step_sizes = Int[]
    for n in steps
        surf = _compute_probability(raster_paths[1:n], mode)
        if prev !== nothing
            diff = surf .- prev
            rmse = sqrt(mean(filter(isfinite, diff .^ 2)))
            push!(rmses, rmse)
            push!(step_sizes, n)
        end
        prev = copy(surf)
    end

    converged = !isempty(rmses) && last(rmses) < CONVERGENCE_THRESHOLD
    final_rmse = isempty(rmses) ? nothing : round(last(rmses); digits=6)
    recommendation = converged ? "Sufficient runs." : "Consider increasing n_runs."

    _write_convergence_report(out_txt;
        analysis=analysis,
        n_runs=n_total,
        threshold=CONVERGENCE_THRESHOLD,
        step_sizes=step_sizes,
        rmses=rmses,
        converged=converged,
        final_rmse=final_rmse,
        recommendation=recommendation,
    )

    return Dict(
        "converged"          => converged,
        "final_rmse"         => final_rmse,
        "n_runs"             => n_total,
        "convergence_report" => out_txt,
        "convergence_plot"   => out_txt,  # legacy alias
        "recommendation"     => recommendation,
        "reason"             => nothing,
        "rmse_trajectory"    => rmses,
        "step_sizes"         => step_sizes,
    )
end

function _write_convergence_stub(path::AbstractString, analysis::AbstractString,
                                  msg::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# Monte Carlo Convergence Report")
        println(io, "analysis: $(uppercase(analysis))")
        println(io, "status:   NOT_APPLICABLE")
        println(io)
        println(io, msg)
    end
    return path
end

function _write_convergence_report(path::AbstractString;
                                    analysis::AbstractString,
                                    n_runs::Integer,
                                    threshold::Real,
                                    step_sizes::AbstractVector{<:Integer},
                                    rmses::AbstractVector{<:Real},
                                    converged::Bool,
                                    final_rmse,
                                    recommendation::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# Monte Carlo Convergence Report")
        println(io, "analysis:        $(uppercase(analysis))")
        println(io, "n_runs:          $n_runs")
        println(io, "threshold:       $threshold  (Lewis 2021 RMSE cutoff)")
        println(io, "converged:       $(converged ? "YES" : "NO")")
        println(io, "final_rmse:      $(final_rmse === nothing ? "NA" : final_rmse)")
        println(io, "recommendation:  $recommendation")
        println(io)
        println(io, "# RMSE trajectory  (one row per consecutive subset comparison)")
        println(io, "step\tn_runs_used\trmse\tbelow_threshold")
        for (i, (n, r)) in enumerate(zip(step_sizes, rmses))
            below = r < threshold ? "yes" : "no"
            @printf(io, "%d\t%d\t%.6f\t%s\n", i, n, r, below)
        end
    end
    return path
end

function _compute_probability(paths::Vector{String}, mode::Symbol)
    arrays = Matrix{Float32}[]
    for p in paths
        arr = _read_raster_for_convergence(p)
        transformed = if mode == :binary
            Float32.(arr .> 0)
        else  # :continuous
            vals = copy(arr)
            finite = filter(isfinite, vals)
            if !isempty(finite)
                p99 = quantile(finite, 0.99)
                if p99 > 0
                    vals = vals ./ Float32(p99)
                end
            end
            vals
        end
        push!(arrays, transformed)
    end
    return mean(arrays)
end

function _read_raster_for_convergence(path::AbstractString)
    if endswith(path, ".asc")
        return _read_ascii_grid_simple(path)
    end
    return read_raster(path).data
end

function _read_ascii_grid_simple(path::AbstractString)
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
