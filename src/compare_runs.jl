# Cross-run pass-frequency matrix.
#
# Reads pass_frequencies_<endpoint>.csv from each of N completed MC runs,
# pivots into a (pass × run) matrix and writes:
#   pass_frequency_matrix.csv    (wide matrix)
#   pass_frequency_long.csv      (long-format)
#   pass_frequency_summary.txt   (per-run top passes + per-pass aggregates)
#   pass_frequency_heatmap.png   (single synthesis figure)
#

using DataFrames
using CSV
using Statistics
using Printf
using Plots

"""
    build_pass_frequency_matrix(run_dirs, output_dir;
                                run_labels=nothing, endpoint_filter=nothing,
                                sort_by="max") -> Union{Nothing,String}

Assemble the cross-run pass-frequency matrix. Returns the matrix CSV path.
"""
function build_pass_frequency_matrix(run_dirs::AbstractVector,
                                     output_dir::AbstractString;
                                     run_labels::Union{Nothing,AbstractVector}=nothing,
                                     endpoint_filter::Union{Nothing,AbstractString}=nothing,
                                     sort_by::AbstractString="max")
    mkpath(output_dir)
    labels = run_labels === nothing ? [basename(rstrip(r, '/')) for r in run_dirs] :
                                       collect(String.(run_labels))
    if length(labels) != length(run_dirs)
        throw(ArgumentError("run_labels length ($(length(labels))) ≠ run_dirs length ($(length(run_dirs)))"))
    end

    long_rows = NamedTuple{(:pass_name, :run, :fraction, :n_crossings, :n_total_runs),
                           Tuple{String,String,Float64,Int,Int}}[]

    for (run_dir, label) in zip(run_dirs, labels)
        csvs = _find_pass_freq_csvs(run_dir)
        if isempty(csvs)
            @info "  [compare_runs] No pass_frequencies CSVs under $run_dir"
            continue
        end
        for csv_path in csvs
            df = try
                CSV.read(csv_path, DataFrame)
            catch err
                @warn "  [compare_runs] read failed" csv_path err
                continue
            end
            isempty(df) && continue
            "pass_name" in names(df) || continue
            ep_name = "endpoint" in names(df) && !isempty(df.endpoint) ?
                      String(first(df.endpoint)) : basename(dirname(csv_path))
            if endpoint_filter !== nothing && ep_name != endpoint_filter
                continue
            end
            col_label = endpoint_filter !== nothing ? label : "$(label)__$(ep_name)"
            for row in eachrow(df)
                push!(long_rows, (
                    pass_name = String(row.pass_name),
                    run = col_label,
                    fraction = Float64(row.fraction),
                    n_crossings = haskey(row, :n_crossings) ? Int(row.n_crossings) : 0,
                    n_total_runs = haskey(row, :n_total_runs) ? Int(row.n_total_runs) : 0,
                ))
            end
        end
    end

    if isempty(long_rows)
        @info "  [compare_runs] No data; nothing written."
        return nothing
    end

    long_df = DataFrame(long_rows)
    matrix_df = unstack(long_df, :pass_name, :run, :fraction; combine=mean)

    pass_names = matrix_df.pass_name
    score = if sort_by == "max"
        [maximum(skipmissing(values(row[Not(:pass_name)])); init=0.0) for row in eachrow(matrix_df)]
    elseif sort_by == "mean"
        [mean(skipmissing(values(row[Not(:pass_name)]))) for row in eachrow(matrix_df)]
    elseif sort_by == "name"
        zeros(length(pass_names))
    else
        throw(ArgumentError("Unknown sort_by: $sort_by (use max, mean, name)"))
    end
    order = sort_by == "name" ? sortperm(pass_names) : sortperm(score; rev=true)
    matrix_df = matrix_df[order, :]

    csv_path = joinpath(output_dir, "pass_frequency_matrix.csv")
    CSV.write(csv_path, matrix_df)
    CSV.write(joinpath(output_dir, "pass_frequency_long.csv"), long_df)

    _write_pass_freq_summary(joinpath(output_dir, "pass_frequency_summary.txt"),
                              matrix_df, long_df)

    try
        _render_heatmap(matrix_df, joinpath(output_dir, "pass_frequency_heatmap.png"))
    catch err
        @warn "  [compare_runs] heatmap render failed (CSV outputs are still complete)" err
    end

    @info "  [compare_runs] wrote $(size(matrix_df)) matrix -> $csv_path"
    return csv_path
end

function _find_pass_freq_csvs(run_dir::AbstractString)
    agg = joinpath(run_dir, "monte_carlo", "aggregated")
    isdir(agg) || return String[]
    out = String[]
    for entry in readdir(agg; join=true)
        isdir(entry) || continue
        for f in readdir(entry; join=true)
            if endswith(f, ".csv") && occursin("pass_frequencies_", basename(f))
                push!(out, f)
            end
        end
    end
    return sort(out)
end

function _write_pass_freq_summary(path::AbstractString,
                                   matrix_df::DataFrame,
                                   long_df::DataFrame)
    mkpath(dirname(path))
    runs = string.(names(matrix_df)[2:end])
    passes = string.(matrix_df.pass_name)
    open(path, "w") do io
        println(io, "# Pass-frequency cross-run summary")
        println(io, "n_runs:   $(length(runs))")
        println(io, "n_passes: $(length(passes))")
        println(io)
        println(io, "# Per-pass aggregates (across all runs)")
        println(io, "pass_name\tmean\tmax\tmin\tn_runs_present")
        for row in eachrow(matrix_df)
            vals = collect(skipmissing(values(row[Not(:pass_name)])))
            if isempty(vals)
                @printf(io, "%s\tNA\tNA\tNA\t0\n", row.pass_name)
            else
                @printf(io, "%s\t%.4f\t%.4f\t%.4f\t%d\n",
                        row.pass_name, mean(vals), maximum(vals), minimum(vals), length(vals))
            end
        end
        println(io)
        println(io, "# Per-run top-5 passes")
        for r in runs
            sub = filter(row -> !ismissing(row[r]), eachrow(matrix_df))
            ranked = sort(sub; by = row -> -Float64(row[r]))
            println(io, "## run: $r")
            println(io, "rank\tpass_name\tfraction")
            for (k, row) in enumerate(ranked[1:min(5, end)])
                @printf(io, "%d\t%s\t%.4f\n", k, row.pass_name, Float64(row[r]))
            end
            println(io)
        end
    end
    return path
end

# The one PNG the framework still produces. This is a synthesis chart
# meant for human eyes (supervisor / thesis figure / decision-makers);
# the CSVs next to it remain the canonical analysable data.
function _render_heatmap(df::DataFrame, out_path::AbstractString)
    isempty(df) && return
    runs = string.(names(df)[2:end])
    passes = string.(df.pass_name)
    n_rows, n_cols = length(passes), length(runs)
    data = Matrix{Float64}(undef, n_rows, n_cols)
    for (i, row) in enumerate(eachrow(df))
        for (j, r) in enumerate(runs)
            v = row[r]
            data[i, j] = ismissing(v) ? NaN : Float64(v)
        end
    end
    fig_w = max(600, 60 * n_cols + 400)
    fig_h = max(400, 40 * n_rows + 200)
    plt = heatmap(runs, passes, data;
                  clims=(0, 1), c=:viridis,
                  xrotation=45,
                  xlabel="Run (configuration × start)",
                  ylabel="Candidate pass",
                  title="Pass-crossing frequency across configurations",
                  size=(fig_w, fig_h))
    savefig(plt, out_path)
    return out_path
end
