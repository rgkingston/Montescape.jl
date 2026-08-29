# not enabled in prototype released with thesis
#concept to be as followed
# Pairwise scenario difference reports.
#
# Given two completed MC runs A and B and a shared endpoint, writes:
#   - LCC probability difference GeoTIFF (B - A)
#   - diff_lcc_<endpoint>_stats.txt     (summary stats of the diff)
#   - bottleneck_<endpoint>.txt         (count + coordinates of bottleneck cells)
#


using Statistics
using Printf

"""
    compare_scenarios(scenario_a_dir, scenario_b_dir, endpoint, output_dir) -> Dict

Returns a Dict with optional keys:
  - `lcc_diff_raster`  : path to the diff GeoTIFF
  - `lcc_diff_report`  : path to the diff stats txt
  - `bottleneck_report`: path to the bottleneck txt
"""
function compare_scenarios(scenario_a_dir::AbstractString,
                           scenario_b_dir::AbstractString,
                           endpoint::AbstractString,
                           output_dir::AbstractString)
    mkpath(output_dir)
    results = Dict{String,Any}()

    prob_a = joinpath(scenario_a_dir, "monte_carlo", "aggregated", endpoint, "lcc_probability_$(endpoint).tif")
    prob_b = joinpath(scenario_b_dir, "monte_carlo", "aggregated", endpoint, "lcc_probability_$(endpoint).tif")

    if isfile(prob_a) && isfile(prob_b)
        ra = read_raster(prob_a)
        rb = read_raster(prob_b)
        a = copy(ra.data); b = copy(rb.data)
        nd = ra.nodata
        a[a .== nd] .= NaN32
        b[b .== nd] .= NaN32
        diff = b .- a

        diff_out = copy(diff)
        diff_out[isnan.(diff_out)] .= nd
        diff_raster = joinpath(output_dir, "diff_lcc_$(endpoint).tif")
        write_raster(diff_raster, diff_out, ra)
        results["lcc_diff_raster"] = diff_raster

        report_path = joinpath(output_dir, "diff_lcc_$(endpoint)_stats.txt")
        _write_diff_report(report_path; diff=diff, endpoint=endpoint,
                            scenario_a=scenario_a_dir, scenario_b=scenario_b_dir)
        results["lcc_diff_report"] = report_path
    end

    mean_b = joinpath(scenario_b_dir, "monte_carlo", "aggregated", endpoint, "cs_mean_current_$(endpoint).tif")
    std_b  = joinpath(scenario_b_dir, "monte_carlo", "aggregated", endpoint, "cs_std_current_$(endpoint).tif")
    if isfile(mean_b) && isfile(std_b)
        rm_ = read_raster(mean_b); rs_ = read_raster(std_b)
        m = copy(rm_.data); s = copy(rs_.data); nd = rm_.nodata
        m[m .== nd] .= NaN32
        s[s .== nd] .= NaN32
        m_finite = filter(isfinite, m); s_finite = filter(isfinite, s)
        isempty(m_finite) && return results
        m_thr = quantile(m_finite, 0.90)
        s_thr = quantile(s_finite, 0.25)

        bottleneck_idxs = Tuple{Int,Int,Float64,Float64}[]
        nrows, ncols = size(m)
        @inbounds for j in 1:ncols, i in 1:nrows
            mv = m[i, j]; sv = s[i, j]
            if isfinite(mv) && isfinite(sv) && mv > m_thr && sv < s_thr
                push!(bottleneck_idxs, (i, j, Float64(mv), Float64(sv)))
            end
        end
        report_path = joinpath(output_dir, "bottleneck_$(endpoint).txt")
        _write_bottleneck_report(report_path; bottlenecks=bottleneck_idxs,
                                  mean_threshold=m_thr, std_threshold=s_thr,
                                  endpoint=endpoint,
                                  scenario_b=scenario_b_dir,
                                  grid_size=(nrows, ncols))
        results["bottleneck_report"] = report_path
    end

    return results
end

function _write_diff_report(path::AbstractString;
                             diff::AbstractMatrix,
                             endpoint::AbstractString,
                             scenario_a::AbstractString,
                             scenario_b::AbstractString)
    mkpath(dirname(path))
    finite = filter(isfinite, diff)
    n_finite = length(finite)
    if n_finite == 0
        open(path, "w") do io
            println(io, "# LCC probability diff (B - A)")
            println(io, "endpoint:   $endpoint")
            println(io, "scenario_a: $scenario_a")
            println(io, "scenario_b: $scenario_b")
            println(io, "status:     EMPTY")
        end
        return path
    end
    pos = count(>(0.05f0), finite)
    neg = count(<(-0.05f0), finite)
    open(path, "w") do io
        println(io, "# LCC probability diff (B - A)")
        println(io, "endpoint:    $endpoint")
        println(io, "scenario_a:  $scenario_a")
        println(io, "scenario_b:  $scenario_b")
        println(io)
        println(io, "key\tvalue")
        @printf(io, "n_finite\t%d\n",                    n_finite)
        @printf(io, "min\t%.6f\n",                       minimum(finite))
        @printf(io, "max\t%.6f\n",                       maximum(finite))
        @printf(io, "mean\t%.6f\n",                      mean(finite))
        @printf(io, "median\t%.6f\n",                    median(finite))
        @printf(io, "std\t%.6f\n",                       std(finite))
        @printf(io, "abs_mean\t%.6f\n",                  mean(abs.(finite)))
        @printf(io, "cells_gained_>0.05\t%d\n",          pos)
        @printf(io, "cells_lost_<-0.05\t%d\n",           neg)
        @printf(io, "fraction_gained\t%.6f\n",           pos / n_finite)
        @printf(io, "fraction_lost\t%.6f\n",             neg / n_finite)
    end
    return path
end

function _write_bottleneck_report(path::AbstractString;
                                   bottlenecks::AbstractVector,
                                   mean_threshold::Real,
                                   std_threshold::Real,
                                   endpoint::AbstractString,
                                   scenario_b::AbstractString,
                                   grid_size::Tuple{Int,Int})
    mkpath(dirname(path))
    nrows, ncols = grid_size
    n = length(bottlenecks)
    sorted = sort(bottlenecks; by = t -> -t[3])  # highest mean current first
    cx = isempty(bottlenecks) ? NaN : mean(b[1] for b in bottlenecks)
    cy = isempty(bottlenecks) ? NaN : mean(b[2] for b in bottlenecks)
    open(path, "w") do io
        println(io, "# Bottleneck cells (high mean current AND low std)")
        println(io, "endpoint:        $endpoint")
        println(io, "scenario_b:      $scenario_b")
        println(io, "grid_rows:       $nrows")
        println(io, "grid_cols:       $ncols")
        @printf(io, "mean_threshold:  %.6f  (p90 of finite mean current)\n", mean_threshold)
        @printf(io, "std_threshold:   %.6f  (p25 of finite std current)\n",  std_threshold)
        println(io, "n_bottlenecks:   $n")
        @printf(io, "fraction_of_grid: %.6f\n", n / (nrows * ncols))
        if !isempty(bottlenecks)
            @printf(io, "centroid_row:    %.2f\n", cx)
            @printf(io, "centroid_col:    %.2f\n", cy)
        end
        println(io)
        println(io, "# Top-50 bottleneck cells by mean current")
        println(io, "rank\trow\tcol\tmean_current\tstd_current")
        for (k, b) in enumerate(sorted[1:min(50, end)])
            @printf(io, "%d\t%d\t%d\t%.6f\t%.6f\n", k, b[1], b[2], b[3], b[4])
        end
    end
    return path
end
