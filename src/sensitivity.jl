# Morris elementary-effects global sensitivity (Campolongo 2007).
#
# Lightweight, defensible global sensitivity analysis for the friction-weight
# inputs. Cheaper than Sobol (~10x a single MC run vs ~1000x) but still
# detects non-linear effects and input interactions, unlike one-at-a-time
# local sensitivity.
#
# The summary statistic is a scalar reduction of the iteration output. By
# default we use the LCP total cost from the iteration result dict; the
# user can pass any (cfg, iter_result) -> Float64 function.

using Random
using Statistics
using Printf

"""
    compute_sensitivity_morris(cfg, run_dir, endpoint;
                                n_trajectories=10,
                                n_levels=4,
                                objective=:lcp_cost,
                                seed=42) -> Union{Nothing,String}

Run Morris elementary-effects on the enabled friction weights for one
endpoint. Returns the CSV path or `nothing` if there is nothing to do.

The sampling design samples `n_trajectories` radial trajectories; each
trajectory perturbs one input at a time by `delta = n_levels/(2*(n_levels-1))`
in the unit hypercube. For k enabled weights we run `n_trajectories * (k+1)`
pipeline iterations. With k = 6 layers and 10 trajectories that's 70 runs.

Objective:
  :lcp_cost              - total cost of the iteration's LCP
  :passes_top_fraction   - max pass-crossing fraction across passes
  custom callable f(cfg, iter_result) -> Float64
"""
function compute_sensitivity_morris(cfg::AbstractDict,
                                     run_dir::AbstractString,
                                     endpoint::AbstractDict;
                                     n_trajectories::Integer=10,
                                     n_levels::Integer=4,
                                     objective=:lcp_cost,
                                     seed::Integer=42)
    base = _enabled_weight_base(cfg)
    isempty(base) && (@info "  [morris] no enabled friction weights, skipping"; return nothing)
    names = collect(keys(base))
    k = length(names)
    delta = n_levels / (2 * (n_levels - 1))

    epname = endpoint["name"]
    out_dir = joinpath(run_dir, "monte_carlo", "aggregated", epname)
    mkpath(out_dir)

    @info "  [morris] $epname: $n_trajectories trajectories x $(k+1) runs = $(n_trajectories*(k+1)) iterations"
    points = read_points(cfg["data"]["sites_vector"];
                          name_field=get(cfg["data"], "site_name_field", "name"))
    src = find_point(points, endpoint["source_name"])
    gnd = find_point(points, endpoint["ground_name"])

    # ---- inherit the same worker pool MC uses --------------------------
    nw = max(1, Int(parse(Int, get(ENV, "MONTESCAPE_MC_WORKERS",
                                   string(Int(get(get(cfg, "monte_carlo", Dict()),
                                               "n_workers", 1)))))))
    if nw > 1 && nworkers() < nw
        _ensure_workers(nw, cfg)   # no-op if workers already up from MC block
    end
    # --------------------------------------------------------------------

    function _run_trajectory(t)
        rng = MersenneTwister(seed + t)
        x0 = rand(rng, k) .* (1 - delta)
        order = randperm(rng, k)
        row = fill(NaN, k)

        weights = _renormalise(Dict(zip(names, x0)))
        y_prev = _objective_value(cfg, src, gnd, weights, objective;
                                   seed=seed * 1000 + t)
        x = copy(x0)
        for step in 1:k
            i = order[step]
            x[i] = clamp(x[i] + delta, 0.0, 1.0)
            weights = _renormalise(Dict(zip(names, x)))
            y_curr = _objective_value(cfg, src, gnd, weights, objective;
                                       seed=seed * 1000 + t * 100 + step)
            row[i] = (y_curr - y_prev) / delta
            y_prev = y_curr
        end
        return t => row
    end

    pairs = if nw > 1
        wp = WorkerPool(workers())
        pmap(wp, 1:n_trajectories) do t
            _run_trajectory(t)
        end
    else
        [_run_trajectory(t) for t in 1:n_trajectories]
    end

    # Reassemble ees matrix
    ees = fill(NaN, n_trajectories, k)
    for (t, row) in pairs
        ees[t, :] = row
    end

    # --- rest of function unchanged (CSV + TXT writing) ----------------
    csv_path = joinpath(out_dir, "sensitivity_morris.csv")
    open(csv_path, "w") do io
        println(io, "input,mean_ee,mean_abs_ee,std_ee,n_trajectories")
        rows = NamedTuple[]
        for (j, n) in enumerate(names)
            col = filter(isfinite, ees[:, j])
            isempty(col) && continue
            push!(rows, (input=n,
                         mu=mean(col), mu_star=mean(abs.(col)),
                         sigma=length(col) > 1 ? std(col) : 0.0,
                         n=length(col)))
        end
        sort!(rows; by=r->r.mu_star, rev=true)
        for r in rows
            @printf(io, "%s,%.6g,%.6g,%.6g,%d\n",
                    r.input, r.mu, r.mu_star, r.sigma, r.n)
        end
    end

    rows = _read_csv_rows(csv_path)
    if !isempty(rows)
        txt_path = joinpath(out_dir, "sensitivity_morris.txt")
        _write_morris_txt(txt_path, rows, epname)
    end

    @info "  [morris] $epname -> $csv_path"
    return csv_path
end
function _write_morris_txt(path::AbstractString, rows, epname::AbstractString)
    labels = [String(r[1]) for r in rows]
    mu_star = [parse(Float64, r[3]) for r in rows]
    mu      = [parse(Float64, r[2]) for r in rows]
    sigma   = [parse(Float64, r[4]) for r in rows]
    max_mu_star = maximum(mu_star; init=0.0)
    label_w = max(8, maximum(length, labels; init=8))
    open(path, "w") do io
        println(io, "# Morris elementary-effects sensitivity")
        println(io, "endpoint: $epname")
        println(io, "inputs ranked by mean |EE| (higher = more influential)")
        println(io)
        @printf(io, "%-*s  %12s  %12s  %12s  %s\n",
                label_w, "input", "mu*", "mu", "sigma", "bar")
        for i in eachindex(labels)
            bar_len = max_mu_star > 0 ? Int(round(40 * mu_star[i] / max_mu_star)) : 0
            @printf(io, "%-*s  %12.6g  %12.6g  %12.6g  %s\n",
                    label_w, labels[i], mu_star[i], mu[i], sigma[i], "#"^bar_len)
        end
    end
    return path
end

# Helpers ---------------------------------------------------------------

function _enabled_weight_base(cfg::AbstractDict)
    fric = cfg["friction"]
    bw = Dict{String,Float64}()
    for k in ("slope", "tri", "aspect", "vegetation", "hydrology")
        sec = get(fric, k, nothing)
        sec === nothing && continue
        if get(sec, "enabled", false)
            bw[k] = Float64(get(sec, "weight", 0.0))
        end
    end
    for layer in get(fric, "custom_layers", Any[])
        if get(layer, "enabled", false)
            bw[String(layer["name"])] = Float64(get(layer, "weight", 0.0))
        end
    end
    return bw
end

function _renormalise(w::AbstractDict)
    total = sum(values(w))
    total <= 0 && return Dict(k => 1.0 / length(w) for k in keys(w))
    return Dict(k => v / total for (k, v) in w)
end

function _objective_value(cfg::AbstractDict, src, gnd, weights::AbstractDict,
                          objective; seed::Integer=0)
    tmp = mktempdir()
    try
        _, friction, ref = build_friction_surface(cfg, tmp;
            weight_override=weights, write=false)
        if objective === :lcp_cost || objective == "lcp_cost"
            res = run_lcp(friction, ref, (src.x, src.y), (gnd.x, gnd.y), tmp;
                          cfg=cfg, name="morris_$(seed)")
            return Float64(get(res, "cost", 0.0))
        elseif objective isa Function
            return Float64(objective(cfg, friction, ref, (src, gnd)))
        else
            throw(ArgumentError("Unsupported objective: $objective"))
        end
    finally
        try; rm(tmp; recursive=true, force=true); catch; end
    end
end

function _read_csv_rows(path::AbstractString)
    lines = readlines(path)
    length(lines) < 2 && return Tuple[]
    rows = Tuple[]
    for ln in lines[2:end]
        toks = split(ln, ',')
        length(toks) >= 4 || continue
        push!(rows, Tuple(toks))
    end
    return rows
end
