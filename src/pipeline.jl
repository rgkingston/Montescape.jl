# Top-level orchestrator: runs the deterministic single-shot analyses and
# then (if enabled) the Monte Carlo ensemble per endpoint, followed by
# post-MC diagnostics (convergence + pass-frequency + Morris sensitivity
# + metadata + viz).

using Logging

"""
    run_pipeline(config_path) -> Dict

End-to-end execution. Returns a dict summarising all written paths.
"""
function run_pipeline(config_path::AbstractString)
    cfg = load_config(config_path)
    run_dir = create_run_dir(cfg)
    @info "==> Run dir: $run_dir"

    custom_files = get(cfg["friction"], "custom_layer_files", String[])
    load_custom_layer_files!(custom_files, cfg["_project_root"])

    summary = Dict{String,Any}("run_dir" => run_dir)

    @info "==> Building friction surface (single-shot)"
    fric_path, friction, ref = build_friction_surface(cfg, run_dir)
    summary["friction"] = fric_path

    points = read_points(cfg["data"]["sites_vector"];
                          name_field=get(cfg["data"], "site_name_field", "name"))
    passes = nothing
    if get(cfg["data"], "passes_vector", nothing) !== nothing
        passes = read_points(cfg["data"]["passes_vector"]; name_field="name")
    end

    summary["endpoints"] = Dict{String,Any}()

    for ep in cfg["endpoints"]
        epname = ep["name"]
        @info "==> Endpoint: $epname"
        ep_summary = Dict{String,Any}()

        src = find_point(points, ep["source_name"])
        gnd = find_point(points, ep["ground_name"])
        wps = Tuple{Float64,Float64}[]
        if passes !== nothing
            for wn in get(ep, "waypoints", String[])
                p = find_point(passes, wn)
                push!(wps, (p.x, p.y))
            end
        end
        buf = Float64(get(ep, "buffer_m", 0.0))

        if get(cfg["lcp"], "enabled", false)
            ep_summary["lcp"] = run_lcp(friction, ref, (src.x, src.y), (gnd.x, gnd.y),
                                          run_dir; cfg=cfg, name=epname,
                                          waypoints=wps, buffer_m=buf)
        end
        if get(cfg["lcc"], "enabled", false)
            ep_summary["lcc"] = run_lcc(friction, ref, (src.x, src.y), (gnd.x, gnd.y),
                                          run_dir; cfg=cfg, name=epname,
                                          buffer_m=buf)
        end
        if get(cfg["circuitscape"], "enabled", false)
            ep_summary["circuitscape"] = run_circuitscape(friction, ref,
                (src.x, src.y), (gnd.x, gnd.y), run_dir; cfg=cfg, name=epname)
        end

        if get(cfg["monte_carlo"], "enabled", false)
            @info "==> Monte Carlo: $epname"
            mc_results = run_monte_carlo(cfg, run_dir, ep)
            ep_summary["monte_carlo"] = aggregate_monte_carlo(cfg, run_dir, ep, mc_results)
        end

        summary["endpoints"][epname] = ep_summary
    end

    if get(cfg["monte_carlo"], "enabled", false)
        _post_mc_diagnostics!(cfg, run_dir, summary)
    end

    try
        @info "==> Generating quantitative reports"
        summary["reports"] = generate_report(run_dir)
    catch err
        @warn "  [reports] failed" err
    end

    try
        summary["metadata"] = save_run_metadata(cfg, run_dir, summary)
    catch err
        @warn "  [meta] failed" err
    end

    @info "==> Done. Outputs under: $run_dir"
    return summary
end

# Convergence + pass-frequency + Morris sensitivity per endpoint.
# Failures are logged but do not abort the run.
function _post_mc_diagnostics!(cfg::AbstractDict, run_dir::AbstractString,
                                summary::AbstractDict)
    conv_out = Dict{String,Any}()
    morris_out = Dict{String,Any}()
    sens_cfg = get(cfg["monte_carlo"], "sensitivity", Dict{String,Any}())
    morris_enabled = Bool(get(sens_cfg, "morris", true))
    morris_traj = Int(get(sens_cfg, "morris_trajectories", 10))
    morris_levels = Int(get(sens_cfg, "morris_levels", 4))

     for ep in cfg["endpoints"]
        ep_mc_dir = joinpath(run_dir, "monte_carlo", ep["name"])
    ep_agg_dir = joinpath(run_dir, "monte_carlo", "aggregated", ep["name"])
        if isdir(ep_mc_dir)
            for analysis in ("lcc", "circuitscape")
                try
                    conv_out["$(ep["name"])__$analysis"] = check_convergence(ep_mc_dir;
                        analysis=analysis,
                        out_dir=ep_agg_dir)
                catch err
                    conv_out["$(ep["name"])__$analysis"] = Dict(
                        "converged" => nothing, "reason" => sprint(showerror, err))
                end
            end
        end  
        if morris_enabled
            try
                morris_out[ep["name"]] = compute_sensitivity_morris(cfg, run_dir, ep;
                    n_trajectories=morris_traj, n_levels=morris_levels)
            catch err
                morris_out[ep["name"]] = nothing
                @warn "  [morris] $(ep["name"]) failed" err
            end
        end
    end
    try
        buf = Float64(get(get(cfg, "analysis", Dict()), "pass_buffer_m", 2000.0))
        summary["pass_frequencies"] = compute_pass_frequencies_all_endpoints(cfg, run_dir;
                                                                              buffer_m=buf)
    catch err
        @warn "  [pass_freq] failed" err
    end
    summary["convergence"] = conv_out
    summary["sensitivity_morris"] = morris_out
    return summary
end
