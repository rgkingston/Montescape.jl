"""
    MovementFramework

A pure-Julia framework for archaeological movement modelling.
Integrates LCP, LCC and Circuitscape with native Monte Carlo uncertainty
quantification.

Main entry point:

    run_pipeline("path/to/config.yml")
"""
module MovementFramework

using Logging

# --- Core types & config ---------------------------------------------------
include("config.jl")
include("io.jl")

# --- Friction & layers -----------------------------------------------------
include("friction/registry.jl")
include("friction/slope.jl")
include("friction/tri.jl")
include("friction/aspect.jl")
include("friction/vegetation.jl")
include("friction/hydrology.jl")
include("friction/friction.jl")

# --- Analyses --------------------------------------------------------------
include("lcp.jl")
include("lcc.jl")
include("circuitscape.jl")

# --- Monte Carlo + aggregation --------------------------------------------
include("dem_error.jl")
include("monte_carlo.jl")
include("aggregation.jl")

# --- Post-MC diagnostics + comparison + sensitivity -----------------------
include("convergence.jl")
include("pass_frequency.jl")
include("sensitivity.jl")
include("compare_runs.jl")
include("compare_scenarios.jl")
include("visualize.jl")
include("metadata.jl")

# --- HPC + pipeline orchestrator ------------------------------------------
include("hpc.jl")
include("pipeline.jl")
include("run_package.jl")
include("cli.jl")
include("runner.jl")


# --- Web GUI (submodule, also serves as PackageCompiler entry) ------------
include("web_gui.jl")

"""
    julia_main() :: Cint

Standalone-app entry point used by PackageCompiler (`build/build_app.jl`).
Launches the web GUI; the browser is the front-end, everything else runs
here. Returns 0 on clean shutdown, 1 on any uncaught error.
"""
function julia_main()::Cint
    try
        WebGUI.main()
    catch err
        @error "MovementFramework crashed" exception=(err, catch_backtrace())
        return Cint(1)
    end
    return Cint(0)
end

# --- Public API ------------------------------------------------------------
export run_pipeline,
       load_config, validate_config, create_run_dir,
       build_friction_surface,
       run_lcp, run_lcc, run_circuitscape,
       run_monte_carlo, aggregate_monte_carlo,
       register_layer!, list_registered_layers,
       compute_slope_friction, compute_tri, compute_aspect_penalty,
       perturb_dem,
       check_convergence,
       compute_pass_frequencies_for_endpoint,
       compute_pass_frequencies_all_endpoints,
       compute_sensitivity_morris,
       build_pass_frequency_matrix,
       compare_scenarios,
       generate_report,
       save_run_metadata,
       write_slurm_script,
       generate_run_package,
       spawn_pipeline_in_terminal,
       read_raster, write_raster, read_points, find_point,
       world_to_pixel, pixel_to_world, pixel_resolution

end # module
