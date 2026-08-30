# Montescape.jl
Master’s Thesis Project: MonteScape.jl - A Multi-Method Computational Framework for Spatial Movement Analysis with Sensitivity Analysis: The Case of Hannibal’s Alpine Crossing (218 BCE)

A pure-Julia framework for archaeological movement modelling.
Integrates **Least-Cost Path** (LCP), **Least-Cost Corridor** (LCC), and **Circuitscape** (circuit-theoretic flow) in a single configurable pipeline with native **Monte Carlo** uncertainty quantification.

Every Monte Carlo iteration rebuilds the friction surface from freshly sampled weights and (optionally) a perturbed DEM, then re-runs LCP / LCC / Circuitscape on that fresh surface. No silently re-baked slopes.

## 

- You **build your YAML config from scratch** in a browser-based form. No need  to learn the YAML schema or copy a pre-made config.
- Three independent actions in the always-visible footer:
  - **Save YAML** ; write the config to a file.
  - **Run locally** ; launch the pipeline in a separate terminal, on this machine, against the saved config.
  - **Generate HPC package** ; produce a self-contained folder with the config, your input rasters/vectors, the source code, a SLURM script and an
    HPC README. You copy that folder to your cluster and `sbatch` it.
- A validation banner lists every problem (missing Start point, weights summing to 0.85 instead of 1.0, DEM file that doesn't exist, etc.) before you ever launch a run.

## Install (no Julia setup required for end users)

Download the installer for your OS from the project Releases page and run it. The installer ships a private Julia runtime plus a precompiled sysimage, so the user never installs Julia themselves.

- **Windows**: `Montescape-Setup.exe`
- **macOS**: `Montescape-macOS.dmg`
- **Linux**: `Montescape-x86_64.AppImage`

Launch the app. A terminal window appears (the Julia HTTP server) and your default browser opens at `http://localhost:8765` showing a 7-tab form. Fill it in, click *Run locally* or *Generate HPC package*.

## Run from source (developers)

```bash
git clone <repo-url> movement_framework
cd movement_framework
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Launch the GUI from source (opens browser automatically)
julia --project=. app/desktop_web.jl

# Or run a config straight from the REPL
julia --project=. -e 'using MovementFramework; run_pipeline("path/to/config.yml")'

# Build the standalone binary + sysimage (15-30 min)
julia --project=. build/build_app.jl

# Run tests
julia --project=. -e 'using Pkg; Pkg.test()'
```

Requires Julia >= 1.10.

## What the form lets you configure

Every field is paired with the YAML key it writes (shown as a small caption underneath), so the form doubles as documentation. The footer shows a live YAML preview you can copy.

1. **Data**; DEM raster, sites vector, sites name column, optional via-points vector, grid resolution. Plain text paths; no file pickers.
2. **Routes** ; Start / End drop-downs auto-populated from the sites file via ArchGDAL; via-point multi-select from the passes file; per-route buffer. A *Quick build* helper turns one Start + one End + N picked
   via-points into N parallel routes.
3. **Friction**; weighted layers: slope (nine cost functions; Pandolf physiology block appears only when slope.method = pandolf), TRI, aspect, vegetation with reclass, hydrology with crossing penalty. Plus a custom
   layer block: list `.jl` files that call `register_layer!`, then add a row per layer with name / enabled / weight / description. Live red/green weight-sum indicator in the footer until weights sum to 1.0.
4. **Analyses**; LCP (enabled, knight moves), LCC (enabled, percent or absolute threshold, knight moves, export-as-polygon), Circuitscape (enabled, scenario, Julia threads, 4-neighbour mode, average resistances,
   unit currents, direct grounds, current maps, voltage maps, log transform).
5. **Monte Carlo** ; enabled, runs, workers, seed, sampling method (Dirichlet with concentration, Latin Hypercube with spread, or random), analyses subset, DEM-error profile (FABDEM / SRTM / lidar... pre-fills
   the published vertical sigma and correlation length), Morris sensitivity (trajectories, levels), aggregation outputs, pass-crossing buffer.
6. **HPC**; user / host / CPUs / memory / walltime / partition / account / scratch base. Feeds *Generate HPC package* only; never submits.
7. **Synthesis**; post-processing of completed runs. Auto-discovers finished runs under a parent folder; cross-run pass-frequency matrix; or pairwise B-A difference report with auto-discovered shared endpoints.

## Adding a custom friction layer

Write a Julia file that calls `register_layer!`:

```julia
# /abs/path/to/my_layer.jl
using MovementFramework: register_layer!

function compute_my_layer(dem, res, cfg)
    # ... your logic, returning a Float32 matrix clipped to [0, 2] ...
end
```

## Associated Master's Thesis
MonteScape was developed as part of my Master's thesis. The thesis will be added to the GitHub after having received my grade.
