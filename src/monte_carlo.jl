# Monte Carlo orchestrator.
# Each iteration:
#   1. Sample new weight vector from the requested distribution.
#   2. Optionally perturb the DEM with spatially-correlated Gaussian noise.
#   3. REBUILD the friction surface from those weights + perturbed DEM.
#   4. Re-run requested analyses on that fresh friction surface.
#
# Iterations are independent and run in parallel via Distributed.pmap.
# Workers are spawned lazily on the first call to run_monte_carlo so a
# laptop run with mc_workers=1 stays single-process, while an HPC job
# with mc_workers=32 fans out across cores.

using Random
using Distributions
using Distributed
using ProgressMeter
using Statistics

"""
    run_monte_carlo(cfg, run_dir, endpoint) -> Dict

Run the Monte Carlo ensemble for one endpoint. Results are written to
`<run_dir>/monte_carlo/<endpoint_name>/iter_<i>/...`. Aggregation is
performed afterwards by `aggregate_monte_carlo`.

Parallelism: iterations are dispatched with `pmap` over a worker pool
of size `monte_carlo.n_workers`. Workers are added on demand; each
worker `using MovementFramework` once so per-iteration overhead is
just the function call + serialised arguments.
"""
function run_monte_carlo(cfg::AbstractDict, run_dir::AbstractString, endpoint::AbstractDict)
    mc = cfg["monte_carlo"]
    get(mc, "enabled", false) || return Dict()

    n      = Int(get(mc, "n_runs", 100))
    seed   = Int(get(mc, "seed", 42))
    ana    = collect(String.(get(mc, "analyses", ["lcp"])))
    nw     = max(1, Int(parse(Int, get(ENV, "MONTESCAPE_MC_WORKERS", string(Int(get(mc, "n_workers", 1)))))))
    epname = endpoint["name"]
    ep_dir = joinpath(run_dir, "monte_carlo", epname)
    mkpath(ep_dir)

    # Load fixed inputs once on the driver. They get serialised to
    # workers as part of the pmap closure (Julia handles transport).
    points = read_points(cfg["data"]["sites_vector"];
                         name_field=get(cfg["data"], "site_name_field", "name"))
    src = find_point(points, endpoint["source_name"])
    gnd = find_point(points, endpoint["ground_name"])

    waypoints = Tuple{Float64,Float64}[]
    if !isempty(get(endpoint, "waypoints", String[])) &&
       get(cfg["data"], "passes_vector", nothing) !== nothing
        passes = read_points(cfg["data"]["passes_vector"]; name_field="name")
        for wn in endpoint["waypoints"]
            p = find_point(passes, wn)
            push!(waypoints, (p.x, p.y))
        end
    end

    base_weights = _collect_base_weights(cfg)
    sampler = _build_sampler(mc, base_weights, seed)

    @info "  [MC] $epname: $n iterations, $nw worker(s) (yaml=$(get(mc, "n_workers", 1)), env=$(get(ENV, "MONTESCAPE_MC_WORKERS", "not set"))), analyses=$(ana), sampler=$(get(mc, "sampling", "dirichlet"))"
    # Bring workers online and load the package on each, but only when
    # we have something to parallelise across.
    if nw > 1
        _ensure_workers(nw, cfg)
    end

    function _run_one(i)
        iter_dir = joinpath(ep_dir, "iter_$(lpad(i, 4, '0'))")
        mkpath(iter_dir)
        rng = MersenneTwister(seed + i)

        weights = sampler(rng)
        dem_ref = read_raster(cfg["data"]["dem_raster"])
        dem_perturbed = perturb_dem(dem_ref.data, pixel_resolution(dem_ref),
                                    get(mc, "dem_error", Dict()), rng)

        _, friction, ref = build_friction_surface(cfg, iter_dir;
                                                   dem_override=dem_perturbed,
                                                   weight_override=weights,
                                                   write=true)

        iter_result = Dict{String,Any}("weights" => weights)
        if "lcp" in ana
            iter_result["lcp"] = run_lcp(friction, ref, (src.x, src.y), (gnd.x, gnd.y),
                                          iter_dir; cfg=cfg, name="iter",
                                          waypoints=waypoints,
                                          dem=dem_perturbed)
        end
        if "lcc" in ana
            iter_result["lcc"] = run_lcc(friction, ref, (src.x, src.y), (gnd.x, gnd.y),
                                          iter_dir; cfg=cfg, name="iter",
                                          dem=dem_perturbed)
        end
        if "circuitscape" in ana
            iter_result["circuitscape"] = run_circuitscape(friction, ref,
                (src.x, src.y), (gnd.x, gnd.y), iter_dir; cfg=cfg, name="iter")
        end
        return i => iter_result
    end

    progress = Progress(n; desc="  [MC] $epname ")
    pairs = if nw > 1
        wp = WorkerPool(workers())
        out = pmap(wp, 1:n) do i
            _run_one(i)
        end
        # Cross-process progress would need ParallelProgressMeter; this is
        # simpler and good enough (bar fills in one burst at the end).
        for _ in 1:n; next!(progress); end
        out
    else
        out = Vector{Pair{Int,Dict{String,Any}}}(undef, n)
        for i in 1:n
            out[i] = _run_one(i)
            next!(progress)
        end
        out
    end

    return Dict{Int,Dict{String,Any}}(pairs)
end

# Add workers up to the requested count and make sure each one has the
# package loaded. Idempotent: extra calls with the same nw are a no-op.
function _ensure_workers(nw::Integer, cfg::AbstractDict)
    needed = nw - nworkers()
    if needed > 0
        @info "  [MC] spawning $needed extra worker process(es) (target total: $nw)"
        addprocs(needed; exeflags="--project=$(dirname(Base.active_project()))")  # ? fix
    end
    # Load MovementFramework + ArchGDAL on all workers
    @everywhere begin
        try
            using MovementFramework
            using ArchGDAL
        catch err
            @warn "worker failed to load MovementFramework" err
        end
    end
    # Load custom layer files so _LAYER_REGISTRY gets populated on each worker
    clf_paths = get(get(cfg, "friction", Dict()), "custom_layer_files", String[])
    if !isempty(clf_paths)
        @everywhere MovementFramework.load_custom_layer_files!($clf_paths)
    end
    return nworkers()
end
# Collect the base weights of all enabled friction layers (slope, tri,
# aspect, vegetation, hydrology, plus every enabled custom layer).
function _collect_base_weights(cfg::AbstractDict)
    fric = cfg["friction"]
    bw = Dict{String,Float64}()
    for key in ("slope", "tri", "aspect", "vegetation", "hydrology")
        sec = get(fric, key, nothing)
        sec === nothing && continue
        if get(sec, "enabled", false)
            bw[key] = Float64(get(sec, "weight", 0.0))
        end
    end
    for layer in get(fric, "custom_layers", Any[])
        if get(layer, "enabled", false)
            bw[String(layer["name"])] = Float64(get(layer, "weight", 0.0))
        end
    end
    return bw
end

function _build_sampler(mc::AbstractDict, base::Dict{String,Float64}, seed::Int)
    method = get(mc, "sampling", "dirichlet")
    names = collect(keys(base))
    base_vec = [base[n] for n in names]

    if method == "dirichlet"
        alpha = Float64(get(mc, "dirichlet_concentration", 6.0)) .* base_vec
        dist = Dirichlet(alpha)
        return rng -> Dict(zip(names, rand(rng, dist)))
    elseif method == "latin_hypercube"
        spread = Float64(get(mc, "lhs_spread", 0.5))
        return function (rng)
            v = [base[n] * (1.0 + spread * (2 * rand(rng) - 1)) for n in names]
            v = max.(v, 0.0)
            s = sum(v)
            s == 0 && return Dict(zip(names, base_vec))
            return Dict(zip(names, v ./ s))
        end
    elseif method == "random"
        return function (rng)
            v = rand(rng, length(names))
            return Dict(zip(names, v ./ sum(v)))
        end
    else
        throw(ArgumentError("Unknown sampling method: $method (use dirichlet, latin_hypercube, random)"))
    end
end
