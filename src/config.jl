# Config loading + validation.

using YAML
using Dates

"""
    load_config(path::AbstractString) -> Dict{String,Any}

Read a YAML config and resolve relative data paths against the repo root
(the directory containing Project.toml, walking up from the config file).
The project root is stored in cfg["_project_root"].
"""
function load_config(path::AbstractString)
    cfg_path = abspath(path)
    cfg = YAML.load_file(cfg_path; dicttype=Dict{String,Any})

    base = _find_project_root(cfg_path)
    cfg["_project_root"] = base
    cfg["_config_path"]  = cfg_path

    # Resolve data paths
    for key in ("dem_raster", "sites_vector", "passes_vector")
        if haskey(cfg["data"], key) && cfg["data"][key] !== nothing
            p = String(cfg["data"][key])
            if !isabspath(p)
                cfg["data"][key] = abspath(joinpath(base, p))
            end
        end
    end
    # Resolve custom layer file paths
    if haskey(cfg["friction"], "custom_layer_files")
        cfg["friction"]["custom_layer_files"] = [
            isabspath(String(p)) ?
                String(p) :
                abspath(joinpath(base, String(p)))
            for p in cfg["friction"]["custom_layer_files"]
        ]
    end
    # Expand $USER / ${USER} / %USERNAME% in output_dir, and fall back to
    # a writable local default when the configured output_dir does not exist
    # and cannot be created (typical on Windows when the YAML targets /scratch).
    if haskey(cfg["project"], "output_dir")
        cfg["project"]["output_dir"] = _expand_and_fallback_output_dir(
            String(cfg["project"]["output_dir"]))
    end

    validate_config(cfg)
    return cfg
end

# Walk up from the config file until we find Project.toml. If none is found
# (e.g. running a stand-alone packaged config), fall back to the config's parent.
function _find_project_root(cfg_path::AbstractString)
    d = dirname(abspath(cfg_path))
    for _ in 1:8
        isfile(joinpath(d, "Project.toml")) && return d
        parent = dirname(d)
        parent == d && break
        d = parent
    end
    return dirname(abspath(cfg_path))
end

"""
    _expand_user_path(p) -> String

Expand \$USER, \${USER}, and %USERNAME% on any platform using the current
user's name (taken from the USER or USERNAME environment variables, with
a safe fallback to "user" if neither is set).
"""
function _expand_user_path(p::AbstractString)
    user = get(ENV, "USER", get(ENV, "USERNAME", "user"))
    out = String(p)
    out = replace(out, "\$USER" => user)
    out = replace(out, "\${USER}" => user)
    out = replace(out, "%USERNAME%" => user)
    return out
end

# If the configured output dir lives under /scratch (or /work, /lustre) but
# we are not on a system where that path exists or can be created, fall back
# to ~/movement_outputs so local sanity-check runs do not blow up. The
# original cluster-targeting path is preserved for run_package generation,
# which rewrites it back to the canonical form before emitting the package.
function _expand_and_fallback_output_dir(p::AbstractString)
    expanded = _expand_user_path(p)
    if Sys.iswindows() &&
       (startswith(expanded, "/scratch") || startswith(expanded, "/work") || startswith(expanded, "/lustre"))
        fallback = joinpath(homedir(), "movement_outputs")
        @warn "Configured output_dir '$p' is an HPC path; using local fallback '$fallback'. The original path will be restored when generating a run package for HPC submission."
        return fallback
    end
    return expanded
end

"""
    validate_config(cfg)

Light schema check + friction weight sum verification. Throws ArgumentError
on missing required keys or inconsistent weights.
"""
function validate_config(cfg::AbstractDict)
    for k in ("project", "data", "grid", "endpoints", "friction")
        haskey(cfg, k) || throw(ArgumentError("Missing required config section: $k"))
    end

    fric = cfg["friction"]
    # `mode` is accepted for backward compatibility with old configs but only
    # the "build" path is supported. Prebuilt friction surfaces were dropped
    # because they can't participate in Monte Carlo weight perturbation.
    mode = get(fric, "mode", "build")
    if mode != "build"
        throw(ArgumentError(
            "friction.mode='$mode' is no longer supported. Only 'build' is " *
            "valid; the friction surface is always assembled from weighted layers " *
            "so Monte Carlo can perturb the weights."))
    end
    total = 0.0
    for key in ("slope", "tri", "aspect", "vegetation", "hydrology")
        sec = get(fric, key, nothing)
        sec === nothing && continue
        if get(sec, "enabled", false)
            total += Float64(get(sec, "weight", 0.0))
        end
    end
    for layer in get(fric, "custom_layers", Any[])
        if get(layer, "enabled", false)
            total += Float64(get(layer, "weight", 0.0))
        end
    end
    if !(abs(total - 1.0) < 0.01)
        throw(ArgumentError(
            "Friction weights must sum to 1.0 (got $(round(total, digits=4))). Adjust weights in config."))
    end
    return cfg
end

"""
    create_run_dir(cfg) -> String

Create the timestamped run directory and write a copy of the config inside it.
"""
function create_run_dir(cfg::AbstractDict)
    base   = cfg["_project_root"]
    outdir = cfg["project"]["output_dir"]
    outdir = isabspath(outdir) ? outdir : abspath(joinpath(base, outdir))
    mkpath(outdir)

    runid = get(cfg["project"], "run_id", "auto")
    if runid == "auto"
        ts = Dates.format(now(), "yyyymmdd_HHMMSS")
        runid = string(cfg["project"]["name"], "_", ts)
    end
    run_dir = joinpath(outdir, runid)
    mkpath(run_dir)

    YAML.write_file(joinpath(run_dir, "config_used.yml"), cfg)
    return run_dir
end
