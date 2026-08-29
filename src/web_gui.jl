# Web GUI launcher for MovementFramework, wrapped as a submodule so it can
# be included from MovementFramework.jl and used as the julia_main entry
# point for PackageCompiler standalone builds.
#
# Serves a single HTML page over HTTP.jl on http://localhost:8765 and opens
# the user's default browser automatically. The browser is just a renderer;
# all logic runs in Julia.
#
# Design principles:
#   - The user BUILDS their own YAML from a blank form. No premade configs.
#   - Every field that affects a run is editable.
#   - File paths are plain text inputs (no native file pickers).
#   - Three independent actions: Save YAML / Run locally / Generate HPC package.

module WebGUI

using HTTP
using JSON
using YAML
using Sockets
using Dates
using ..MovementFramework

const PORT = 8765

# Resolve the HTML file. In dev (running from the repo) it sits at
# <repo>/app/web/index.html. @__DIR__ here is <repo>/src, so dirname(@__DIR__)
# is the repo root. In a PackageCompiler bundle the source tree is not
# shipped, so we also look next to the binary.
function _html_file()
    candidates = String[
        joinpath(dirname(@__DIR__), "app", "web", "index.html"),
        joinpath(Sys.BINDIR, "..", "share", "MovementFramework", "web", "index.html"),
        joinpath(Sys.BINDIR, "..", "web", "index.html"),
    ]
    for c in candidates
        isfile(c) && return c
    end
    return candidates[1]  # let the read error surface meaningfully
end

# ============================================================================
# Default state
# ============================================================================

function default_state()
    user = get(ENV, "USER", get(ENV, "USERNAME", "user"))
    home_outputs = joinpath(homedir(), "montescape_runs")
    Dict{String,Any}(
        # project
        "project_name"        => "my_route_study",
        "output_dir"          => home_outputs,
        "seed"                => 42,
        "load_yaml_path"      => "",
        # data + grid
        "resolution_m"        => 90,
        "dem_raster"          => "",
        "sites_vector"        => "",
        "passes_vector"       => "",
        "site_name_field"     => "name",
        # friction high-level
        "friction_output_name"=> "friction_surface.tif",
        # NOTE: custom_layer_files is derived from friction.custom_layers[*].file_path
        # at save time, so the UI only stores it on the rows themselves.
        # LCP
        "lcp_enabled"         => true,
        "lcp_knight"          => true,
        # LCC
        "lcc_enabled"         => true,
        "lcc_threshold_method"=> "percent",
        "lcc_threshold_value" => 5.0,
        "lcc_knight"          => true,
        "lcc_export_vector"   => true,
        # Circuitscape
        "cs_enabled"          => true,
        "cs_scenario"         => "advanced",
        "cs_julia_threads"    => 4,
        "cs_four_neighbours"  => false,
        "cs_avg_resistances"  => true,
        "cs_unit_currents"    => true,
        "cs_direct_grounds"   => true,
        "cs_write_cur_maps"   => true,
        "cs_write_volt_maps"  => false,
        "cs_log_transform"    => false,
        # Monte Carlo
        "mc_enabled"          => true,
        "mc_runs"             => 100,
        "mc_workers"          => 8,
        "mc_sampling"         => "dirichlet",
        "mc_seed"             => 42,
        "mc_dirichlet_conc"   => 6.0,
        "mc_lhs_spread"       => 0.5,
        "mc_analyses"         => ["lcp", "lcc", "circuitscape"],
        "pass_buffer_m"       => 2000,
        # DEM error
        "dem_err_enabled"     => true,
        "dem_profile"         => "",
        "dem_err_m"           => "",
        "dem_corr_len_m"      => "",
        "dem_method"          => "smoothed_gaussian",
        # Morris sensitivity
        "morris_enabled"      => true,
        "morris_trajectories" => 15,
        "morris_levels"       => 4,
        # Aggregation
        "agg_probability"     => true,
        "agg_mean_current"    => true,
        "agg_confidence"      => 95,
        # HPC (independent of local mc_workers / cs_julia_threads)
        "hpc_user"            => user,
        "hpc_host"            => "hpc.cluster.example",
        "cpus"                => 96,
        "mem_gb"              => 600,
        "walltime"            => "24:00:00",
        "partition"           => "",
        "account"             => "",
        "scratch_base"        => "/scratch/\$USER",
        # Synthesis
        "synthesis_parent"    => home_outputs,
        "synthesis_out_dir"   => joinpath(home_outputs, "_synthesis"),
        "scenario_a"          => "",
        "scenario_b"          => "",
        "scenario_endpoint"   => "",
        # Nested
        "endpoints"           => Any[Dict{String,Any}(
            "name"=>"my_route", "source_name"=>"", "ground_name"=>"",
            "buffer_m"=>1000, "waypoints"=>String[])],
        "friction"            => _default_friction(),
    )
end

function _default_friction()
    Dict{String,Any}(
        "slope" => Dict{String,Any}(
            "enabled"=>true, "weight"=>0.30, "method"=>"pandolf",
            "body_weight_kg"=>80.0, "load_kg"=>25.0,
            "terrain_factor"=>1.0, "speed_ms"=>1.3,
        ),
        "tri"        => Dict{String,Any}("enabled"=>true, "weight"=>0.05),
        "aspect"     => Dict{String,Any}("enabled"=>true, "weight"=>0.05, "preferred_direction"=>180),
        "vegetation" => Dict{String,Any}("enabled"=>false, "weight"=>0.0, "raster"=>"", "reclass"=>Dict{String,Any}()),
        "hydrology"  => Dict{String,Any}("enabled"=>false, "weight"=>0.0, "raster"=>"", "crossing_penalty"=>8),
        "custom_layers" => Any[],
    )
end

const STATE = Ref{Dict{String,Any}}(Dict{String,Any}())

# ============================================================================
# YAML (flat state <-> nested YAML)
# ============================================================================

function build_yaml_dict(s::AbstractDict)
    fric_out = deepcopy(s["friction"])
    fric_out["mode"] = "build"
    fric_out["output_name"] = s["friction_output_name"]
    layer_files = String[]
    cleaned_layers = Any[]
    for layer in get(fric_out, "custom_layers", Any[])
        fp = String(get(layer, "file_path", ""))
        if !isempty(fp) && !(fp in layer_files)
            push!(layer_files, fp)
        end
        l = Dict{String,Any}()
        for (k, v) in layer
            k == "file_path" && continue
            l[k] = v
        end
        push!(cleaned_layers, l)
    end
    fric_out["custom_layer_files"] = layer_files
    fric_out["custom_layers"] = cleaned_layers

    return Dict{String,Any}(
        "project" => Dict(
            "name" => s["project_name"], "output_dir" => s["output_dir"],
            "run_id" => "auto", "seed" => s["seed"]),
        "analysis" => Dict("pass_buffer_m" => Float64(s["pass_buffer_m"])),
        "data" => Dict(
            "dem_raster" => s["dem_raster"],
            "sites_vector" => s["sites_vector"],
            "site_name_field" => s["site_name_field"],
            "passes_vector" => isempty(s["passes_vector"]) ? nothing : s["passes_vector"]),
        "grid" => Dict("resolution_m" => s["resolution_m"], "snap_to" => "dem"),
        "endpoints" => deepcopy(s["endpoints"]),
        "friction" => fric_out,
        "lcp" => Dict("enabled" => s["lcp_enabled"], "knight" => s["lcp_knight"]),
        "lcc" => Dict(
            "enabled" => s["lcc_enabled"],
            "threshold" => Dict(
                "method" => s["lcc_threshold_method"],
                "value"  => Float64(s["lcc_threshold_value"])),
            "knight" => s["lcc_knight"],
            "export_vector" => s["lcc_export_vector"]),
        "circuitscape" => Dict(
            "enabled" => s["cs_enabled"],
            "scenario" => s["cs_scenario"],
            "julia_threads" => Int(s["cs_julia_threads"]),
            "connect_four_neighbors_only" => s["cs_four_neighbours"],
            "use_average_resistances"     => s["cs_avg_resistances"],
            "use_unit_currents"           => s["cs_unit_currents"],
            "use_direct_grounds"          => s["cs_direct_grounds"],
            "write_cur_maps"              => s["cs_write_cur_maps"],
            "write_volt_maps"             => s["cs_write_volt_maps"],
            "log_transform_maps"          => s["cs_log_transform"]),
        "monte_carlo" => let mc = Dict{String,Any}(
                "enabled"    => s["mc_enabled"],
                "n_runs"     => Int(s["mc_runs"]),
                "n_workers"  => Int(s["mc_workers"]),
                "sampling"   => s["mc_sampling"],
                "seed"       => Int(s["mc_seed"]),
                "dem_error"  => Dict(
                    "enabled"              => s["dem_err_enabled"],
                    "dem_profile"          => s["dem_profile"] == "" ? nothing : s["dem_profile"],
                    "error_m"              => s["dem_err_m"] == "" ? 0.0 : Float64(s["dem_err_m"]),
                    "correlation_length_m" => s["dem_corr_len_m"] == "" ? 0.0 : Float64(s["dem_corr_len_m"]),
                    "method"               => s["dem_method"]),
                "analyses"   => collect(String, s["mc_analyses"]),
                "sensitivity" => Dict(
                    "morris"               => s["morris_enabled"],
                    "morris_trajectories"  => Int(s["morris_trajectories"]),
                    "morris_levels"        => Int(s["morris_levels"])),
                "aggregate"  => Dict(
                    "probability_surface"  => s["agg_probability"],
                    "mean_current"         => s["agg_mean_current"],
                    "confidence_interval"  => s["agg_confidence"] / 100.0),
            )
            if s["mc_sampling"] == "dirichlet"
                mc["dirichlet_concentration"] = Float64(s["mc_dirichlet_conc"])
            elseif s["mc_sampling"] == "latin_hypercube"
                mc["lhs_spread"] = Float64(s["mc_lhs_spread"])
            end
            mc
        end,
    )
end

function apply_yaml!(s::AbstractDict, cfg::AbstractDict)
    proj = get(cfg, "project", Dict()); data = get(cfg, "data", Dict())
    s["project_name"]    = String(get(proj, "name", "my_route_study"))
    s["output_dir"]      = String(get(proj, "output_dir", s["output_dir"]))
    s["seed"]            = Int(get(proj, "seed", 42))
    s["dem_raster"]      = String(get(data, "dem_raster", ""))
    s["sites_vector"]    = String(get(data, "sites_vector", ""))
    pv = get(data, "passes_vector", "")
    s["passes_vector"]   = pv === nothing ? "" : String(pv)
    s["site_name_field"] = String(get(data, "site_name_field", "name"))
    s["resolution_m"]    = Int(get(get(cfg, "grid", Dict()), "resolution_m", 90))
    eps = collect(Any, get(cfg, "endpoints", Any[]))
    s["endpoints"]       = isempty(eps) ? Any[Dict{String,Any}(
        "name"=>"my_route", "source_name"=>"", "ground_name"=>"",
        "buffer_m"=>1000, "waypoints"=>String[])] : eps

    incoming = get(cfg, "friction", Dict{String,Any}())
    merged = deepcopy(_default_friction())
    for (k, v) in incoming
        if k in ("slope", "tri", "aspect", "vegetation", "hydrology")
            merge!(merged[k], Dict{String,Any}(v))
        elseif k == "custom_layers"
            merged["custom_layers"] = collect(Any, v)
        end
    end
    yaml_files = collect(String, get(incoming, "custom_layer_files", String[]))
    yaml_layers = collect(Any, get(incoming, "custom_layers", Any[]))
    stitched = Any[]
    for (i, layer) in enumerate(yaml_layers)
        d = Dict{String,Any}(layer)
        d["file_path"] = i <= length(yaml_files) ? yaml_files[i] : ""
        push!(stitched, d)
    end
    merged["custom_layers"] = stitched
    s["friction"] = merged
    s["friction_output_name"] = String(get(incoming, "output_name", "friction_surface.tif"))

    lcp = get(cfg, "lcp", Dict()); lcc = get(cfg, "lcc", Dict()); cs = get(cfg, "circuitscape", Dict())
    s["lcp_enabled"] = Bool(get(lcp, "enabled", true)); s["lcp_knight"] = Bool(get(lcp, "knight", true))
    s["lcc_enabled"] = Bool(get(lcc, "enabled", true))
    thr = get(lcc, "threshold", Dict())
    s["lcc_threshold_method"] = String(get(thr, "method", "percent"))
    s["lcc_threshold_value"]  = Float64(get(thr, "value", 5.0))
    s["lcc_knight"]           = Bool(get(lcc, "knight", true))
    s["lcc_export_vector"]    = Bool(get(lcc, "export_vector", true))

    s["cs_enabled"]          = Bool(get(cs, "enabled", true))
    s["cs_scenario"]         = String(get(cs, "scenario", "advanced"))
    s["cs_julia_threads"]    = Int(get(cs, "julia_threads", 4))
    s["cs_four_neighbours"]  = Bool(get(cs, "connect_four_neighbors_only", false))
    s["cs_avg_resistances"]  = Bool(get(cs, "use_average_resistances", true))
    s["cs_unit_currents"]    = Bool(get(cs, "use_unit_currents", true))
    s["cs_direct_grounds"]   = Bool(get(cs, "use_direct_grounds", true))
    s["cs_write_cur_maps"]   = Bool(get(cs, "write_cur_maps", true))
    s["cs_write_volt_maps"]  = Bool(get(cs, "write_volt_maps", false))
    s["cs_log_transform"]    = Bool(get(cs, "log_transform_maps", false))

    mc = get(cfg, "monte_carlo", Dict())
    s["mc_enabled"]      = Bool(get(mc, "enabled", true))
    s["mc_runs"]         = Int(get(mc, "n_runs", 100))
    s["mc_workers"]      = Int(get(mc, "n_workers", 8))
    s["mc_sampling"]     = String(get(mc, "sampling", "dirichlet"))
    s["mc_seed"]         = Int(get(mc, "seed", s["seed"]))
    s["mc_dirichlet_conc"] = Float64(get(mc, "dirichlet_concentration", 6.0))
    s["mc_lhs_spread"]   = Float64(get(mc, "lhs_spread", 0.5))
    s["mc_analyses"]     = collect(String, get(mc, "analyses", ["lcp", "lcc", "circuitscape"]))
    s["pass_buffer_m"]   = Int(round(Float64(get(get(cfg, "analysis", Dict()), "pass_buffer_m", 2000.0))))

    de = get(mc, "dem_error", Dict())
    s["dem_err_enabled"] = Bool(get(de, "enabled", true))
    dp = get(de, "dem_profile", "")
    s["dem_profile"]     = dp === nothing ? "" : String(dp)
    let v = get(de, "error_m", "")
        s["dem_err_m"] = (v === nothing || v == "") ? "" : Float64(v)
    end
    let v = get(de, "correlation_length_m", "")
        s["dem_corr_len_m"] = (v === nothing || v == "") ? "" : Int(round(Float64(v)))
    end
    s["dem_method"]      = String(get(de, "method", "smoothed_gaussian"))

    sens = get(mc, "sensitivity", Dict())
    s["morris_enabled"]      = Bool(get(sens, "morris", true))
    s["morris_trajectories"] = Int(get(sens, "morris_trajectories", 10))
    s["morris_levels"]       = Int(get(sens, "morris_levels", 4))

    agg = get(mc, "aggregate", Dict())
    s["agg_probability"]  = Bool(get(agg, "probability_surface", true))
    s["agg_mean_current"] = Bool(get(agg, "mean_current", true))
    s["agg_confidence"]   = Int(round(100 * Float64(get(agg, "confidence_interval", 0.95))))
    return s
end

function current_weight_sum(s::AbstractDict)
    fric = s["friction"]
    total = 0.0
    for k in ("slope", "tri", "aspect", "vegetation", "hydrology")
        sec = get(fric, k, nothing); sec === nothing && continue
        get(sec, "enabled", false) && (total += Float64(get(sec, "weight", 0.0)))
    end
    for layer in get(fric, "custom_layers", Any[])
        get(layer, "enabled", false) && (total += Float64(get(layer, "weight", 0.0)))
    end
    return total
end

# Plain-English validation. Empty list = ready to run.
function validate_state(s::AbstractDict)
    issues = String[]
    isempty(strip(String(s["project_name"]))) && push!(issues, "Run name is empty.")
    isempty(strip(String(s["output_dir"]))) && push!(issues, "Output directory is empty.")
    isempty(strip(String(s["dem_raster"]))) && push!(issues, "DEM raster path is empty.")
    isempty(strip(String(s["sites_vector"]))) && push!(issues, "Sites vector path is empty.")
    if !isempty(s["dem_raster"]) && !isfile(s["dem_raster"])
        push!(issues, "DEM raster file does not exist: $(s["dem_raster"])")
    end
    if !isempty(s["sites_vector"]) && !isfile(s["sites_vector"])
        push!(issues, "Sites vector file does not exist: $(s["sites_vector"])")
    end
    if !isempty(s["passes_vector"]) && !isfile(s["passes_vector"])
        push!(issues, "Via-points vector file does not exist: $(s["passes_vector"])")
    end
    eps = s["endpoints"]
    isempty(eps) && push!(issues, "No routes configured. Add at least one route on the Routes tab.")
    for (i, ep) in enumerate(eps)
        isempty(strip(String(get(ep, "name", "")))) &&
            push!(issues, "Route $i has no name.")
        isempty(strip(String(get(ep, "source_name", "")))) &&
            push!(issues, "Route $i has no Start point.")
        isempty(strip(String(get(ep, "ground_name", "")))) &&
            push!(issues, "Route $i has no End point.")
    end
    ws = current_weight_sum(s)
    if !(abs(ws - 1.0) < 0.01)
        push!(issues, "Friction layer weights sum to $(round(ws, digits=3)) but must sum to 1.00.")
    end
    # Custom layer rows: warn if the file is missing, if the name is blank, or
    # if the typed name doesn't match any register_layer! call in the file.
    for (i, layer) in enumerate(get(s["friction"], "custom_layers", Any[]))
        get(layer, "enabled", false) || continue
        fp = String(get(layer, "file_path", ""))
        nm = strip(String(get(layer, "name", "")))
        if isempty(fp)
            push!(issues, "Custom layer row $i has no Julia file path.")
            continue
        elseif !isfile(fp)
            push!(issues, "Custom layer row $i file does not exist: $fp")
            continue
        end
        if isempty(nm)
            push!(issues, "Custom layer row $i has no name (must match the string in register_layer!).")
            continue
        end
        available = _registered_names_in_file(fp)
        if !isempty(available) && !(nm in available)
            available_str = join(available, ", ")
            push!(issues,
                "Custom layer row $i: name '$nm' is not registered by $fp. " *
                "That file registers: $available_str.")
        end
    end
    return issues
end

# ============================================================================
# HTTP handlers
# ============================================================================

function _json_response(obj; status::Int=200)
    HTTP.Response(status, ["Content-Type" => "application/json"],
                   body=JSON.json(obj))
end

function handle_index(::HTTP.Request)
    body = read(_html_file(), String)
    return HTTP.Response(200, ["Content-Type" => "text/html; charset=utf-8"], body=body)
end

function handle_get_state(::HTTP.Request)
    return _json_response(Dict(
        "state"      => STATE[],
        "weight_sum" => current_weight_sum(STATE[]),
        "issues"     => validate_state(STATE[]),
    ))
end

function handle_set_state(req::HTTP.Request)
    try
        incoming = JSON.parse(String(req.body))
        for (k, v) in incoming
            STATE[][k] = v
        end
        return _json_response(Dict(
            "ok"         => true,
            "weight_sum" => current_weight_sum(STATE[]),
            "issues"     => validate_state(STATE[]),
        ))
    catch err
        return _json_response(Dict("ok" => false, "error" => sprint(showerror, err)); status=400)
    end
end

# Scan a Julia file for register_layer!("name", ...) calls and return the names.
function _registered_names_in_file(path::AbstractString)
    names = String[]
    isfile(path) || return names
    text = try; read(path, String) catch; "" end
    for m in eachmatch(r"register_layer!\s*\(\s*\"([^\"]+)\""s, text)
        push!(names, String(m.captures[1]))
    end
    return names
end

function handle_inspect_layer_file(req::HTTP.Request)
    try
        body = JSON.parse(String(req.body))
        path = String(get(body, "path", ""))
        if isempty(path)
            return _json_response(Dict("names" => String[], "exists" => false))
        end
        return _json_response(Dict(
            "names"  => _registered_names_in_file(path),
            "exists" => isfile(path)))
    catch err
        return _json_response(Dict("names" => String[], "error" => sprint(showerror, err)))
    end
end

function handle_list_vector_names(req::HTTP.Request)
    try
        body = JSON.parse(String(req.body))
        path  = String(get(body, "path", ""))
        field = String(get(body, "field", "name"))
        isempty(path) && return _json_response(Dict("names" => String[]))
        isfile(path) || return _json_response(Dict("names" => String[],
                                                     "error" => "file not found"))
        pts = MovementFramework.read_points(path; name_field=field)
        names = sort(unique(String[String(p.name) for p in pts]))
        return _json_response(Dict("names" => names))
    catch err
        return _json_response(Dict("names" => String[], "error" => sprint(showerror, err)))
    end
end

function handle_preview_yaml(::HTTP.Request)
    try
        cfg = build_yaml_dict(STATE[])
        io = IOBuffer()
        YAML.write(io, cfg)
        return _json_response(Dict("yaml" => String(take!(io))))
    catch err
        return _json_response(Dict("yaml" => "", "error" => sprint(showerror, err)))
    end
end

function handle_validate(::HTTP.Request)
    return _json_response(Dict(
        "issues"     => validate_state(STATE[]),
        "weight_sum" => current_weight_sum(STATE[]),
    ))
end

function handle_save_yaml(req::HTTP.Request)
    try
        body = JSON.parse(String(req.body))
        path = String(get(body, "path", ""))
        if isempty(path)
            s = STATE[]
            path = joinpath(s["output_dir"], string(s["project_name"], ".yml"))
        end
        cfg = build_yaml_dict(STATE[])
        mkpath(dirname(abspath(path)))
        YAML.write_file(path, cfg)
        return _json_response(Dict("ok" => true, "path" => abspath(path)))
    catch err
        return _json_response(Dict("ok" => false, "error" => sprint(showerror, err)); status=500)
    end
end

function handle_run_local(::HTTP.Request)
    try
        issues = validate_state(STATE[])
        isempty(issues) || return _json_response(Dict(
            "ok" => false,
            "error" => "Fix these problems first:\n  - " * join(issues, "\n  - ")); status=400)
        cfg = build_yaml_dict(STATE[])
        path = joinpath(tempdir(), "mf_$(Dates.format(now(), "yyyymmddHHMMSS")).yml")
        YAML.write_file(path, cfg)
        proj = dirname(Base.active_project())
        MovementFramework.spawn_pipeline_in_terminal(path; project_dir=proj)
        return _json_response(Dict(
            "ok" => true,
            "message" => "Pipeline launched in a new terminal",
            "config_path" => path))
    catch err
        return _json_response(Dict("ok" => false, "error" => sprint(showerror, err)); status=500)
    end
end

function handle_make_slurm(req::HTTP.Request)
    try
        issues = validate_state(STATE[])
        isempty(issues) || return _json_response(Dict(
            "ok" => false,
            "error" => "Fix these problems first:\n  - " * join(issues, "\n  - ")); status=400)
        body = try JSON.parse(String(req.body)) catch; Dict() end
        path = String(get(body, "yaml_path", ""))
        s = STATE[]
        if isempty(path)
            path = joinpath(s["output_dir"], string(s["project_name"], ".yml"))
            cfg = build_yaml_dict(s)
            mkpath(dirname(abspath(path)))
            YAML.write_file(path, cfg)
        end
        sbatch = MovementFramework.write_slurm_script(path;
            cpus = s["cpus"], mem_gb = s["mem_gb"], walltime = s["walltime"],
            partition = isempty(s["partition"]) ? nothing : s["partition"],
            account   = isempty(s["account"])   ? nothing : s["account"],
            project_path = dirname(Base.active_project()))
        return _json_response(Dict(
            "ok" => true,
            "yaml_path"   => abspath(path),
            "sbatch_path" => abspath(sbatch),
            "message"     => "Wrote SLURM script. From a cluster login node:  sbatch $(basename(sbatch))"))
    catch err
        return _json_response(Dict("ok" => false, "error" => sprint(showerror, err)); status=500)
    end
end

function handle_generate_package(::HTTP.Request)
    try
        issues = validate_state(STATE[])
        isempty(issues) || return _json_response(Dict(
            "ok" => false,
            "error" => "Fix these problems first:\n  - " * join(issues, "\n  - ")); status=400)
        s = STATE[]
        cfg = build_yaml_dict(s)
        tmp_cfg = joinpath(tempdir(), "mf_pkg_$(Dates.format(now(), "yyyymmddHHMMSS")).yml")
        YAML.write_file(tmp_cfg, cfg)
        target = MovementFramework.generate_run_package(tmp_cfg;
            output_dir = s["output_dir"],
            package_name = "$(s["project_name"])_runpkg",
            cpus = s["cpus"], mem_gb = s["mem_gb"], walltime = s["walltime"],
            partition = isempty(s["partition"]) ? nothing : s["partition"],
            account   = isempty(s["account"])   ? nothing : s["account"],
            hpc_user  = s["hpc_user"], hpc_host = s["hpc_host"],
            scratch_base = s["scratch_base"])
        return _json_response(Dict("ok" => true, "path" => target))
    catch err
        return _json_response(Dict("ok" => false, "error" => sprint(showerror, err)); status=500)
    end
end

function handle_load_yaml(req::HTTP.Request)
    try
        body = JSON.parse(String(req.body))
        path = String(get(body, "path", ""))
        isempty(path) && return _json_response(Dict("ok" => false,
            "error" => "no path provided"); status=400)
        isfile(path) || return _json_response(Dict("ok" => false,
            "error" => "file not found: $path"); status=400)
        cfg = MovementFramework.load_config(path)
        apply_yaml!(STATE[], cfg)
        return _json_response(Dict(
            "ok" => true, "state" => STATE[],
            "weight_sum" => current_weight_sum(STATE[]),
            "issues" => validate_state(STATE[]),
        ))
    catch err
        return _json_response(Dict("ok" => false, "error" => sprint(showerror, err)); status=500)
    end
end

function handle_discover_runs(req::HTTP.Request)
    try
        body = JSON.parse(String(req.body))
        parent = String(get(body, "parent", ""))
        isdir(parent) || return _json_response(Dict("runs" => String[],
            "error" => "parent folder not found"))
        runs = String[]
        for child in readdir(parent; join=true)
            isdir(child) && isdir(joinpath(child, "monte_carlo", "aggregated")) && push!(runs, child)
        end
        return _json_response(Dict("runs" => sort(runs)))
    catch err
        return _json_response(Dict("runs" => String[], "error" => sprint(showerror, err)))
    end
end

function handle_shared_endpoints(req::HTTP.Request)
    try
        body = JSON.parse(String(req.body))
        a = String(get(body, "a", ""))
        b = String(get(body, "b", ""))
        endpoints_for(d) = begin
            agg = joinpath(d, "monte_carlo", "aggregated")
            isdir(agg) || return Set{String}()
            Set(String[name for name in readdir(agg) if isdir(joinpath(agg, name))])
        end
        eps = collect(intersect(endpoints_for(a), endpoints_for(b)))
        return _json_response(Dict("endpoints" => sort(eps)))
    catch err
        return _json_response(Dict("endpoints" => String[], "error" => sprint(showerror, err)))
    end
end

function handle_synthesis_matrix(::HTTP.Request)
    try
        s = STATE[]; parent = s["synthesis_parent"]
        isdir(parent) || return _json_response(Dict("ok" => false,
            "error" => "Parent folder not found: $parent"); status=400)
        runs = String[]
        for child in readdir(parent; join=true)
            isdir(child) && isdir(joinpath(child, "monte_carlo", "aggregated")) && push!(runs, child)
        end
        isempty(runs) && return _json_response(Dict("ok" => false,
            "error" => "No completed runs found in $parent"); status=400)
        path = MovementFramework.build_pass_frequency_matrix(runs, s["synthesis_out_dir"])
        return _json_response(Dict("ok" => true,
            "message" => path === nothing ? "Matrix: no data" : "Matrix written: $path"))
    catch err
        return _json_response(Dict("ok" => false, "error" => sprint(showerror, err)); status=500)
    end
end

function handle_synthesis_diff(::HTTP.Request)
    try
        s = STATE[]
        a = s["scenario_a"]; b = s["scenario_b"]; ep = s["scenario_endpoint"]
        (isempty(a) || isempty(b) || isempty(ep)) &&
            return _json_response(Dict("ok" => false,
                "error" => "Need scenario A, scenario B, and endpoint"); status=400)
        res = MovementFramework.compare_scenarios(a, b, ep, s["synthesis_out_dir"])
        return _json_response(Dict("ok" => true,
            "message" => isempty(res) ? "No outputs (check inputs exist)" :
                          join(values(res), "; ")))
    catch err
        return _json_response(Dict("ok" => false, "error" => sprint(showerror, err)); status=500)
    end
end

# ============================================================================
# Router + entry
# ============================================================================

function make_router()
    router = HTTP.Router()
    HTTP.register!(router, "GET",  "/",                       handle_index)
    HTTP.register!(router, "GET",  "/api/state",              handle_get_state)
    HTTP.register!(router, "POST", "/api/state",              handle_set_state)
    HTTP.register!(router, "POST", "/api/list_vector_names",  handle_list_vector_names)
    HTTP.register!(router, "POST", "/api/inspect_layer_file", handle_inspect_layer_file)
    HTTP.register!(router, "GET",  "/api/preview_yaml",       handle_preview_yaml)
    HTTP.register!(router, "GET",  "/api/validate",           handle_validate)
    HTTP.register!(router, "POST", "/api/save_yaml",          handle_save_yaml)
    HTTP.register!(router, "POST", "/api/load_yaml",          handle_load_yaml)
    HTTP.register!(router, "POST", "/api/run_local",          handle_run_local)
    HTTP.register!(router, "POST", "/api/make_slurm",         handle_make_slurm)
    HTTP.register!(router, "POST", "/api/generate_package",   handle_generate_package)
    HTTP.register!(router, "POST", "/api/discover_runs",      handle_discover_runs)
    HTTP.register!(router, "POST", "/api/shared_endpoints",   handle_shared_endpoints)
    HTTP.register!(router, "POST", "/api/synthesis_matrix",   handle_synthesis_matrix)
    HTTP.register!(router, "POST", "/api/synthesis_diff",     handle_synthesis_diff)
    return router
end

function open_browser(url::AbstractString)
    try
        if Sys.iswindows()
            run(Cmd(["cmd", "/c", "start", "", url]); wait=false)
        elseif Sys.isapple()
            run(Cmd(["open", url]); wait=false)
        else
            run(Cmd(["xdg-open", url]); wait=false)
        end
    catch err
        @warn "Could not launch browser; open the URL manually" url exception=err
    end
end

function main()
    STATE[] = default_state()
    router = make_router()
    url = "http://localhost:$PORT"
    println("="^60)
    println(" Montescape GUI")
    println(" Open in your browser: $url")
    println(" (this terminal will keep showing logs - leave it open)")
    println(" Press Ctrl+C to quit.")
    println("="^60)
    @async (sleep(1.0); open_browser(url))
    HTTP.serve(router, Sockets.localhost, PORT)
end

end # module WebGUI
