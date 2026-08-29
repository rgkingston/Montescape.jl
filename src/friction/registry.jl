# Plugin-style registry for custom DEM-derived layers.
#
# Users register a function (dem, res, cfg) -> Matrix{Float32} under a name,
# and reference that name in their YAML config under friction.custom_layers.
#
# No case-study code lives in the core package.

const _LAYER_REGISTRY = Dict{String, Function}()

"""
    register_layer!(name::AbstractString, fn::Function)

Register a custom friction layer. The function must accept
`(dem::Matrix{Float32}, res::Float64, cfg::AbstractDict)` and return a
`Matrix{Float32}` normalised so 1.0 represents the reference difficulty
(typically clipped to [0, 2]).

Example:
    register_layer!("my_layer", (dem, res, cfg) -> Float32.(...))
"""
function register_layer!(name::AbstractString, fn::Function)
    _LAYER_REGISTRY[String(name)] = fn
    return name
end

"""
    lookup_layer(name) -> Function

Throws ArgumentError listing available layers if not found.
"""
function lookup_layer(name::AbstractString)
    haskey(_LAYER_REGISTRY, name) && return _LAYER_REGISTRY[name]
    throw(ArgumentError(
        "Unknown custom layer '$name'. Available: $(sort(collect(keys(_LAYER_REGISTRY)))). " *
        "Make sure the .jl file that calls register_layer! is listed in " *
        "friction.custom_layer_files."))
end

"""
    list_registered_layers() -> Vector{String}
"""
list_registered_layers() = sort(collect(keys(_LAYER_REGISTRY)))

"""
    load_custom_layer_files!(paths)

Include each .jl file; each is expected to call register_layer!.
Paths can be absolute or relative to the project root.
"""
function load_custom_layer_files!(paths, project_root::AbstractString=pwd())
    for p in paths
        p === nothing && continue
        full = isabspath(String(p)) ? String(p) : abspath(joinpath(project_root, String(p)))
        isfile(full) || throw(ArgumentError("custom_layer_files entry not found: $full"))
        Base.include(Main, full)
    end
end
