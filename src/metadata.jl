# Per-run reproducibility metadata (port of framework/metadata.py).
#
# Writes <run_dir>/run_metadata.json capturing the config used, the
# results summary, host/platform information, the Julia version, and

using JSON
using Dates

"""
    save_run_metadata(cfg, run_dir, results) -> String

Write run_metadata.json into `run_dir` and return its path.
"""
function save_run_metadata(cfg::AbstractDict, run_dir::AbstractString, results::AbstractDict)
    meta = Dict{String,Any}(
        "timestamp"   => string(now()),
        "platform"    => string(Sys.MACHINE, " / ", Sys.KERNEL),
        "julia"       => string(VERSION),
        "config"      => _sanitize(cfg),
        "results"     => _sanitize(results),
        "output_dir"  => run_dir,
    )
    out = joinpath(run_dir, "run_metadata.json")
    open(out, "w") do io
        JSON.print(io, meta, 2)
    end
    @info "  [meta] Saved: $out"
    return out
end

_sanitize(obj) = obj
_sanitize(::Nothing) = nothing
_sanitize(t::Tuple) = collect(_sanitize.(t))
_sanitize(v::AbstractVector) = [_sanitize(x) for x in v]
function _sanitize(d::AbstractDict)
    out = Dict{String,Any}()
    for (k, v) in d
        ks = String(string(k))
        startswith(ks, "_") && continue   # skip internal keys like _project_root
        out[ks] = _sanitize(v)
    end
    return out
end
