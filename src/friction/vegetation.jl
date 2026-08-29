# WARNING
# has not been fully tested 
# Vegetation
# / land cover layer. Reclassifies integer codes to friction values.

function compute_vegetation(arr::AbstractMatrix, cfg::AbstractDict)
    mapping = get(cfg, "reclass", Dict{Any,Any}())
    out = zeros(Float32, size(arr))
    if isempty(mapping)
        return out
    end
    intmap = Dict{Int,Float32}()
    for (k, v) in mapping
        intmap[parse_or_int(k)] = Float32(v)
    end
    @inbounds for i in eachindex(arr)
        out[i] = get(intmap, Int(round(arr[i])), 0.0f0)
    end
    return out
end

parse_or_int(k::Integer) = Int(k)
parse_or_int(k::AbstractString) = parse(Int, k)
