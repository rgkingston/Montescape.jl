# # WARNING
# has not been fully tested and completed
# Hydrology layer: a fixed crossing penalty wherever the input raster > 0.

function compute_hydrology(arr::AbstractMatrix, cfg::AbstractDict)
    penalty = Float32(get(cfg, "crossing_penalty", 8.0))
    return Float32.(ifelse.(arr .> 0, penalty, 0.0f0))
end
