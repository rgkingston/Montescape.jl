# Custom Elephant terrain traversability layer.
#
# Mirrors movement_framework/src/layers/custom_elephant_terrain.py on main.
# Reclassifies slope into a 1-10 resistance scale following
# Wall et al. (2006) and Teixeira et al. (2025) Permissive Friction.
#
# Default thresholds (degrees) -> resistance (1-10):
#     <=15  -> 1   (flat to gentle)
#     15-25 -> 2   (moderate)
#     25-30 -> 4   (steep)
#     30-40 -> 6   (very steep)
#     40-50 -> 8   (near-limit)
#     50-60 -> 9   (extreme)
#     >60   -> 10  (impassable)
#
# Thresholds can be overridden via the `slope_thresholds` cfg entry.
# Output is normalised against `reference_value` (default 10) and clipped
# to [0, 2] so it composes correctly with the other friction layers.

using MovementFramework: register_layer!

function compute_elephant_terrain(dem::AbstractMatrix{<:Real}, res::Real, cfg::AbstractDict)
    raw = get(cfg, "slope_thresholds", Dict{Any,Any}(
        15 => 2, 25 => 4, 30 => 6, 40 => 8, 50 => 9, 60 => 10))
    thresholds = sort([(Float64(_to_num(k)), Float64(_to_num(v))) for (k, v) in raw];
                      by = first)
    ref = Float64(get(cfg, "reference_value", 10.0))

    rows, cols = size(dem)
    out = ones(Float32, rows, cols)
    r = Float64(res)
    @inbounds for i in 2:rows-1, j in 2:cols-1
        gy = (Float64(dem[i+1, j]) - Float64(dem[i-1, j])) / (2r)
        gx = (Float64(dem[i, j+1]) - Float64(dem[i, j-1])) / (2r)
        slope_deg = rad2deg(atan(sqrt(gx^2 + gy^2)))
        # apply thresholds in ascending order; last matching threshold wins
        v = 1.0
        for (thr, val) in thresholds
            if slope_deg > thr
                v = val
            end
        end
        out[i, j] = Float32(v)
    end
    return clamp.(out ./ Float32(ref), 0.0f0, 2.0f0)
end

_to_num(x::Number) = x
_to_num(x::AbstractString) = parse(Float64, x)

register_layer!("custom_elephant_terrain", compute_elephant_terrain)
