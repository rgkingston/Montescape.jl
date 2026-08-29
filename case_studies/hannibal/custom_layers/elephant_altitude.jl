# Custom Elephant altitude stress layer.
#
# Mirrors movement_framework/src/layers/custom_elephant_altitude.py on main.
# Independent of slope: based on elevation thresholds capturing hypoxia,
# cold exposure, reduced forage, and snow-line effects on war elephants in
# late-autumn alpine crossings (Polybius accounts, ~218 BCE).
#
# Default thresholds (metres) -> resistance (1-10):
#     <=1000     -> 1   (comfortable range)
#     1000-1500  -> 2   (mild stress)
#     1500-2000  -> 4   (moderate: cold, reduced forage)
#     2000-2500  -> 7   (high: snow line, hypoxia begins)
#     2500-3000  -> 9   (severe: near-lethal exposure)
#     >3000      -> 10  (extreme)
#
# Thresholds can be overridden via the `altitude_thresholds` cfg entry.

using MovementFramework: register_layer!

function compute_elephant_altitude(dem::AbstractMatrix{<:Real}, res::Real, cfg::AbstractDict)
    raw = get(cfg, "altitude_thresholds", Dict{Any,Any}(
        1000 => 2, 1500 => 4, 2000 => 7, 2500 => 9, 3000 => 10))
    thresholds = sort([(Float64(_to_num_alt(k)), Float64(_to_num_alt(v))) for (k, v) in raw];
                      by = first)
    ref = Float64(get(cfg, "reference_value", 10.0))

    rows, cols = size(dem)
    out = ones(Float32, rows, cols)
    @inbounds for i in 1:rows, j in 1:cols
        z = Float64(dem[i, j])
        v = 1.0
        for (thr, val) in thresholds
            if z > thr
                v = val
            end
        end
        out[i, j] = Float32(v)
    end
    return clamp.(out ./ Float32(ref), 0.0f0, 2.0f0)
end

_to_num_alt(x::Number) = x
_to_num_alt(x::AbstractString) = parse(Float64, x)

register_layer!("custom_elephant_altitude", compute_elephant_altitude)
