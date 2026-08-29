# Aspect-based directional preference penalty.
# Returns 0 where the slope faces the preferred direction, 1 where it
# faces the opposite direction.

function compute_aspect_penalty(dem::AbstractMatrix{<:Real}, res::Real, cfg::AbstractDict)
    rows, cols = size(dem)
    gy = zeros(Float64, rows, cols)
    gx = zeros(Float64, rows, cols)
    r = Float64(res)
    @inbounds for i in 2:rows-1, j in 2:cols-1
        gy[i, j] = (Float64(dem[i+1, j]) - Float64(dem[i-1, j])) / (2 * r)
        gx[i, j] = (Float64(dem[i, j+1]) - Float64(dem[i, j-1])) / (2 * r)
    end
    # Aspect: 0 = north, increasing clockwise
    aspect = mod.(rad2deg.(atan.(-gx, gy)), 360.0)
    preferred = Float64(get(cfg, "preferred_direction", 180.0))
    diff = abs.(aspect .- preferred)
    diff = min.(diff, 360.0 .- diff)        # shortest angular distance
    return Float32.(diff ./ 180.0)
end
