# Terrain Ruggedness Index (Riley et al. 1999).
# Mean absolute difference between centre cell and its 8 neighbours.
# Normalised against Riley's "extremely rugged" reference of 80 m.

function compute_tri(dem::AbstractMatrix{<:Real})
    rows, cols = size(dem)
    out = zeros(Float32, rows, cols)
    @inbounds for i in 2:rows-1, j in 2:cols-1
        c = Float64(dem[i, j])
        s = 0.0
        s += abs(Float64(dem[i-1, j-1]) - c)
        s += abs(Float64(dem[i-1, j  ]) - c)
        s += abs(Float64(dem[i-1, j+1]) - c)
        s += abs(Float64(dem[i  , j-1]) - c)
        s += abs(Float64(dem[i  , j+1]) - c)
        s += abs(Float64(dem[i+1, j-1]) - c)
        s += abs(Float64(dem[i+1, j  ]) - c)
        s += abs(Float64(dem[i+1, j+1]) - c)
        out[i, j] = Float32(s / 8.0)
    end
    # Riley 1999: TRI=80m is "extremely rugged"
    return clamp.(out ./ 80.0f0, 0.0f0, 2.0f0)
end
