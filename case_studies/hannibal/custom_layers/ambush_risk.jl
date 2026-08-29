# Custom Ambush risk layer for Hannibal's crossing.
#
# West-to-east risk gradient: peak risk on the western edge
#  low risk in the east  allied
#  Extra penalty in low-elevation western valleys where
# ambushes were physically easiest to stage.

using MovementFramework: register_layer!

function compute_ambush_risk(dem::AbstractMatrix{<:Real}, res::Real, cfg::AbstractDict)
    west_frac    = Float64(get(cfg, "west_fraction", 0.0))
    east_frac    = Float64(get(cfg, "east_fraction", 0.65))
    peak         = Float64(get(cfg, "peak_penalty", 8.0))
    safe         = Float64(get(cfg, "safe_penalty", 1.0))
    valley_thr   = Float64(get(cfg, "valley_elevation_threshold", 300.0))
    valley_bonus = Float64(get(cfg, "valley_bonus", 2.0))
    sigma_px     = Float64(get(cfg, "smoothing_sigma", 5.0))
    ref          = Float64(get(cfg, "reference_value", 10.0))

    rows, cols = size(dem)
    out = zeros(Float32, rows, cols)
    midpoint = cld(cols, 2)
    @inbounds for i in 1:rows, j in 1:cols
        f = (j - 1) / max(cols - 1, 1)
        t = if f <= west_frac
            0.0
        elseif f >= east_frac
            1.0
        else
            (f - west_frac) / (east_frac - west_frac)
        end
        cost = peak + (safe - peak) * t
        z = Float64(dem[i, j])
        if 0 < z < valley_thr && j < midpoint
            cost += valley_bonus
        end
        out[i, j] = Float32(clamp(cost, safe, peak + valley_bonus))
    end
    return clamp.(_gauss_smooth_2d(out, sigma_px) ./ Float32(ref), 0.0f0, 2.0f0)
end

function _gauss_smooth_2d(arr::AbstractMatrix{<:Real}, sigma_px::Real)
    radius = max(Int(ceil(3 * sigma_px)), 1)
    kernel = [exp(-(x^2) / (2 * sigma_px^2)) for x in -radius:radius]
    kernel ./= sum(kernel)
    rows, cols = size(arr)
    tmp = zeros(Float64, rows, cols)
    @inbounds for i in 1:rows, j in 1:cols
        s = 0.0
        for k in -radius:radius
            jj = clamp(j + k, 1, cols)
            s += arr[i, jj] * kernel[k + radius + 1]
        end
        tmp[i, j] = s
    end
    out = zeros(Float32, rows, cols)
    @inbounds for i in 1:rows, j in 1:cols
        s = 0.0
        for k in -radius:radius
            ii = clamp(i + k, 1, rows)
            s += tmp[ii, j] * kernel[k + radius + 1]
        end
        out[i, j] = Float32(s)
    end
    return out
end

register_layer!("custom_ambush_risk", compute_ambush_risk)
