# Spatially-correlated DEM error fields for Monte Carlo perturbation.
# Two methods:
#   - smoothed_gaussian : white noise * sigma, Gaussian-smoothed to the
#                         correlation length, variance-corrected to preserve sigma.
#   - spectral          : FFT-based unconditional Gaussian simulation with
#                         exponential covariance (closer to Heuvelink 1998).

using FFTW
using Random
using Statistics

"""
    perturb_dem(dem, res, cfg, rng) -> Matrix{Float32}

Return a copy of `dem` plus a spatially-correlated Gaussian error field.
Does nothing (returns the input) if cfg["enabled"] is false.
"""
function perturb_dem(dem::AbstractMatrix{<:Real}, res::Real,
                      cfg::AbstractDict, rng::AbstractRNG)
    get(cfg, "enabled", false) || return Float32.(dem)

    σ = Float64(get(cfg, "error_m", 5.0))
    L = Float64(get(cfg, "correlation_length_m", 150.0))
    method = get(cfg, "method", "smoothed_gaussian")

    err = method == "spectral" ?
        _spectral_field(size(dem), Float64(res), σ, L, rng) :
        _smoothed_gaussian_field(size(dem), Float64(res), σ, L, rng)

    return Float32.(dem) .+ Float32.(err)
end

function _smoothed_gaussian_field(shape, res, σ, L, rng)
    rows, cols = shape
    noise = randn(rng, rows, cols)
    # Gaussian kernel sized to L (in pixels)
    kernel_sigma_px = max(L / res, 1.0)
    smoothed = _gaussian_smooth(noise, kernel_sigma_px)
    # Smoothing reduces variance; rescale to recover requested σ
    s = std(smoothed)
    s == 0 && return zeros(rows, cols)
    return σ .* smoothed ./ s
end

function _gaussian_smooth(arr::AbstractMatrix{<:Real}, sigma_px::Real)
    rows, cols = size(arr)
    radius = max(Int(ceil(3 * sigma_px)), 1)
    kernel = [exp(-(x^2) / (2 * sigma_px^2)) for x in -radius:radius]
    kernel ./= sum(kernel)

    # 1D horizontal pass
    tmp = zeros(Float64, rows, cols)
    @inbounds for i in 1:rows, j in 1:cols
        s = 0.0
        for k in -radius:radius
            jj = clamp(j + k, 1, cols)
            s += arr[i, jj] * kernel[k + radius + 1]
        end
        tmp[i, j] = s
    end
    # 1D vertical pass
    out = zeros(Float64, rows, cols)
    @inbounds for i in 1:rows, j in 1:cols
        s = 0.0
        for k in -radius:radius
            ii = clamp(i + k, 1, rows)
            s += tmp[ii, j] * kernel[k + radius + 1]
        end
        out[i, j] = s
    end
    return out
end

function _spectral_field(shape, res, σ, L, rng)
    rows, cols = shape
    # Exponential covariance: C(h) = σ² * exp(-h/L)
    # PSD = FFT of covariance. We sample white noise in spectral domain
    # with amplitudes = sqrt(PSD), then IFFT.
    ky = fftfreq(rows, 1.0 / res)
    kx = fftfreq(cols, 1.0 / res)
    kk = [sqrt(kxi^2 + kyi^2) for kyi in ky, kxi in kx]
    # PSD of exponential covariance in 2D ~ 2π L² / (1 + (2π L k)²)^1.5
    psd = (2π * L^2) ./ (1.0 .+ (2π .* L .* kk) .^ 2) .^ 1.5
    spectrum = sqrt.(psd) .* (randn(rng, rows, cols) .+ im .* randn(rng, rows, cols))
    field = real(ifft(spectrum))
    s = std(field)
    s == 0 && return zeros(rows, cols)
    return σ .* field ./ s
end
