# Slope-based cost functions.
#
# All functions return Float32 array normalised against the cost at a
# reference slope of 30 deg, so values are comparable across study areas:
#   0.0 = flat / free
#   1.0 = reference 30 deg cost
#   up to 2.0 = beyond reference (clipped)

using Statistics

"""
    compute_slope_friction(dem, res, cfg) -> Matrix{Float32}

Dispatches on `cfg["method"]`. Supported methods:
  tobler, tobler_offpath, pandolf, minetti, herzog,
  llobera_sluckin, langmuir, irmischer_clarke, linear, exponential.
"""
function compute_slope_friction(dem::AbstractMatrix{<:Real}, res::Real, cfg::AbstractDict)
    method = get(cfg, "method", "tobler")
    slope_rad, slope_deg = _slope_arrays(dem, Float64(res))
    ref_rad = deg2rad(30.0)

    cost, ref_cost = _dispatch_cost(method, slope_rad, slope_deg, ref_rad, cfg)

    normalized = Float32.(cost ./ ref_cost)
    return clamp.(normalized, 0.0f0, 2.0f0)
end

function _slope_arrays(dem::AbstractMatrix{<:Real}, res::Float64)
    rows, cols = size(dem)
    gy = zeros(Float64, rows, cols)
    gx = zeros(Float64, rows, cols)
    @inbounds for i in 2:rows-1, j in 2:cols-1
        gy[i, j] = (Float64(dem[i+1, j]) - Float64(dem[i-1, j])) / (2 * res)
        gx[i, j] = (Float64(dem[i, j+1]) - Float64(dem[i, j-1])) / (2 * res)
    end
    slope_rad = atan.(sqrt.(gx .^ 2 .+ gy .^ 2))
    slope_deg = rad2deg.(slope_rad)
    return slope_rad, slope_deg
end

function _dispatch_cost(method, slope_rad, slope_deg, ref_rad, cfg)
    if method == "tobler"
        return _tobler(slope_rad, ref_rad, false)
    elseif method == "tobler_offpath"
        return _tobler(slope_rad, ref_rad, true)
    elseif method == "pandolf"
        return _pandolf(slope_rad, ref_rad, cfg)
    elseif method == "minetti"
        return _minetti(slope_rad, ref_rad)
    elseif method == "herzog"
        return _herzog(slope_rad, ref_rad)
    elseif method == "llobera_sluckin"
        return _llobera_sluckin(slope_rad, ref_rad)
    elseif method == "langmuir"
        return _langmuir(slope_rad, ref_rad)
    elseif method == "irmischer_clarke"
        return _irmischer_clarke(slope_rad, ref_rad, cfg)
    elseif method == "linear"
        return slope_deg ./ 90.0, 30.0 / 90.0
    elseif method == "exponential"
        k = Float64(get(cfg, "k", 0.1))
        return exp.(slope_deg .* k) .- 1.0, exp(30.0 * k) - 1.0
    else
        throw(ArgumentError("Unknown slope method: $method. " *
            "Valid: tobler, tobler_offpath, pandolf, minetti, herzog, " *
            "llobera_sluckin, langmuir, irmischer_clarke, linear, exponential."))
    end
end

# Tobler 1993. speed = 6 * exp(-3.5 * |tan(slope)+0.05|), off-path = *0.6.
function _tobler(slope_rad, ref_rad, offpath::Bool)
    f = offpath ? 0.6 : 1.0
    
    # 1. Core Tobler equation yielding speed in km/h
    speed_kmh = @. 6.0 * exp(-3.5 * abs(tan(slope_rad) + 0.05)) * f
    speed_kmh = max.(speed_kmh, 0.01)
    
    # 2. Convert km/h to m/s, then invert to get seconds per meter
    speed_ms = speed_kmh ./ 3.6
    cost = 1.0 ./ speed_ms

    # --- Reference Cost ---
    ref_speed_kmh = 6.0 * exp(-3.5 * abs(tan(ref_rad) + 0.05)) * f
    ref_speed_kmh = max(ref_speed_kmh, 0.01)
    
    ref_speed_ms = ref_speed_kmh / 3.6
    ref_cost = 1.0 / ref_speed_ms
    
    return cost, ref_cost
end

# Pandolf 1977 + Santee 2003 downhill correction.
# Uphill term uses +G; downhill subtracts Santee correction.
function _pandolf(slope_rad, ref_rad, cfg)
    W = Float64(get(cfg, "body_weight_kg", 80.0))
    L = Float64(get(cfg, "load_kg", 25.0))
    n = Float64(get(cfg, "terrain_factor", 1.0))
    V = Float64(get(cfg, "speed_ms", 1.3))

    # Gradients must remain SIGNED (% grade)
    slope_pct = tan.(slope_rad) .* 100.0   
    G_ref = tan(ref_rad) * 100.0

    # 1. Calculate Base Pandolf using signed slope
    # (For negative slopes, the 0.35 * V * G term correctly decreases cost)
    cost = @. 1.5 * W + 2.0 * (W + L) * (L / W)^2 +
              n * (W + L) * (1.5 * V^2 + 0.35 * V * slope_pct)

    # 2. Calculate the Santee 2003 Downhill Correction Factor (CF)
    cf = @. n * ((slope_pct * (W + L) * V) / 3.5 - (W + L) * (slope_pct + 6.0)^2 / W + 25.0 - V^2)
    
    # 3. Conditionally subtract CF ONLY for downhill steps (slope_pct < 0)
    cost = @. ifelse(slope_pct < 0.0, cost - cf, cost)
    cost = @. max(cost, 0.01)

    # --- Reference Cost Calculation ---
    ref_base = 1.5 * W + 2.0 * (W + L) * (L / W)^2 +
               n * (W + L) * (1.5 * V^2 + 0.35 * V * G_ref)
    
    ref_cf = n * ((G_ref * (W + L) * V) / 3.5 - (W + L) * (G_ref + 6.0)^2 / W + 25.0 - V^2)
    
    ref_cost = G_ref < 0.0 ? ref_base - ref_cf : ref_base
    ref_cost = max(ref_cost, 0.01)

    return cost, ref_cost
end

# Minetti 2002 walking metabolic polynomial (J kg-1 m-1).
function _minetti(slope_rad, ref_rad)
    g = tan.(slope_rad)
    cost = @. 280.5*g^5 - 58.7*g^4 - 76.8*g^3 + 51.9*g^2 + 19.6*g + 2.5
    cost = max.(cost, 0.01)
    gr = tan(ref_rad)
    ref_cost = 280.5*gr^5 - 58.7*gr^4 - 76.8*gr^3 + 51.9*gr^2 + 19.6*gr + 2.5
    return cost, ref_cost
end

# Herzog 2014 6th-degree polynomial fit to Minetti data (CORRECT version).
function _herzog(slope_rad, ref_rad)
    s = tan.(slope_rad)
    cost = @. 1337.8*s^6 + 278.19*s^5 - 517.39*s^4 - 78.199*s^3 +
              93.419*s^2 + 19.825*s + 1.64
    cost = max.(cost, 0.01)
    sr = tan(ref_rad)
    ref_cost = 1337.8*sr^6 + 278.19*sr^5 - 517.39*sr^4 - 78.199*sr^3 +
               93.419*sr^2 + 19.825*sr + 1.64
    return cost, ref_cost
end

# Llobera & Sluckin 2007 4th-degree.
function _llobera_sluckin(slope_rad, ref_rad)
    t = tan.(slope_rad)
    cost = @. 2.635 + 17.37*t + 42.37*t^2 - 21.43*t^3 + 14.93*t^4
    cost = max.(cost, 0.01)
    tr = tan(ref_rad)
    ref_cost = 2.635 + 17.37*tr + 42.37*tr^2 - 21.43*tr^3 + 14.93*tr^4
    return cost, ref_cost
end

# Langmuir 1984 piecewise time-in-seconds-per-metre.
# Naismith + Langmuir corrections, the r.walk formula but under our control.
function _langmuir(slope_rad, ref_rad)
    a = 0.72              # s / m horizontal base
    b_up = 6.0            # s per m of uphill gain
    c_gentle_dh = -2.0    # standard r.walk value (~ -1.9998)
    d_steep_dh  = 2.0     # standard r.walk value (~ 1.9998)

    s = tan.(slope_rad)
    sd = rad2deg.(slope_rad)
    cost = similar(s)

    @inbounds for i in eachindex(s)
        si = s[i]
        di = sd[i]
        
        if di >= 0.0
            # Uphill: Base + (uphill coefficient * fractional rise)
            cost[i] = a + b_up * si
        elseif di >= -5.0
            # Flat/negligible downhill (0 to -5 degrees): GRASS treats as flat base 'a'
            cost[i] = a
        elseif di >= -12.0
            # Gentle downhill (-5 to -12 degrees): Speed up.
            # We use abs(si) so the negative coefficient explicitly SUBTRACTS time.
            cost[i] = a + c_gentle_dh * abs(si)
        else
            # Steep downhill (< -12 degrees): Slow down.
            # Time decreases for the first 12 degrees, then increases for the remainder.
            # GRASS r.walk uses a continuous piecewise function:
            # For slopes steeper than -12 deg (tan(-12) approx -0.21256), it scales up
            s_12 = tan(deg2rad(-12.0))
            cost[i] = a + c_gentle_dh * abs(s_12) + d_steep_dh * (abs(si) - abs(s_12))
        end
    end
    cost = max.(cost, 0.01)

    # Fixed reference cost handling for the scalar input
    sr = tan(ref_rad)
    dr = rad2deg(ref_rad)
    if dr >= 0.0
        ref_cost = a + b_up * sr
    elseif dr >= -5.0
        ref_cost = a
    elseif dr >= -12.0
        ref_cost = a + c_gentle_dh * abs(sr)
    else
        s_12 = tan(deg2rad(-12.0))
        ref_cost = a + c_gentle_dh * abs(s_12) + d_steep_dh * (abs(sr) - abs(s_12))
    end
    ref_cost = max(ref_cost, 0.01)

    return cost, ref_cost
end
# Irmischer & Clarke 2018 on-path empirical fit.
function _irmischer_clarke(slope_rad, ref_rad, cfg)
    onpath = get(cfg, "on_path", true)
    female = get(cfg, "female", false)
    scale = female ? 0.95 : 1.0
    amp = onpath ? 1.0 : 0.67

    # Gradient in percent
    s_pct = tan.(slope_rad) .* 100.0
    
    # Speed calculated in km/h (0.95 scales the entire parenthesis block)
    speed_kmh = @. scale * (0.11 + amp * exp(-((s_pct + 5.0)^2) / (2.0 * 30.0^2)))
    speed_kmh = max.(speed_kmh, 0.001)
    
    # Convert km/h to m/s, then invert to get seconds per meter
    speed_ms = speed_kmh ./ 3.6
    cost = 1.0 ./ speed_ms

    # --- Reference Cost ---
    sr_pct = tan(ref_rad) * 100.0
    ref_speed_kmh = scale * (0.11 + amp * exp(-((sr_pct + 5.0)^2) / (2.0 * 30.0^2)))
    ref_speed_kmh = max(ref_speed_kmh, 0.001)
    
    ref_speed_ms = ref_speed_kmh / 3.6
    ref_cost = 1.0 / ref_speed_ms

    return cost, ref_cost
end
"""
    edge_slope_cost(slope::Float64, cfg::AbstractDict) -> Float64

Directional (signed) slope cost multiplier for a single edge.
slope = ?h / (dist * res), positive = uphill, negative = downhill.
Dispatches on cfg["method"], matching compute_slope_friction methods.
Returns cost relative to flat ground.
"""
function edge_slope_cost(slope::Float64, cfg::AbstractDict)
    method = get(cfg, "method", "tobler")

    if method == "tobler" || method == "tobler_offpath"
        f = method == "tobler_offpath" ? 0.6 : 1.0
        speed_kmh     = 6.0 * exp(-3.5 * abs(slope + 0.05)) * f
        ref_speed_kmh = 6.0 * exp(-3.5 * abs(0.0  + 0.05)) * f
        return max(ref_speed_kmh, 0.01) / max(speed_kmh, 0.01)

    elseif method == "pandolf"
        W = Float64(get(cfg, "body_weight_kg", 80.0))
        L = Float64(get(cfg, "load_kg", 25.0))
        n = Float64(get(cfg, "terrain_factor", 1.0))
        V = Float64(get(cfg, "speed_ms", 1.3))
        G = slope * 100.0
        cost = 1.5*W + 2.0*(W+L)*(L/W)^2 + n*(W+L)*(1.5*V^2 + 0.35*V*G)
        if G < 0.0
            cf = n*((G*(W+L)*V)/3.5 - (W+L)*(G+6.0)^2/W + 25.0 - V^2)
            cost -= cf
        end
        cost_ref = 1.5*W + 2.0*(W+L)*(L/W)^2 + n*(W+L)*(1.5*V^2)
        return max(cost, 0.01) / max(cost_ref, 0.01)

    elseif method == "minetti"
        g = slope
        cost     = 280.5*g^5 - 58.7*g^4 - 76.8*g^3 + 51.9*g^2 + 19.6*g + 2.5
        cost_ref = 2.5
        return max(cost, 0.01) / cost_ref

    elseif method == "herzog"
        s = slope
        cost     = 1337.8*s^6 + 278.19*s^5 - 517.39*s^4 - 78.199*s^3 +
                   93.419*s^2 + 19.825*s + 1.64
        cost_ref = 1.64
        return max(cost, 0.01) / cost_ref

    elseif method == "llobera_sluckin"
        t = slope
        cost     = 2.635 + 17.37*t + 42.37*t^2 - 21.43*t^3 + 14.93*t^4
        cost_ref = 2.635
        return max(cost, 0.01) / cost_ref

    elseif method == "langmuir"
        a           = 0.72
        b_up        = 6.0
        c_gentle_dh = -2.0
        d_steep_dh  =  2.0
        sd = rad2deg(atan(abs(slope)))
        if slope >= 0.0
            cost = a + b_up * slope
        elseif sd <= 5.0
            cost = a
        elseif sd <= 12.0
            cost = a + c_gentle_dh * abs(slope)
        else
            s_12 = tan(deg2rad(12.0))
            cost = a + c_gentle_dh * abs(s_12) + d_steep_dh * (abs(slope) - s_12)
        end
        cost_ref = a
        return max(cost, 0.01) / cost_ref

    elseif method == "irmischer_clarke"
        onpath = get(cfg, "on_path", true)
        female = get(cfg, "female", false)
        scale  = female ? 0.95 : 1.0
        amp    = onpath ? 1.0 : 0.67
        s_pct  = slope * 100.0
        speed_kmh     = scale * (0.11 + amp * exp(-((s_pct + 5.0)^2) / (2.0 * 30.0^2)))
        ref_speed_kmh = scale * (0.11 + amp * exp(-(( 0.0  + 5.0)^2) / (2.0 * 30.0^2)))
        return max(ref_speed_kmh, 0.001) / max(speed_kmh, 0.001)

    elseif method == "linear"
        cost_ref = 30.0 / 90.0
        return max(abs(slope) * 100.0 / 90.0, 0.01) / cost_ref

    elseif method == "exponential"
        k = Float64(get(cfg, "k", 0.1))
        cost_ref = exp(30.0 * k) - 1.0
        return max(exp(abs(slope) * 100.0 * k) - 1.0, 0.01) / max(cost_ref, 0.01)

    else
        throw(ArgumentError("Unknown slope method: $method"))
    end
end