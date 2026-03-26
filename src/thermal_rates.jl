# Copyright (c) 2026 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

# Thermal rates and effective potential for HNL production
# PDF §2.1, Eqs. (2.3)–(2.6)
#
# c_α(T) from numerical data in hep-ph/0612182 (Ref. [13], data from Ref. [18])
# V_α(T,p) from Eq. (2.4)
# d_α(T) ≡ -V_α / Γ_α from Eq. (2.5) — momentum-independent

# ─── Load numerical c_α(T) data from hep-ph/0612182 ───

data_from_hep_ph_0612182_directory = joinpath(external_data_directory, "from_hep-ph.0612182") |> ensure_directory

const _hatIQ_files = [
    joinpath(data_from_hep_ph_0612182_directory, "hatIQ_M001_alpha1_Fermi.dat"),
    joinpath(data_from_hep_ph_0612182_directory, "hatIQ_M001_alpha2_Fermi.dat"),
    joinpath(data_from_hep_ph_0612182_directory, "hatIQ_M001_alpha3_Fermi.dat"),
]

"""
    _load_hatIQ(filepath) -> (T_MeV, qoverT, hatIQ) matrices

Load hat{I}_Q data from hep-ph/0612182 data files. Returns raw arrays.
"""
function _load_hatIQ(filepath::AbstractString)
    lines = readlines(filepath)
    # Skip comment lines starting with #
    data_lines = filter(l -> !startswith(strip(l), "#") && !isempty(strip(l)), lines)

    T_vals = Float64[]
    q_vals = Float64[]
    I_vals = Float64[]

    for line ∈ data_lines
        parts = split(strip(line))
        length(parts) ≥ 3 || continue
        push!(T_vals, parse(Float64, parts[1]))  # T in MeV
        push!(q_vals, parse(Float64, parts[2]))  # q/T
        push!(I_vals, parse(Float64, parts[3]))  # hat{I}_Q = c_α
    end

    return T_vals, q_vals, I_vals
end

"""
    _build_cα_interpolation(filepath; xp=3.0)

Build a c_α(T) interpolation from hat{I}_Q data at fixed q/T = xp.
Returns a function T::EnergyUnit → c_α (dimensionless Float64).
"""
function _build_cα_interpolation(filepath::AbstractString; xp::Float64=3.0)
    T_vals, q_vals, I_vals = _load_hatIQ(filepath)

    # Find the q/T value closest to xp
    # Maybe interpolate in q/T if needed, but for now just pick the closest slice
    unique_q = sort(unique(q_vals))
    _, idx = findmin(abs.(unique_q .- xp))
    qT_target = unique_q[idx]

    # Filter to the selected q/T slice
    mask = q_vals .== qT_target
    T_slice = T_vals[mask]   # in MeV
    c_slice = I_vals[mask]

    # Sort by T (ascending for interpolation)
    perm = sortperm(T_slice)
    T_sorted = T_slice[perm]
    c_sorted = c_slice[perm]

    # Build log-T interpolation
    interp = LinearInterpolation(
        c_sorted, log.(T_sorted);
        extrapolation=ExtrapolationType.Constant, cache_parameters=true
    )

    return interp
end

# Build interpolations for each flavor at q/T = 3
const _cα_interps = [_build_cα_interpolation(f; xp=3.0) for f ∈ _hatIQ_files]

export cα
"""
    cα(α::Int, T::EnergyUnit) -> Float64

Neutrino production rate coefficient c_α for flavor α ∈ {1,2,3} = {e,μ,τ}
at temperature T, evaluated at q/T = 3.

Γ_α = c_α * G_F² * p * T⁴   [Eq. (2.3)]
"""
function cα(α::Int, T::EnergyUnit)
    @assert 1 ≤ α ≤ 3 "Flavor index must be 1 (e), 2 (μ), or 3 (τ)"
    log_T_in_MeV = T |> EUval(MeV) |> log
    return _cα_interps[α](log_T_in_MeV)
end

# To AI agent:
# We can apply the implementation in `refs/code-2507.12270` for computing ρ_ℓ and ρ_ν.
# By Quan-feng Wu, 2026-03-26

# ─── Energy density of charged lepton ℓ_α ───

"""
    ρ_ℓ(α::Int, T::EnergyUnit) -> EnergyUnit (dim 4)

Energy density of charged lepton ℓ_α (particle + antiparticle) at temperature T.
Uses numerical integration for arbitrary m_ℓ/T ratio.

ρ_ℓ = (4 / 2π²) ∫₀^∞ u² √(u² + x²) / (e^√(u²+x²) + 1) du × T⁴

where x = m_ℓ/T, and g_ℓ = 4 (2 spin × particle/antiparticle).
"""
function ρ_ℓ(α::Int, T::EnergyUnit)
    x = EUval(m_ℓ[α] / T)  # m_ℓ / T, dimensionless

    # For very heavy leptons, Boltzmann-suppressed
    x > 50 && return zero(EU) * EU(1, 3)  # return 0 with dim 4

    # Numerical integration via simple Gauss-Laguerre-like quadrature
    # Transform: u from 0 to ∞
    # Use adaptive Simpson or simple trapezoidal on [0, u_max]
    u_max = max(30.0, 3 * x + 30.0)
    N_pts = 500
    du = u_max / N_pts

    integral = 0.0
    for i in 1:N_pts
        u = (i - 0.5) * du
        E_u = sqrt(u^2 + x^2)
        integral += u^2 * E_u / (exp(E_u) + 1) * du
    end

    g_ℓ = 4.0  # 2 spin × 2 (particle + antiparticle)
    return EU(g_ℓ / (2π^2) * integral, 0) * T^4
end

# ─── Energy density of neutrino ν_α ───

"""
    ρ_ν(T::EnergyUnit) -> EnergyUnit (dim 4)

Energy density of a single neutrino species (ν_α + ν̄_α) at temperature T.
Massless approximation: ρ_ν = (7/8) × (2/g_total) × (π²/30) × T⁴

g_ν = 2 (left-handed ν + right-handed ν̄), so ρ_ν = 2 × (7/8) × (π²/30) × T⁴
"""
function ρ_ν(T::EnergyUnit)
    return EU(2 * (7 / 8) * π^2 / 30, 0) * T^4
end

# ─── Effective potential V_α ───

export Vα
"""
    Vα(α::Int, T::EnergyUnit, p::EnergyUnit) -> EnergyUnit (dim 1)

Effective matter potential for ν_α-N mixing at temperature T and momentum p.
Eq. (2.4): V_α = -(8√2/3) G_F (ρ_{ν_α}/m_Z² + ρ_{ℓ_α}/m_W²) p
"""
function Vα(α::Int, T::EnergyUnit, p::EnergyUnit)
    return -(8sqrt(2) / 3) * G_F * (ρ_ν(T) / m_Z^2 + ρ_ℓ(α, T) / m_W^2) * p
end

# ─── Ratio d_α ≡ -V_α / Γ_α ───

export dα
"""
    dα(α::Int, T::EnergyUnit) -> Float64

Ratio d_α ≡ -V_α / Γ_α.  Eq. (2.5).
Independent of momentum p (cancels in ratio).

d_α = (8√2 / 3) × (ρ_{ν_α}/m_Z² + ρ_{ℓ_α}/m_W²) / (c_α G_F T⁴)
"""
function dα(α::Int, T::EnergyUnit)
    c = cα(α, T)
    rho_term = ρ_ν(T) / m_Z^2 + ρ_ℓ(α, T) / m_W^2
    return (8sqrt(2) / 3) * EUval(rho_term / (G_F * T^4)) / c
end

# ─── Production rate Γ_α ───

export Γα
"""
    Γα(α::Int, T::EnergyUnit, p::EnergyUnit) -> EnergyUnit (dim 1)

Active neutrino production rate Γ_α = c_α G_F² p T⁴.  Eq. (2.3).
"""
function Γα(α::Int, T::EnergyUnit, p::EnergyUnit)
    return EU(cα(α, T), 0) * G_F^2 * p * T^4
end

export ρ_ℓ, ρ_ν
