# Copyright (c) 2026 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

# Physical constants used across all notebooks
# Reference: PDG 2024, Ref. [18] (hep-ph/0612182)

export EU, NU, M_Pl
const EU = GeV
const NU = NaturalUnit(EU)
const M_Pl = NU.M_Pl

# Fermi constant
export G_F
const G_F = GeV(1.1663785e-5, -2)

# Gauge boson masses
export m_W, m_Z
const m_W = GeV(80.3692)
const m_Z = GeV(91.1880)

# Lepton masses
export m_e, m_μ, m_τ, m_ℓ
const m_e = MeV(0.51099895)
const m_μ = MeV(105.6583755)
const m_τ = MeV(1776.93)
# Lepton masses array indexed by flavor: 1→e, 2→μ, 3→τ
const m_ℓ = [m_e, m_μ, m_τ]

# QCD scale
export Λ_QCD
const Λ_QCD = MeV(200)

# Nucleon masses and parameters
export m_p, m_n, Q_np, τ_n
const m_p = MeV(938.27208816)
const m_n = MeV(939.5654205)
const Q_np = m_n - m_p
const τ_n = 878.4 * NU.s

# Fine-structure constant
export α_EM
const α_EM = inv(137.035999177)

# Baryon-to-photon ratio
export η_B
const η_B = 6.1e-10

