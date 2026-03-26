# Copyright (c) 2026 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using DelimitedFiles
using Integrals

include(joinpath(dirname(@__DIR__), "src", "utilities.jl"))
include(joinpath(dirname(@__DIR__), "src", "directories.jl"))

# ─── Read parameters with defaults ───

function _read_param(prompt::String, default)
    print(prompt, " [", default, "]: ")
    input = strip(readline(stdin))
    return isempty(input) ? default : parse(typeof(default), input)
end

# ─── Precompute massive Fermi-Dirac energy density integral ───

const _FD_integral_file = joinpath(integrals_directory, "integrals_massive_FD.csv")

function _generate_FD_integrals(; Tₘ_start=1.5e-3, Tₘ_end=1e15, n_pts=10000)
    @info "Precomputing massive Fermi-Dirac integrals ($(n_pts) points, T/m ∈ [$(Tₘ_start), $(Tₘ_end)])..."

    f_FD(pₘ, Tₘ) = inv(exp(sqrt(pₘ^2 + 1) / Tₘ) + 1)
    integrand_ρ_FD(pₘ, Tₘ) = sqrt(pₘ^2 + 1) * pₘ^2 * f_FD(pₘ, Tₘ) / (2π^2)

    Tₘ_list = geomspace(Tₘ_start, Tₘ_end, n_pts)
    ρ_FD_list = zeros(length(Tₘ_list))

    for (i, Tₘ) ∈ enumerate(Tₘ_list)
        prob = IntegralProblem(integrand_ρ_FD, (0.0, Inf), Tₘ)
        sol = solve(prob, QuadGKJL())
        ρ_FD_list[i] = sol.u
    end

    open(_FD_integral_file, "w") do io
        println(io, "# T_over_m,energy_density_FD")
        writedlm(io, hcat(Tₘ_list, ρ_FD_list), ',')
    end

    @info "Saved to $(_FD_integral_file)"

    return nothing
end

isfile(_FD_integral_file) || _generate_FD_integrals()
