########################################################################
#                                                                      #
# DOT_RydSim/src/Schrödinger.mod.jl                                    #
#                                                                      #
# (c) Dirk Oliver Theis 2023                                           #
#                                                                      #
# License:                                                             #
#                                                                      #
#             Apache 2.0                                               #
#                                                                      #
########################################################################


"""
Module `Schrödinger`

Simulation of multi-qubit quantum evolution under time-dependent Hamiltonian.

# Exports

## General
* Function [`schröd!()`](@ref)

"""
module Schrödinger
export schröd!

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 0. ToC  +
#                                                                                                             |
#  Table of Contents                                                                                          |
#  —————————————————                                                                                          |
#                                                                                                             |
#                                                                                                             |
#  1. Imports & Helpers                                                                                       |
#                                                                                                             |
#     1.1. X,N matrices                                                                                       |
#     1.2. Shit                                                                                               |
#                                                                                                             |
#                                                                                                             |
#  2. Types                                                                                                   |
#                                                                                                             |
#     2.1. Exceptions                                                                                         |
#                                                                                                             |
#                                                                                                             |
#  3. Work horses                                                                                             |
#                                                                                                             |
#     3.1. Function `timestep!()`                                                                             |
#     3.2. Function `schröd!()`                                                                               |
#                                                                                                             |
#—————————————————————————————————————————————————————————————————————————————————————————————————————————————+


# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1. Imports & Helpers

import ..μs_t, ..Rad_per_μs_t, ..Radperμs_per_μs_t
import ..Pulse, ..phase, ..𝑎𝑣𝑔, ..𝑠𝑡𝑒𝑝
using  ..DOT_NiceMath

using LinearAlgebra: Hermitian, I as Id


using Unitful
using Unitful: μs

using Logging

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.1: X,Z matrices
#

@doc raw"""
Function `N₁(A ::Integer, ℂ ::Type{<Complex})` — Rydberg |r⟩⟨r| operator for the 1st atom

Letter is customarily `N` as it's kinda a number operator, except of course, here, |g⟩=|1⟩
and |r⟩=|2⟩ — so it ain't exactly no *number* operator.

## Arguments:
1. `A` — number of atoms
2. `ℂ` — complex number type to use as `eltype`.

### Returns:
* 2ᴬ×2ᴬ Hermitian matrix; type `Hermitian{ℂ,Matrix{ℂ}}`
"""
N₁(A ::Integer, ℂ::Type{<:Complex}) ::Hermitian{ℂ,Matrix{ℂ}} = begin
    @assert A ≥ 1
    let 𝙽 = [  ℂ(0)      0
                0       +1  ]

        𝙽 ⊗ Id(2^(A-1)) |> Hermitian
    end
end

(
    X₁(A ::Integer ; γ::ℂ)   ::Hermitian{ℂ,Matrix{ℂ}}
) where{ℂ} =
    begin
        @assert A ≥ 1
        let 𝚇 = [   0     γ
                    γ'    0   ]

            𝚇 ⊗ Id(2^(A-1)) |> Hermitian
        end
    end

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.2: Shit
#

log_of_pow2(x ::Integer) = begin @assert ispow2(x) ; trailing_zeros(x) end


# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2. Types

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.1: Ctrl_Exception

struct Ctrl_Exception <: Exception
    msg ::String
end

import Base: showerror
showerror(io::IO, e::Ctrl_Exception) = print(io, "schröd!(): Bad quantum ctrl data: ",e.msg)


# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 3. Work horses

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 3.1: timestep!()
function  timestep!(ψ  ::Vector{ℂ},
                    𝛥𝑡 ::μs_t{ℝ}
                    ;
                    𝜔  ::Rad_per_μs_t{ℝ},
                    X  ::Hermitian{ℂ,𝕄_t},
                    𝛿  ::Rad_per_μs_t{ℝ},
                    N  ::Hermitian{ℂ,𝕄_t},
                    R  ::Hermitian{ℂ,𝕄_t}            ) ::Nothing      where{ℝ,ℂ,𝕄_t}

    let ωΔt ::ℝ              = ustrip(NoUnits, 𝜔⋅𝛥𝑡),
        δΔt ::ℝ              = ustrip(NoUnits, 𝛿⋅𝛥𝑡),
        Δt  ::ℝ              = ustrip(μs, 𝛥𝑡),
        A   ::Hermitian{ℂ,𝕄_t} = ωΔt⋅X - δΔt⋅N + Δt⋅R

        ψ .= cis(A)'ψ
    end
    nothing
end #^ timestep!()

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 3.2: schröd!()
@doc raw"""
Function `schröd!(ψ  ::Vector{ℂ},  𝑇 ::μs_t{ℝ} ; ...)   where{ℝ,ℂ,𝕄_t}`

Simulates time evolution under time-dependent Hamiltonian
```math
    H(t)/\hbar = \omega(t) Xᵧ - \delta(t) N + R
```

where `Xᵧ` is defined as |1⟩⟨2|⋅γ + |2⟩⟨1|⋅γ̄.  The phase, γ, is not time-dependent.

The number of atoms, `A`, is derived from the dimensions of ``R``, which is assumed to be
a 2ᴬ×2ᴬ matrix.  The matrices ``X,N`` are constructed to act on the "first" atom, i.e.,
X⊗1⊗...⊗1 and N⊗1⊗...⊗1; see the function `⊗` above.

## Arguments

Note the convention: *Italic variables are unitful*. That also applies to functions: The
returned value of an italic-type function must be unitful.

### Position variables
1. `ψ` — state vector at time 𝑡=0μs; will be modified in place.
2. `𝑇` — end-time of evolution, in μs.

### Returns
... nothing.

### Keyword arguments

#### Mandatory keyword arguments: Defining the Hamiltonian

* `Ω ::Pulse` — Defines the Rabi-frequency pulse shape, including the phase, γ
* `Δ ::Pulse` — Defines the detuning pulse shape.
* `R ::Hermitian{ℂ,𝕄_t}` — Rydberg interaction term for all atoms.

#### Optional keyword arguments
* `ε ::ℝ` — (`\varepsilon`) Simulation accuracy; determines size of time steps.
"""
function schröd!(ψ  ::Vector{ℂ},
                 𝑇  ::μs_t{ℝ}
                 ;
                 Ω  ::P₁,
                 Δ  ::P₂,
                 R  ::Hermitian{ℂ,𝕄_t},
                 ε  ::ℝ                   = ℝ(1e-3) ) ::Nothing   where{ℝ,ℂ,𝕄_t, P₁<:Pulse, P₂<:Pulse}

    @assert ε > 0   "ε ≤ 0"

    A    = log_of_pow2( length(ψ) )       ; @assert A ≥ 1               "Need at least one atom, i.e., length ψ ≥ 2."
    𝟐ᴬ   = length(ψ)                      ; @assert 2^A == 𝟐ᴬ           "Crazy bug #1"
    N    = N₁(A,ℂ)                        ; @assert size(N) == size(R)  "Sizes of `ψ` and `R` don't match."
    X    = X₁(A;γ=phase(Ω))               ; @assert size(X) == size(N)  "Crazy bug #2"


    warn_RWA_count = 0
    warn_RWA_out   = (1:3)∪(10:10:30)∪(100:100:300)∪(1000:1000:3000)∪(10000:10000:30000)

    𝑡 ::μs_t{ℝ} = 0μs

    while 𝑡  <  𝑇 - 1e-50μs

        𝑠Ω = 𝑠𝑡𝑒𝑝(Ω, 𝑡 ; ε )
        𝑠Δ = 𝑠𝑡𝑒𝑝(Δ, 𝑡 ; ε )

        Ω_𝛥𝑡 ::μs_t{ℝ} = min(𝑇-𝑡, 𝑠Ω)
        Δ_𝛥𝑡 ::μs_t{ℝ} = min(𝑇-𝑡, 𝑠Δ)

        𝛥𝑡 = min( Ω_𝛥𝑡, Δ_𝛥𝑡 )

        Ω_𝜇 ::Rad_per_μs_t{ℝ} = 𝑎𝑣𝑔(Ω, 𝑡; 𝛥𝑡)
        Δ_𝜇 ::Rad_per_μs_t{ℝ} = 𝑎𝑣𝑔(Δ, 𝑡; 𝛥𝑡)



        𝑠Ω > 0μs ||
            throw(Ctrl_Exception("Time-step at 𝑡=$(BigFloat(𝑡)) for Ω is non-positive: \
                                          $(BigFloat(𝑠Ω)) ≤ 0μs"))
        𝑠Δ > 0μs ||
            throw(Ctrl_Exception("Time-step at 𝑡=$(BigFloat(𝑡)) for Δ is non-positive: \
                                          $(BigFloat(𝑠Δ)) ≤ 0μs"))
        Ω_𝜇 ≥ 0/μs ||
            throw(Ctrl_Exception("At time 𝑡=$(BigFloat(𝑡)) Ω is negative: \
                                          $(BigFloat(Ω_𝜇s)) < 0/μs"))
        if abs(Δ_𝜇) > 1e-6/μs && Ω_𝜇 < 1e-10/μs
            warn_RWA_count += 1
            if n_warn_RWA ∈ warn_RWA_out
                @warn   "\
                        schröd!(): RWA break-down[$(n_warn_RWA)]: \
                        Δ is non-zero, but Ω is very small: \
                        " 𝑡=BigFloat(𝑡) Δ=BigFloat(Δ_𝜇) Ω=BigFloat(Ω_𝜇)
            end
        end

        timestep!(ψ, 𝛥𝑡 ; 𝜔=Ω_𝜇, 𝛿=Δ_𝜇,
                          X, N, R)

        𝑡 += 𝛥𝑡

    end #^ while 𝑡
    nothing
end #^ schröd!()


end # module Schrödinger

#
# Thoughts about the algorithm
#
# Matrix exponential
#
# For matrices up to 4×4:       cis( StaticMatrix )           ~1 μs  (vs ~2μs)
#
# For matrices from 16×16       cis( Hermitian(Matrix) )     ~26 μs (vs ∞)
#
#
#EOF
