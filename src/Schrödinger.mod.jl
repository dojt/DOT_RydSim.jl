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
* (nothing yet)

## Number definitions (in sub-modules `Numbers`𝑥𝑦𝑧)
* (also nothing)

# Sub-modules

Sub-module names are not exported by `DOT_RydSim`.

* `DOT_NiceMath.` — nope...
"""
module Schrödinger
export no_thing

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 0. TOC

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1. Imports & Helpers

import ..μs_t
import ..Rad_per_μs_t
using  ..DOT_NiceMath

using LinearAlgebra: Hermitian, I as Id

using Unitful
using Unitful: μs


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.1: Z,X ops
#

N(A ::Int, ℂ::Type{<:Complex}) = begin
    @assert A ≥ 1
    let 𝙽 = [  ℂ(0)      0
                0       +1  ]

        𝙽 ⊗ Id(2^(A-1)) |> Hermitian
    end
end

X(A ::Int, γ::ℂ) where{ℂ} = begin
    @assert A ≥ 1
    let 𝚇 = [   0     γ
                γ'    0   ]

        𝚇 ⊗ Id(2^(A-1)) |> Hermitian
    end
end


# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2. Work horses

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.1: timestep!()
function  timestep!(ψ  ::Vector{ℂ},
                    𝛥𝑡 ::μs_t{ℝ}
                    ;
                    𝜔  ::Rad_per_μs_t{ℝ},
                    X  ::Hermitian{ℂ,𝕄_t},
                    𝛿  ::Rad_per_μs_t{ℝ},
                    Z  ::Hermitian{ℂ,𝕄_t},
                    R  ::Hermitian{ℂ,𝕄_t}            ) ::Nothing      where{ℝ,ℂ,𝕄_t}

    let ωΔt ::ℝ              = ustrip(NoUnits, 𝜔⋅𝛥𝑡),
        δΔt ::ℝ              = ustrip(NoUnits, 𝛿⋅𝛥𝑡),
        Δt  ::ℝ              = ustrip(μs, 𝛥𝑡),
        A   ::Hermitian{ℂ,𝕄_t} = ωΔt⋅X - δΔt⋅Z + Δt⋅R

        ψ .= cis(A)'ψ
    end
    nothing
end #^ timestep!()


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.2: schröd!()
@doc raw"""
Function `schröd!(ψ  ::Vector{ℂ},  𝑇 ::μs_t{ℝ},  γ ::ℂ ; ...)   where{ℝ,ℂ,𝕄_t}`

Simulates time evolution under time-dependent Hamiltonian
```math
    H(t)/\hbar = \omega(t) Xᵧ - \delta(t) Z + R
```

where `Xᵧ` is defined as |1⟩⟨2|⋅γ + |2⟩⟨1|⋅γ̄.  The phase, γ, is not time-dependent.

The number of atoms, `A`, is derived from the dimensions of ``R``, which is assumed to be
a 2ᴬ×2ᴬ matrix.  The matrices ``X,Z`` are constructed to act on the "first" atom, i.e.,
X⊗1⊗...⊗1 and Z⊗1⊗...⊗1; see the function `⊗` above.

## Arguments

Note the convention: *Italic variables are unitful*. That also applies to functions: The
returned value of an italic-type function must be unitful.

#### Place variables
* `ψ` — state vector at time 𝑡=0μs; will be modified in place.
* `𝑇` — end-time of evolution, in μs.

#### Keyword arguments
* `𝜔`  — function of time (in `μs{ℝ}`) , with return value of type `Rad_per_μs_t{ℝ}`
* `𝛿`  — function of time (in `μs{ℝ}`) , with return value of type `Rad_per_μs_t{ℝ}`
* `R` — Hermitian operators, of type `Hermitian{ℂ,𝕄_t}`

#### Returns
...nothing.



"""
function schröd!(ψ  ::Vector{ℂ},
                 𝑇  ::μs_t{ℝ},
                 γ  ::ℂ
                 ;
                 𝜔  ::Function,
                 𝛿  ::Function,
                 R  ::Hermitian{ℂ,𝕄_t}               ) ::Nothing      where{ℝ,ℂ,𝕄_t}




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
