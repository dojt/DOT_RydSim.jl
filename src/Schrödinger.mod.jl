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
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 0. Imports

import ..μs_t
import ..Rad_per_μs_t
using  ..DOT_NiceMath

using LinearAlgebra: Hermitian
using Unitful
using Unitful: μs

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.1: timestep()
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
        A   ::Hermitian{ℂ,𝕄_t} = ωΔt⋅X + δΔt⋅Z + Δt⋅R

        ψ .= cis(A)'ψ
    end
    nothing
end #^ timestep!()

function schröd!(ψ  ::Vector{ℂ},
                 𝑇  ::μs_t{ℝ}
                 ;
                 𝜔  ::Function, # with values in Rad_per_μs_t{ℝ},
                 X  ::Hermitian{ℂ,𝕄_t},
                 𝛿  ::Function, # with values in Rad_per_μs_t{ℝ},
                 Z  ::Hermitian{ℂ,𝕄_t},
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
