########################################################################
#                                                                      #
# DOT_RydSim.jl                                                        #
#                                                                      #
# (c) Dirk Oliver Theis 2023                                           #
#                                                                      #
# License:                                                             #
#                                                                      #
#             Apache 2.0                                               #
#                                                                      #
########################################################################

"""
Module `DOT_RydSim`

Quantum simulation of (small!!) arrays of Rydberg atoms.

# Exports

## General
* (nothing yet)

## Number definitions (in sub-modules `Numbers`𝑥𝑦𝑧)
* (also nothing)

# Sub-modules

Sub-module names are not exported.

* `Schrödinger` — Simulation of quantum evolution
"""
module DOT_RydSim
export Fn_Select, schröd!


# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 0. Packages

using DOT_NiceMath
using DOT_NiceMath.NumbersF64

using  Unitful: Quantity, μs, 𝐓, Unit, FreeUnits
import Unitful

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1. Types, Units, Helpers

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.1. Unit types

const μs_t{              𝕂<:Real } =                                                                # 1→ `μs_t`
    Quantity{𝕂, 𝐓,
             FreeUnits{ (Unit{:Second,𝐓}(-6, 1//1),),
                        𝐓,
                        nothing }
             }

const Rad_per_μs_t{      𝕂<:Real } =                                                                # 1→ `Rad_per_μs_t`
    Quantity{𝕂, 𝐓^(-1//1),
             FreeUnits{ (Unit{:Second,𝐓}(-6,-1//1),),
                        𝐓^(-1//1),
                        nothing }
             }

const Radperμs_per_μs_t{ 𝕂<:Real } =                                                                # 1→ `Rad_per_μs_t`
    Quantity{𝕂, 𝐓^(-2//1),
             FreeUnits{ (Unit{:Second,𝐓}(-6,-2//1),),
                        𝐓^(-2//1),
                        nothing }
             }


const GHz_t{             𝕂<:Real } =                                                                # 1→ `Hz_t`
    Quantity{𝕂, 𝐓^(-1//1),
             FreeUnits{ (Unit{:Hertz,𝐓^(-1//1)}(9,1//1),),
                        𝐓^(-1//1),
                        nothing }
             }


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.2. Types: Fn_Select

@doc raw"""
Abstract type `Pulse`

This abstract type establishes an interface for querying properties of pulses.

With ``f`` denoting the pulse shape (as a function of time in μs and with values in radians), syntax
and semantics of the interface are as follows (note the italics function names, in line with the
convention *unitful iff italics*).

Methdos for the functions `𝑎𝑣𝑔()`, `𝑠𝑡𝑒𝑝()`, `phase()`, `plot!()` must be defined:

  * `𝑎𝑣𝑔(p::Pulse,  𝑡 ::μs_t ; 𝛥𝑡 ::μs_t) ::Rad_per_μs_t` — returns

    ```math
    \mu_{t,Δ\!t} := \tfrac{1}{\Delta\!t} \int_t^{t+\Delta\!t} f(s) \,ds
    ```

  * `𝑠𝑡𝑒𝑝(p::Pulse, 𝑡 ::μs_t ; ε ::ℝ) ::μs_t` — returns the largest ``\Delta\!t``
    such that:

    ```math
    \int_t^{t+\Delta!t} |f(s) - \mu_{t,\Delta!t} |\,ds \le \varepsilon
    ```

    with ``\mu_{.,.}`` as above.  The `ε` is not italics as radians is, strictly speaking, not a
    unit.

  * `phase(p::Pulse) ::ℝ` — returns the phase, which must be time-independent.

  * `plot!(plotobject, p::Pulse)` — plots the pulse into the given plot object.
"""
abstract type Pulse end

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.3. Sub-module Schrödinger

include("Schrödinger.mod.jl")

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.4. Helper: Rounding

function δround( x ::𝕂₁
                 ;
                 δ ::Rational{ℤ}                    ) ::Rational{ℤ}                    where{𝕂₁, ℤ}

    δ ⋅ rationalize(ℤ,
                    floor(x/δ +1//2)  )
end

function δround( x ::Quantity{𝕂₁,T₁,F₁}
                 ;
                 δ ::Quantity{Rational{ℤ} ,T₂,F₂}   ) ::Quantity{Rational{ℤ},T₂,F₂}    where{𝕂₁,T₁,F₁, ℤ,T₂,F₂}

    δ ⋅ rationalize(ℤ,
                    floor(x/δ +1//2)  )
end

# Do I need any of these:
#
# _ratdiv(num,den) = rat( NoUnits( div( num, den , RoundNearest) ) )
# _round(x::Frequency ; δ::Frequency) = MHz(δ)⋅_ratdiv( MHz(x),MHz(δ) )
# _round(x::Time      ; δ::Time     ) =  μs(δ)⋅_ratdiv(  μs(x), μs(δ) )
# _round(x::Length    ; δ::Length   ) =  μm(δ)⋅_ratdiv(  μm(x), μm(δ) )
#
# round_Ω(dev::HW_Descr,  𝛺::Frequency) = _round(𝛺;δ=dev.𝛺ᵣₑₛ)
# round_Δ(dev::HW_Descr,  𝛥::Frequency) = _round(𝛥;δ=dev.𝛥ᵣₑₛ)
# round_t(dev::HW_Descr,  𝑡::Time     ) = _round(𝑡;δ=dev.𝑡ᵣₑₛ)
#
# round_xy(dev::HW_Descr, 𝑧::Length   ) = _round(𝑧;δ=dev.lattice.posᵣₑₛ)
#
# φπ(dev::HW_Descr) = dev.𝜑ᵣₑₛ⋅_ratdiv(π,dev.𝜑ᵣₑₛ)


# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2. Pulse constructors


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.2. Ω_BangBang


@doc raw"""
* wait before pulse
* event 1
* ramp up
* event 2
* plateau
* event 3
* ramp down
* event 4
* wait after pulse
* end of time
"""
struct Pulse__Ω_BangBang{ℚ,ℝ} <: Pulse
    γ      ::Complex{ℝ}
    𝑒𝑣𝑒𝑛𝑡𝑠 ::NTuple{5,ℚ}
end

function Pulse__Ω_BangBang{ℚ,ℝ}( 𝑡₀       ::μs_t{ℚ},
                                 𝑡₁       ::μs_t{ℚ},
                                 𝑇        ::μs_t{ℚ},
                                 𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 ::Rad_per_μs_t{ℚ}
                                 ;
                                 𝛺ₘₐₓ       ::Rad_per_μs_t{ℚ},
                                 𝛺ᵣₑₛ       ::Rad_per_μs_t{ℚ},
                                 𝛺_𝑚𝑎𝑥_𝑠𝑙𝑒𝑤 ::Radperμs_per_μs_t{ℚ},

                                 φᵣₑₛ       ::ℚ,

                                 𝑡ₘₐₓ       ::μs_t{ℚ},
                                 𝑡ᵣₑₛ       ::μs_t{ℚ},
                                 𝛥𝑡ₘᵢₙ      ::μs_t{ℚ}                ) ::Pulse__Ω_BangBang{ℚ,ℝ}

    ℂ = Complex{ℝ}

    @assert -𝛺ₘₐₓ < 𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 < +𝛺ₘₐₓ
    @assert 0μs ≤ 𝑡₀ < 𝑡₁ ≤ 𝑇

    γ ::ℂ = cis( δround(ℝ(π);δ=ϕᵣₑₛ) )

    𝑒𝑣𝑒𝑛𝑡𝑠 ::NTuple{5,ℚ} =
        (
                #   wait before pulse
            𝑡₀, #1              ⌝
                #   ramp up     |
            1//1,#2             |
                #   plateau     |  pulse, incl. ramp-down
            1//1,#3             |
                #   ramp down   |
            𝑡₁, #4              ⌟
                #   wait after pulse
            𝑇   #5
        )


    return Pulse__Ω_BangBang(γ,𝑒𝑣𝑒𝑛𝑡𝑠)
end

(  phase(Ω::Pulse__Ω_BangBang{ℚ,ℂ}) ::ℂ  ) where{ℚ,ℂ}        = Ω.γ

function 𝑎𝑣𝑔(Ω::Pulse__Ω_BangBang{ℚ,ℂ},  𝑡 ::μs_t{𝕂} ;  𝛥𝑡 ::μs_t{𝕂}) ::𝕂   where{ℚ,ℂ,𝕂}
    
end

function 𝑠𝑡𝑒𝑝(Ω::Pulse__Ω_BangBang{ℚ,ℂ}, 𝑡 ::μs_t{𝕂} ;  ε  ::μs_t{𝕂}) ::𝕂   where{ℚ,ℂ,𝕂}
    blubb
end


#function plot!(plothere, Ω::Pulse__Ω_BangBang ; kwargs...)
#    Plots.plot!(plothere, [ 𝑡 for 𝑡 ∈ Ω.𝑒𝑣𝑒𝑛𝑡𝑠 ], ; kwargs...)
#end


#     for near-square:
#
#     phases/events
# t[0]=0
#         wait_before
# 𝑡[1]
#         ramp_up
# 𝑡[2]
#         hold
# 𝑡[3]
#         ramp_down
# 𝑡[4]
#         wait_after
# 𝑡[5]=𝑇
#
#
#     return: 𝜔 ::Function, γ::ℂ


   end


end # module DOT_RydSim
# EOF
