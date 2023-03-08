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

const μs_t{              ℝ<:Real } =                                                                # 1→ `μs_t`
    Quantity{ℝ, 𝐓,
             FreeUnits{ (Unit{:Second,𝐓}(-6, 1//1),),
                        𝐓,
                        nothing }
             }

const Rad_per_μs_t{      ℝ<:Real } =                                                                # 1→ `Rad_per_μs_t`
    Quantity{ℝ, 𝐓^(-1//1),
             FreeUnits{ (Unit{:Second,𝐓}(-6,-1//1),),
                        𝐓^(-1//1),
                        nothing }
             }

const Radperμs_per_μs_t{ ℝ<:Real } =                                                                # 1→ `Rad_per_μs_t`
    Quantity{ℝ, 𝐓^(-2//1),
             FreeUnits{ (Unit{:Second,𝐓}(-6,-2//1),),
                        𝐓^(-2//1),
                        nothing }
             }


const GHz_t{             ℝ<:Real } =                                                                # 1→ `Hz_t`
    Quantity{ℝ, 𝐓^(-1//1),
             FreeUnits{ (Unit{:Hertz,𝐓^(-1//1)}(9,1//1),),
                        𝐓^(-1//1),
                        nothing }
             }


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.2. Types: Fn_Select

"""
Module `Fn_Select`

This sub-module of `DOT_RydSim` provides the dummy types
* `STEP` and
* `AVG`
to select the method for the pulse-shape functions.  Used in the Function Maker functions below.

"""
module Fn_Select
    export STEP, AVG

    abstract type Function_Selection end

    struct STEP <: Function_Selection end
    struct AVG  <: Function_Selection end

end

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.3. Sub-module Schrödinger

include("Schrödinger.mod.jl")

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.4. Helper: Rounding

using DOT_NiceMath.NumbersF64: ℤ, ℚ

continue here: Copy round fn from Pluto notebook

......................................................................................



# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2. Function maker


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.2. Ω_BangBang

function make_ctrl_fn__Ω_bangbang( 𝑡₀       ::μs_t{ℚ},
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
                                   𝛥𝑡ₘᵢₙ      ::μs_t{ℚ}                  ) ::Function

    @assert -𝛺ₘₐₓ < 𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 < +𝛺ₘₐₓ
    @assert 0μs ≤ 𝑡₀ < 𝑡₁ ≤ 𝑇

    𝑒𝑣𝑒𝑛𝑡𝑠 ::NTuple{5,ℚ} =
        (
                #   wait before pulse
            𝑡₀, #1              ⌝
                #   ramp up     |
            ,   #2              |
                #   plateau     |  pulse, incl. ramp-down
            ,   #3              |
                #   ramp down   |
            𝑡₁  #4              ⌟
                #   wait after pulse
            𝑇   #5
        )


    function
    𝜔(::Type{Fn_Select.AVG},  𝑡 ::μs_t{𝐑} ;  𝛥𝑡 ::μs_t{𝐑}) ::𝐑   where{𝐑<:Real}
        blah
    end
    function
    𝜔(::Type{Fn_Select.Step}, 𝑡 ::μs_t{𝐑} ;  ε  ::μs_t{𝐑}) ::𝐑   where{𝐑<:Real}
        blubb
    end


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
