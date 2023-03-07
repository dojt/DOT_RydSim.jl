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
export no_thing


# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 0. Packages

using DOT_NiceMath
using DOT_NiceMath.NumbersF64

using  Unitful: Quantity, μs, 𝐓, Unit, FreeUnits
import Unitful

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1. Types & Units

const μs_t{              REAL<:Real } =                                                             # 1→ `μs_t`
    Quantity{REAL, 𝐓,
             FreeUnits{ (Unit{:Second,𝐓}(-6, 1//1),),
                        𝐓,
                        nothing }
             }

const Rad_per_μs_t{      REAL<:Real } =                                                             # 1→ `Rad_per_μs_t`
    Quantity{REAL, 𝐓^(-1//1),
             FreeUnits{ (Unit{:Second,𝐓}(-6,-1//1),),
                        𝐓^(-1//1),
                        nothing }
             }

const Radperμs_per_μs_t{ REAL<:Real } =                                                             # 1→ `Rad_per_μs_t`
    Quantity{REAL, 𝐓^(-2//1),
             FreeUnits{ (Unit{:Second,𝐓}(-6,-2//1),),
                        𝐓^(-2//1),
                        nothing }
             }


include("Schrödinger.mod.jl")

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2. Function maker

function make_Ctrl_fn(;) ::Function

    for near-square:

    phases/events
t[0]=0
        wait_before
𝑡[1]
        ramp_up
𝑡[2]
        hold
𝑡[3]
        ramp_down
𝑡[4]
        wait_after
𝑡[5]=𝑇


    return: 𝜔 ::Function, γ::ℂ

   end


end # module DOT_RydSim
# EOF
