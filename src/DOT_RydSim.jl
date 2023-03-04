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

const μs_t{         REAL<:Real } = Quantity{REAL, 𝐓,                                                 # 1→ `μs_t`
                                            FreeUnits{ (Unit{:Second,𝐓}(-6, 1//1),),
                                                       𝐓,
                                                       nothing }
                                            }

const Rad_per_μs_t{ REAL<:Real } = Quantity{REAL, 𝐓^(-1//1),                                         # 1→ `Rad_per_μs_t`
                                            FreeUnits{ (Unit{:Second,𝐓}(-6,-1//1),),
                                                       𝐓^(-1//1),
                                                       nothing }
                                            }



include("Schrödinger.mod.jl")






# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2. Shots & Shooting

struct Halfint                                                                                      # 2→ struct `Halfint`
end

end # module DOT_RydSim
# EOF
