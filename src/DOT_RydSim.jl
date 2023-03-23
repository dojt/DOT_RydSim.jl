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
export schröd!
export Pulse__Ω_BangBang
export plotpulse


# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 0. Packages

using Base: @kwdef       # remove once Julia 1.9 comes out

using DOT_NiceMath

using  Unitful: Quantity, μs, 𝐓, Unit, FreeUnits
import Unitful

using Logging

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1. Types, Units, Helpers

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.1. Unit types

const μs_t{              𝕂<:Real } =                                                                #(1.1) `μs_t`
    Quantity{𝕂, 𝐓,
             FreeUnits{ (Unit{:Second,𝐓}(-6, 1//1),),
                        𝐓,
                        nothing }
             }

const Rad_per_μs_t{      𝕂<:Real } =                                                                #(1.1) `Rad_per_μs_t`
    Quantity{𝕂, 𝐓^(-1//1),
             FreeUnits{ (Unit{:Second,𝐓}(-6,-1//1),),
                        𝐓^(-1//1),
                        nothing }
             }

const Radperμs_per_μs_t{ 𝕂<:Real } =                                                                #(1.1) `Rad_per_μs_t`
    Quantity{𝕂, 𝐓^(-2//1),
             FreeUnits{ (Unit{:Second,𝐓}(-6,-2//1),),
                        𝐓^(-2//1),
                        nothing }
             }


const GHz_t{             𝕂<:Real } =                                                                #(1.1) `Hz_t`
    Quantity{𝕂, 𝐓^(-1//1),
             FreeUnits{ (Unit{:Hertz,𝐓^(-1//1)}(9,1//1),),
                        𝐓^(-1//1),
                        nothing }
             }



# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.3. Helper: Unitful bug fix
# Fix for bug in Unitful
import Base.==
import Base.:≤
import Base.:<
(   ( x::μs_t{𝕂₁} == y::μs_t{𝕂₂} ) ::Bool   ) where{𝕂₁,𝕂₂}      = x.val == y.val
(   ( x::μs_t{𝕂₁} ≤  y::μs_t{𝕂₂} ) ::Bool   ) where{𝕂₁,𝕂₂}      = x.val ≤  y.val
(   ( x::μs_t{𝕂₁} <  y::μs_t{𝕂₂} ) ::Bool   ) where{𝕂₁,𝕂₂}      = x.val <  y.val


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

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2. Pulse constructors


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.1. Types: Pulse base type

@doc raw"""
Abstract type `Pulse`

This abstract type establishes an interface for querying properties of pulses.

With ``f`` denoting the pulse shape (as a function of time in μs and with values in radians per μs
aka MHz/2π), syntax and semantics of the interface are as follows (note the italics function names, in
line with the convention *unitful iff italics*).

Methdos for the functions `𝑎𝑣𝑔()`, `𝑠𝑡𝑒𝑝()`, `phase()`, `plotpulse()` must be defined:

  * `𝑎𝑣𝑔(p::Pulse,  𝑡 ::μs_t ; 𝛥𝑡 ::μs_t) ::Rad_per_μs_t` — returns

    ```math
    \mu_{t,Δ\!t} := \tfrac{1}{\Delta\!t} \int_t^{t+\Delta\!t} f(s) \,ds
    ```

  * `𝑠𝑡𝑒𝑝(p::Pulse, 𝑡 ::μs_t ; ε ::ℝ) ::μs_t` — returns the largest ``\Delta\!t``
    such that:

    ```math
    \int_t^{t+\Delta!t} |f(s) - \mu_{t,\Delta!t} |\,ds \le \varepsilon
    ```

    with ``\mu_{.,.}`` as above.  (`ε` is `\varepsilon`.)

    !!! note "Note."
        To simplify implementation, returning a *lower bound* on that maximum is considerd conformant
        with the interface.

  * `phase(p::Pulse) ::ℝ` — returns the phase, which must be time-independent.

  * `plotpulse(p::Pulse) :: @NamedTuple{x⃗::Vector,y⃗::Vector}` — returns x- and y-data for plotting.
"""
abstract type Pulse end

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.2. Ω_BangBang Pulse

@doc raw"""
Struct `Pulse__Ω_BangBang` `<:` `Pulse`

## Constructor
```julia
Pulse__Ω_BangBang{ℚ,ℝ}( 𝑡₀       ::μs_t{ℚ},
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
```

## Implementation

### The struct `Pulse__Ω_BangBang{ℚ,ℝ}`

#### Semantics of `𝑒𝑣`
The tuple `𝑒𝑣` holds times of events between phases:
  * 0μs   event: beginning of time
  * —     phase: wait before pulse
  * `𝑒𝑣[1]`
  * —     phase: ramp up
  * `𝑒𝑣[2]`
  * —     phase: plateau
  * `𝑒𝑣[3]`
  * —     phase: ramp down
  * `𝑒𝑣[4]`
  * —     phase: wait after pulse
  * `𝑒𝑣[5]` event: end of time

Implied in this: Entries are increasing with index.

#### Docs of other fields:
See source!
"""
struct Pulse__Ω_BangBang{ℚ,ℝ} <: Pulse                                                              #(2.2) struct Pulse__Ω_BangBang
    γ    ::Complex{ℝ}                   # phase
    𝑒𝑣   ::NTuple{5, μs_t{ℚ} }          # events
    𝑟ꜛ   ::Radperμs_per_μs_t{ℚ}         # up-ramp rate
    𝛺    ::Rad_per_μs_t{ℚ}              # top plateau value
    𝑟ꜜ   ::Radperμs_per_μs_t{ℚ}         # down-ramp rate
end

function Pulse__Ω_BangBang{ℚ,ℝ}(𝑡ᵒⁿ      ::μs_t{ℚ},                                                 #(2.2) constructor Pulse__Ω_BangBang
                                𝑡ᵒᶠᶠ     ::μs_t{ℚ},
                                𝑇        ::μs_t{ℚ},
                                𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 ::Rad_per_μs_t{ℚ}
                                ;
                                𝛺ₘₐₓ       ::Rad_per_μs_t{ℚ},
                                𝛺ᵣₑₛ       ::Rad_per_μs_t{ℚ},
                                𝛺_𝑚𝑎𝑥_𝑠𝑙𝑒𝑤 ::Radperμs_per_μs_t{ℚ},

                                φᵣₑₛ       ::ℚ,                     # "\varphi"

                                𝑡ₘₐₓ       ::μs_t{ℚ},
                                𝑡ᵣₑₛ       ::μs_t{ℚ},
                                𝛥𝑡ₘᵢₙ      ::μs_t{ℚ}                ) ::
                                                          Pulse__Ω_BangBang{ℚ,ℝ}   where{ℚ,ℝ}

    ℂ = Complex{ℝ}

    @assert 
    @assert 0μs ≤ 𝑡ᵒⁿ < 𝑡ᵒᶠᶠ ≤ 𝑇

    𝛺ₘₐₓ > 0/μs              || throw(ArgumentError("𝛺ₘₐₓ must be positive."))
    𝛺_𝑚𝑎𝑥_𝑠𝑙𝑒𝑤 > 0/μs^2      || throw(ArgumentError("Max slew rate 𝛺_𝑚𝑎𝑥_𝑠𝑙𝑒𝑤 must \
                                                     be positive."))
    𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 % 𝛺ᵣₑₛ == 0/μs  || throw(ArgumentError("𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 ($(𝛺_𝑡𝑎𝑟𝑔𝑒𝑡)) is not integer \
                                                     multiple of 𝛺ᵣₑₛ ($(𝛺ᵣₑₛ))."))
    -𝛺ₘₐₓ ≤ 𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 ≤ +𝛺ₘₐₓ || throw(ArgumentError("𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 ($(𝛺_𝑡𝑎𝑟𝑔𝑒𝑡)) is not in \
                                                     range [-𝛺ₘₐₓ,+𝛺ₘₐₓ] ($(𝛺ₘₐₓ))."))
    𝑡ᵒⁿ + 𝛥𝑡ₘᵢₙ ≤ 𝑡ᵒᶠᶠ       || throw(ArgumentError("Gap 𝑡ᵒⁿ → 𝑡ᵒᶠᶠ ($(𝑡ᵒᶠᶠ-𝑡ᵒⁿ)) \
                                                     smaller than 𝛥𝑡ₘᵢₙ ($(𝛥𝑡ₘᵢₙ))."))
    𝑡ᵒⁿ ≤ 0μs || 𝑡ᵒⁿ > 𝛥𝑡ₘᵢₙ || throw(ArgumentError("Gap 0μs → 𝑡ᵒⁿ ($(𝑡ᵒⁿ)) \
                                                     smaller than 𝛥𝑡ₘᵢₙ ($(𝛥𝑡ₘᵢₙ))."))
    𝑡ᵒⁿ %  𝑡ᵣₑₛ == 0μs       || throw(ArgumentError("𝑡ᵒⁿ ($(𝑡ᵒⁿ)) is not integer multiple \
                                                     of 𝑡ᵣₑₛ ($(𝑡ᵣₑₛ)$."))
    𝑡ᵒᶠᶠ % 𝑡ᵣₑₛ == 0μs       || throw(ArgumentError("𝑡ᵒᶠᶠ ($(𝑡ᵒᶠᶠ)) is not integer multiple \
                                                     of 𝑡ᵣₑₛ ($(𝑡ᵣₑₛ))."))

    γ::ℂ =
        if 𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 < 0/μs
            𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 = -𝛺_𝑡𝑎𝑟𝑔𝑒𝑡
            cis( δround(ℝ(π);δ=φᵣₑₛ) )
        else
            ℂ(0)
        end


    𝑟ꜛ     = 𝛺_𝑚𝑎𝑥_𝑠𝑙𝑒𝑤
    𝑟ꜜ     = 𝛺_𝑚𝑎𝑥_𝑠𝑙𝑒𝑤
    𝑡₀₋ₜₐᵣ = 𝛺_𝑡𝑎𝑟𝑔𝑒𝑡/ 𝑟ꜛ
    𝑡ₜₐᵣ₋₀ = 𝛺_𝑡𝑎𝑟𝑔𝑒𝑡/ 𝑟ꜜ

    𝑒𝑣 ::NTuple{5, μs_t{ℚ} } = if   𝑡₀+𝑡₀₋ₜₐᵣ  ≤  𝑡₁-𝑡ₜₐᵣ₋₀
        (
            # wait before pulse
            𝑡₀,
            # ramp up
            𝑡₀+𝑡₀₋ₜₐᵣ,
            # plateau
            𝑡₁-𝑡ₜₐᵣ₋₀,
            #   ramp down
            𝑡₁,
            # wait after pulse
            𝑇
        )
    else let 𝛥𝑡 = 𝑡₁ - 𝑡₀
        # Solve  𝑠⋅𝑟ꜛ = (𝛥𝑡-𝑠)⋅𝑟ꜜ   for 𝑠:
        𝑠 = 𝑟ꜜ/( 𝑟ꜛ+𝑟ꜜ )⋅𝛥𝑡
        (
            # wait before pulse
            𝑡₀,
            # ramp up
            𝑡₀+𝑠,
            # plateau is empty!!
            𝑡₁-(𝛥𝑡-𝑠),
            #   ramp down
            𝑡₁,
            # wait after pulse
            𝑇
        )
    end end

    𝑒𝑣[4] ≤ 𝑒𝑣[5]    || throw(ArgumentError("𝛺_BangBang pulse shape doesn't fit: \
                                             gap between 𝑡₁=$(𝑡₁) and 𝑇=$(𝑇) \
                                             too small for slew rate $(𝛺_𝑚𝑎𝑥_𝑠𝑙𝑒𝑤)."))

    return Pulse__Ω_BangBang(γ, 𝑒𝑣, 𝑟ꜛ, 𝛺_𝑡𝑎𝑟𝑔𝑒𝑡, 𝑟ꜜ)
end

function phase(Ω::Pulse__Ω_BangBang{ℚ,ℝ}) ::Complex{ℝ}      where{ℚ,ℝ}                              #(2.2) phase() Pulse__Ω_BangBang
    # let's take the opportunity to run some checks

    issorted(Ω.𝑒𝑣)   ||  throw(ErrorException("Pulse__Ω_BangBang: \
                                               𝑒𝑣=$(Ω.𝑒𝑣) not sorted. This is a bug."))y
    Ω.𝛺 ≥ 0/μs       ||  throw(ErrorException("Pulse__Ω_BangBang: \
                                               negative 𝛺==$(Ω.𝛺). This is a bug."))

    return Ω.γ
end

#
# This function is to demonstrate the pulse shape data, and maybe for plotting or whatnot.
#
function (Ω::Pulse__Ω_BangBang{ℚ,ℝ})(𝑡 ::μs_t{𝕂}) ::Rad_per_μs_t{𝕂}   where{ℚ,ℝ,𝕂}                  #(2.2) callable Pulse__Ω_BangBang

    (;𝑒𝑣, 𝑟ꜛ, 𝑟ꜜ, 𝛺) = Ω

    β = (2^30+1)//2^30
    if            𝑡 < 0μs            throw(DomainError(𝑡,"Time cannot be negative."))
    elseif  0μs   ≤ 𝑡 ≤ 𝑒𝑣[1]        return 𝕂(0)/μs
    elseif  𝑒𝑣[1] ≤ 𝑡 ≤ 𝑒𝑣[2]        return ( 𝑡-𝑒𝑣[1] )⋅𝑟ꜛ
    elseif  𝑒𝑣[2] ≤ 𝑡 ≤ 𝑒𝑣[3]        return 𝛺
    elseif  𝑒𝑣[3] ≤ 𝑡 ≤ 𝑒𝑣[4]        return ( 𝑒𝑣[4]-𝑡 )⋅𝑟ꜜ
    elseif  𝑒𝑣[4] ≤ 𝑡 ≤ 𝑒𝑣[5]⋅β      return 𝕂(0)/μs
    else                             throw(DomainError(𝑡,"Time exceeds upper bound, 𝑇=$(𝑒𝑣[5])."))
    end
end #^ callable Pulse__Ω_BangBang

function 𝑎𝑣𝑔(Ω ::Pulse__Ω_BangBang{ℚ,ℝ},                                                            #(2.2) 𝑎𝑣𝑔() Pulse__Ω_BangBang
             𝑡 ::μs_t{𝕂}
             ;
             𝛥𝑡 ::μs_t{𝕂}               ) ::Rad_per_μs_t{𝕂}       where{ℚ,ℝ,𝕂}

    (;𝑒𝑣) = Ω
    𝑡ᵉⁿᵈ  = 𝑡+𝛥𝑡
    sum   = 𝕂(0)
    for j = 1 : length(𝑒𝑣)-1
        if 𝑡 < 𝑒𝑣[j+1] && 𝑒𝑣[j] < 𝑡ᵉⁿᵈ
            𝑠ⱼ = max(𝑒𝑣[j], 𝑡)
            𝑡ⱼ = min(𝑡ᵉⁿᵈ, 𝑒𝑣[j+1])
            if 𝑠ⱼ < 𝑡ⱼ
                𝜔ₛ = Ω(𝑠ⱼ)
                𝜔ₜ = Ω(𝑡ⱼ)
                sum += (𝑡ⱼ-𝑠ⱼ)⋅( 𝜔ₛ + (𝜔ₜ-𝜔ₛ)/2 )
            end
        end
    end
    return sum/𝛥𝑡
end #^ 𝑎𝑣𝑔()

function 𝑠𝑡𝑒𝑝(Ω::Pulse__Ω_BangBang{ℚ,ℝ},                                                            #(2.2) 𝑠𝑡𝑒𝑝() Pulse__Ω_BangBang
              𝑡 ::μs_t{𝕂}
              ;
              ε  ::𝕂                    ) ::μs_t{𝕂}   where{ℚ,ℝ,𝕂}

    (;𝑒𝑣, 𝑟ꜛ, 𝑟ꜜ) = Ω

    # Lazy: We compare not to the average but to the value

    β = (2^30+1)//2^30
    if            𝑡 < 0μs            throw(DomainError(𝑡,"Time cannot be negative."))
    elseif  0μs   ≤ 𝑡 < 𝑒𝑣[1]        return                 𝑒𝑣[1]-𝑡
    elseif  𝑒𝑣[1] ≤ 𝑡 < 𝑒𝑣[2]        return min( √(2ε/𝑟ꜛ) , 𝑒𝑣[2]-𝑡 )
    elseif  𝑒𝑣[2] ≤ 𝑡 < 𝑒𝑣[3]        return                 𝑒𝑣[3]-𝑡
    elseif  𝑒𝑣[3] ≤ 𝑡 < 𝑒𝑣[4]        return min( √(2ε/𝑟ꜜ) , 𝑒𝑣[4]-𝑡 )
    elseif  𝑒𝑣[4] ≤ 𝑡 ≤ 𝑒𝑣[5]⋅β      return max(            𝑒𝑣[5]-𝑡 , 0μs)
    else                             throw(DomainError(𝑡,"Time exceeds upper bound, 𝑇=$(𝑒𝑣[5])."))
    end

end #^ 𝑠𝑡𝑒𝑝()


function plotpulse(Ω::Pulse__Ω_BangBang) ::NamedTuple                                               #(2.2) plotpulse() Pulse__Ω_BangBang

    𝑋 = Iterators.flatten( [ [(0//1)μs], (𝑡 for 𝑡 ∈ Ω.𝑒𝑣) ] )

    return (  x⃗ = collect(𝑋),
              y⃗ = [ Ω(𝑥) for 𝑥 ∈ 𝑋 ]  )
end

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 3. Sub-module Schrödinger

include("Schrödinger.mod.jl")

end # module DOT_RydSim
# EOF
