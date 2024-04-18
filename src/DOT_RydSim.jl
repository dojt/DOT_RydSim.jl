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

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 0. ToC  +
#                                                                                                             |
#  Table of Contents                                                                                          |
#  -----------------                                                                                          |
#                                                                                                             |
#    1.  Module header & imports                                                                              |
#                                                                                                             |
#        1.1.  Exports                                                                                        |
#        1.2.  Imports                                                                                         |
#                                                                                                             |
#                                                                                                             |
#    2.  Types, Units, Helpers                                                                                |
#                                                                                                             |
#        2.1.  Unit types                                                                                     |
#        2.2.  -%-                                                                                            |
#        2.3.  Helper: Unitful bug fix                                                                        |
#        2.4.  Helper: Rounding                                                                               |
#                                                                                                             |
#                                                                                                             |
#    3.  Pulse construction                                                                                   |
#                                                                                                             |
#        3.1.  Pulse base type & interface description                                                        |
#        3.2.  Δ_BangBang Pulse                                                                               |
#        3.3.  Ω_BangBang Pulse                                                                               |
#                                                                                                             |
#                                                                                                             |
#    4.  Sub-module `Schrödinger` (include)                                                                   |
#                                                                                                             |
#    5.  Sub-module `HW_Descriptions` (include)                                                               |
#                                                                                                             |
#—————————————————————————————————————————————————————————————————————————————————————————————————————————————+


# ******************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1. Module header & imports

"""
Module `DOT_RydSim`

Quantum simulation of (small!!) arrays of Rydberg atoms.

# Exports

  * Function [`schröd!`](@ref)`()`

  * Abstract type [`Pulse`](@ref), with sub-types

    * [`Pulse__Ω_BangBang`](@ref),
    * [`Pulse__Δ_BangBang`](@ref)
    * ... (tbc)

  * Functions for using Pulses: [`phase`](@ref), [`𝑎𝑣𝑔`](@ref), [`𝑠𝑡𝑒𝑝`](@ref), [`plotpulse`](@ref)

  * Helper function [`δround`](@ref) and friends, incl. [`is_δrounded`](@ref)`()`

  * Helper functions
     - `(    𝔑(a::Integer, A ::Integer, ␣::ℂ)  ::Hermitian{ℂ,Matrix{ℂ}}    )where{ℂ}`   and
     - `(    𝔛(a::Integer, A ::Integer, γ::ℂ)  ::Hermitian{ℂ,Matrix{ℂ}}    )where{ℂ}`
    These functions return the Rabi term operator γ |g⟩⟨r|ₐ + γ̄ |r⟩⟨g|ₐ and dephasing term operator
    |r⟩⟨r|ₐ for atom `a` of the Rydberg Hamiltonian.  The 3rd argument (complex number) determines
    the type `ℂ`: its value is ignored in `𝔑`, and gives the frame in `𝔛`.

# Sub-modules

Sub-module names are not exported.

* [`Schrödinger`](@ref) — Simulation of quantum evolution
* [`HW_Descriptions`](@ref) — Types for defining and functions for reading hardware
  descriptions.
"""
module DOT_RydSim


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.1. Exports

export schröd!
export Pulse, phase, 𝑎𝑣𝑔, 𝑠𝑡𝑒𝑝, plotpulse
export δround, δround_down, δround_up, δround_to0,  is_δrounded
export Pulse__Ω_BangBang, Pulse__Δ_BangBang
export 𝔑, 𝔛



# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.2. Imports

using DOT_NiceMath

using  Unitful: Quantity, μs, 𝐓, Unit, FreeUnits
import Unitful

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2. Types, Units, Helpers

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.1. Unit types

const μs_t{              𝕂<:Real } =                                                                #(2.1) `μs_t`
    Quantity{𝕂, 𝐓,
             FreeUnits{ (Unit{:Second,𝐓}(-6, 1//1),),
                        𝐓,
                        nothing }
             }

const Rad_per_μs_t{      𝕂<:Real } =                                                                #(2.1) `Rad_per_μs_t`
    Quantity{𝕂, 𝐓^(-1//1),
             FreeUnits{ (Unit{:Second,𝐓}(-6,-1//1),),
                        𝐓^(-1//1),
                        nothing }
             }

const Radperμs_per_μs_t{ 𝕂<:Real } =                                                                #(2.1) `Rad_per_μs_t`
    Quantity{𝕂, 𝐓^(-2//1),
             FreeUnits{ (Unit{:Second,𝐓}(-6,-2//1),),
                        𝐓^(-2//1),
                        nothing }
             }


const GHz_t{             𝕂<:Real } =                                                                #(2.1) `Hz_t`
    Quantity{𝕂, 𝐓^(-1//1),
             FreeUnits{ (Unit{:Hertz,𝐓^(-1//1)}(9,1//1),),
                        𝐓^(-1//1),
                        nothing }
             }



# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.3. Helper: Unitful bug fix
# Fix for bug in Unitful
import Base.==
import Base.:≤
import Base.:<
(  ( x::μs_t{𝕂₁} == y::μs_t{𝕂₂} ) ::Bool  ) where{𝕂₁               ,𝕂₂          } = x.val == y.val
(  ( x::μs_t{𝕂₁} ≤  y::μs_t{𝕂₂} ) ::Bool  ) where{𝕂₁               ,𝕂₂          } = x.val ≤  y.val
(  ( x::μs_t{𝕂₁} <  y::μs_t{𝕂₂} ) ::Bool  ) where{𝕂₁               ,𝕂₂          } = x.val <  y.val

(  ( x::μs_t{𝕂₁} ≤  y::μs_t{𝕂₂} ) ::Bool  ) where{𝕂₁<:AbstractFloat,𝕂₂<:Rational}      =
    begin
        x.val ⋅ y.val.den  ≤  y.val.num
    end
(  ( x::μs_t{𝕂₁} <  y::μs_t{𝕂₂} ) ::Bool  ) where{𝕂₁<:AbstractFloat,𝕂₂<:Rational}      =
    begin
        x.val ⋅ y.val.den  <  y.val.num
    end

(  ( x::μs_t{𝕂₁} ≤  y::μs_t{𝕂₂} ) ::Bool  ) where{𝕂₁<:Rational     ,𝕂₂<:AbstractFloat} =
    begin
        x.val.num  ≤  x.val.den ⋅ y.val
    end
(  ( x::μs_t{𝕂₁} <  y::μs_t{𝕂₂} ) ::Bool  ) where{𝕂₁<:Rational     ,𝕂₂<:AbstractFloat} =
    begin
        x.val.num  <  x.val.den ⋅ y.val
    end


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.4. Helper: Rounding

import Base: rationalize

rationalize( ::Type{I}, x ::Rational{I} ) where{I<:Integer}    = x

@doc raw"""
Functions
```julia
     δround(     x ::         𝕂      ; δ ::         ℚ     ) ::         ℚ
     δround(     𝑥 ::Quantity{𝕂,...} ; 𝛿 ::Quantity{ℚ,...}) ::Quantity{ℚ,...}
     δround_down(𝑥 ::Quantity{𝕂,...} ; 𝛿 ::Quantity{ℚ,...}) ::Quantity{ℚ,...}
     δround_up(  𝑥 ::Quantity{𝕂,...} ; 𝛿 ::Quantity{ℚ,...}) ::Quantity{ℚ,...}
```

Rounds `x` (`𝑥`, resp.) to a multiple of `δ` (`𝛿`, resp.).
"""
function δround( x ::𝕂₁
                 ;
                 δ ::Rational{ℤ}                    ) ::Rational{ℤ}     where{𝕂₁,ℤ}

    δ ⋅ rationalize(ℤ,
                    floor(x/δ +1//2)  )
end

function δround( 𝑥 ::Quantity{𝕂,T₁,F₁}
                 ;
                 𝛿 ::Quantity{Rational{ℤ},T₂,F₂}   ) ::
                                         Quantity{Rational{ℤ},T₂,F₂}    where{𝕂,T₁,F₁, ℤ,T₂,F₂}

    𝛿 ⋅ rationalize(ℤ,     floor(𝑥/𝛿 +1//2)    )
end


@doc raw"""
Functions
```julia
     δround_down(𝑥 ::Quantity{𝕂,...} ; 𝛿 ::Quantity{ℚ,...}) ::Quantity{ℚ,...}
```

Rounds `𝑥` down ("floor") to the closest multiple of `𝛿`.
"""
function δround_down( 𝑥 ::Quantity{𝕂,T₁,F₁}
                      ;
                      𝛿 ::Quantity{Rational{ℤ},T₂,F₂}   ) ::
                                         Quantity{Rational{ℤ},T₂,F₂}    where{𝕂,T₁,F₁, ℤ,T₂,F₂}

    𝛿 ⋅ rationalize(ℤ,     floor(𝑥/𝛿)          )
end

@doc raw"""
Functions
```julia
     δround_down(𝑥 ::Quantity{𝕂,...} ; 𝛿 ::Quantity{ℚ,...}) ::Quantity{ℚ,...}
```

Rounds `𝑥` towards zero ("trunc") to the closest multiple of `𝛿`.
"""
function δround_to0( 𝑥 ::Quantity{𝕂,T₁,F₁}
                     ;
                     𝛿 ::Quantity{Rational{ℤ},T₂,F₂}   ) ::
                                         Quantity{Rational{ℤ},T₂,F₂}    where{𝕂,T₁,F₁, ℤ,T₂,F₂}

    𝛿 ⋅ rationalize(ℤ,     trunc(𝑥/𝛿)          )
end

@doc raw"""
Functions
```julia
     δround_up(𝑥 ::Quantity{𝕂,...} ; 𝛿 ::Quantity{ℚ,...}) ::Quantity{ℚ,...}
```

Rounds `𝑥` up ("ceil") to the closest multiple of `𝛿`.
"""
function δround_up( 𝑥 ::Quantity{𝕂,T₁,F₁}
                    ;
                    𝛿 ::Quantity{Rational{ℤ},T₂,F₂}   ) ::
                                         Quantity{Rational{ℤ},T₂,F₂}    where{𝕂,T₁,F₁, ℤ,T₂,F₂}

    𝛿 ⋅ rationalize(ℤ,     ceil(𝑥/𝛿)           )
end


@doc raw"""
Functions
```julia
    is_δrounded(𝑥 ::Quantity{𝕂,...} ; 𝛿 ::Quantity{ℚ,...}) ::Bool
```

Returns `true` iff `𝑥` is an integer multiple of `𝛿`.
"""
function is_δrounded( 𝑥 ::Quantity{Rational{ℤ},T,F₁}
                      ;
                      𝛿 ::Quantity{Rational{ℤ},T,F₂} ) :: Bool      where{ℤ, T, F₁,F₂}
    @assert 𝛿 ≠ 0 "Nice try."
    return isinteger( 𝑥/𝛿 )
end

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 3. Pulse constructors


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 3.1. Pulse base type & interface description

@doc raw"""
Abstract type `Pulse`

This abstract type establishes an interface for querying properties of pulses.

With ``f`` denoting the pulse shape (as a function of time in μs and with values in radians per μs
aka MHz/2π), syntax and semantics of the interface are as follows (note the italics function names, in
line with the convention *unitful iff italics*).

Methdos for the functions [`𝑎𝑣𝑔`](@ref)`()`, [`𝑠𝑡𝑒𝑝`](@ref)`()`, [`phase`](@ref)`()`, and
[`plotpulse`](@ref)`()` must be defined; see their docs.
"""
abstract type Pulse end


@doc raw"""
Function `phase(p::Pulse) ::ℂ`

Returns the phase (which is be time-independent).

Only Ω (Rabi frequency) pulse shapes need / allow a phase: The Rabi frequency technically
cannot be negative, but a phase of ``e^{i\pi}`` allows "virtually" negative values.  Detuning
(i.e., Δ) pulse shapes don't have a phase, as, unlike the Rabi frequency, the detuning can be
negative.
"""
function phase end

@doc raw"""
Function `𝑎𝑣𝑔(p::Pulse,  𝑡 ::μs_t ; 𝛥𝑡 ::μs_t) ::Rad_per_μs_t`

Returns

```math
\mu_{t,Δ\!t} := \tfrac{1}{\Delta\!t} \int_t^{t+\Delta\!t} f(s) \,ds
```
"""
function 𝑎𝑣𝑔 end

@doc raw"""
`𝑠𝑡𝑒𝑝(p::Pulse, 𝑡 ::μs_t ; ε ::ℝ) ::μs_t` — returns the largest ``\Delta\!t``
such that:

```math
\int_t^{t+\Delta!t} |f(s) - \mu_{t,\Delta!t} |\,ds \le \varepsilon
```

with ``\mu_{.,.}`` as in the docs for [`𝑎𝑣𝑔`](@ref)`()`.  (`ε` is `\varepsilon`.)

!!! note "Note."
    To simplify implementation, returning a (non-trivial) *lower bound* on that maximum is
    considerd conformant with the interface.
"""
function 𝑠𝑡𝑒𝑝 end

@doc raw"""
Function `plotpulse(p::Pulse) :: @NamedTuple{x⃗::Vector,y⃗::Vector}`

Returns x- and y-data for plotting.
"""
function plotpulse end

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 3.2. Δ_BangBang Pulse

@doc raw"""
Struct `Pulse__Δ_BangBang` `<:` `Pulse`

!!! note "Note!"

    Detuning (i.e., "Δ") pulse shapes don't have a phase: The phase is used only for fixing the
    sign of the pulse shape, but unlike the Rabi frequency, the detuning can be negative.

## Constructor
```julia
Pulse__Δ_BangBang{ℚ}( 𝑡ᵒⁿ      ::μs_t{ℚ},
                      𝑡ᵒᶠᶠ     ::μs_t{ℚ},
                      𝑇        ::μs_t{ℚ},
                      𝛥_𝑡𝑎𝑟𝑔𝑒𝑡 ::Rad_per_μs_t{ℚ}
                      ;
                      𝛥ₘₐₓ           ::Rad_per_μs_t{ℚ},
                      𝛥ᵣₑₛ           ::Rad_per_μs_t{ℚ},
                      𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤   ::Radperμs_per_μs_t{ℚ},
                      𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 ::Radperμs_per_μs_t{ℚ} = 𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤,

                      𝑡ₘₐₓ           ::μs_t{ℚ},
                      𝑡ᵣₑₛ           ::μs_t{ℚ},
                      𝛥𝑡ₘᵢₙ          ::μs_t{ℚ}               ) ::Pulse__Δ_BangBang{ℚ}
```

Note: No ℝ type-parameter (because there's no phase).

## Implementation

### The struct `Pulse__Δ_BangBang{ℚ}`

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
struct Pulse__Δ_BangBang{ℚ} <: Pulse                                                                #(3.2) struct Pulse__Δ_BangBang
    𝑒𝑣   ::NTuple{5, μs_t{ℚ} }          # events
    𝑟ꜛ   ::Radperμs_per_μs_t{ℚ}         # up-ramp rate
    𝛥    ::Rad_per_μs_t{ℚ}              # top plateau value
    𝑟ꜜ   ::Radperμs_per_μs_t{ℚ}         # down-ramp rate
end

function Pulse__Δ_BangBang{ℚ}(𝑡ᵒⁿ      ::μs_t{ℚ},                                                   #(3.2) constructor Pulse__Δ_BangBang
                              𝑡ᵒᶠᶠ     ::μs_t{ℚ},
                              𝑇        ::μs_t{ℚ},
                              𝛥_𝑡𝑎𝑟𝑔𝑒𝑡 ::Rad_per_μs_t{ℚ}
                              ;
                              𝛥ₘₐₓ           ::Rad_per_μs_t{ℚ},
                              𝛥ᵣₑₛ           ::Rad_per_μs_t{ℚ},
                              𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤   ::Radperμs_per_μs_t{ℚ},
                              𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 ::Radperμs_per_μs_t{ℚ}
                                               = 𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤,

                              𝑡ₘₐₓ           ::μs_t{ℚ},
                              𝑡ᵣₑₛ           ::μs_t{ℚ},
                              𝛥𝑡ₘᵢₙ          ::μs_t{ℚ}                ) ::
                                                          Pulse__Δ_BangBang{ℚ}   where{ℚ}



    0μs ≤ 𝑡ᵒⁿ                || throw(ArgumentError("Need  0μs ≤ 𝑡ᵒⁿ           ."))
          𝑡ᵒⁿ < 𝑡ᵒᶠᶠ         || throw(ArgumentError("Need        𝑡ᵒⁿ < 𝑡ᵒᶠᶠ    ."))
                𝑡ᵒᶠᶠ ≤ 𝑇     || throw(ArgumentError("Need              𝑡ᵒᶠᶠ ≤ 𝑇."))

    𝛥ₘₐₓ > 0/μs              || throw(ArgumentError("𝛥ₘₐₓ must be positive."))
    𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤 > 0/μs^2    || throw(ArgumentError("Max slew rate 𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤 must \
                                                     be positive."))
    𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 > 0/μs^2  || throw(ArgumentError("Max slew rate 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 must \
                                                     be positive."))
    𝛥_𝑡𝑎𝑟𝑔𝑒𝑡 % 𝛥ᵣₑₛ == 0/μs  || throw(ArgumentError("𝛥_𝑡𝑎𝑟𝑔𝑒𝑡 ($(𝛥_𝑡𝑎𝑟𝑔𝑒𝑡)) is not integer \
                                                     multiple of 𝛥ᵣₑₛ ($(𝛥ᵣₑₛ))."))
    -𝛥ₘₐₓ ≤ 𝛥_𝑡𝑎𝑟𝑔𝑒𝑡 ≤ +𝛥ₘₐₓ || throw(ArgumentError("𝛥_𝑡𝑎𝑟𝑔𝑒𝑡 ($(𝛥_𝑡𝑎𝑟𝑔𝑒𝑡)) is not in \
                                                     range [-𝛥ₘₐₓ,+𝛥ₘₐₓ] ($(𝛥ₘₐₓ))."))
    𝑡ᵒⁿ + 𝛥𝑡ₘᵢₙ ≤ 𝑡ᵒᶠᶠ       || throw(ArgumentError("Gap 𝑡ᵒⁿ → 𝑡ᵒᶠᶠ ($(𝑡ᵒᶠᶠ-𝑡ᵒⁿ)) \
                                                     smaller than 𝛥𝑡ₘᵢₙ ($(𝛥𝑡ₘᵢₙ))."))
    𝑡ᵒⁿ ≤ 0μs || 𝑡ᵒⁿ ≥ 𝛥𝑡ₘᵢₙ || throw(ArgumentError("Gap 0μs → 𝑡ᵒⁿ ($(𝑡ᵒⁿ)) \
                                                     smaller than 𝛥𝑡ₘᵢₙ ($(𝛥𝑡ₘᵢₙ))."))
    𝑡ᵒⁿ %  𝑡ᵣₑₛ == 0μs       || throw(ArgumentError("𝑡ᵒⁿ ($(𝑡ᵒⁿ)) is not integer multiple \
                                                     of 𝑡ᵣₑₛ ($(𝑡ᵣₑₛ)."))
    𝑡ᵒᶠᶠ % 𝑡ᵣₑₛ == 0μs       || throw(ArgumentError("𝑡ᵒᶠᶠ ($(𝑡ᵒᶠᶠ)) is not integer multiple \
                                                     of 𝑡ᵣₑₛ ($(𝑡ᵣₑₛ))."))



    # Warning! 𝛥_𝑡𝑎𝑟𝑔𝑒𝑡 can be negative!                                       𝗪𝗮𝗿𝗻𝗶𝗻𝗴!
    #          Signs of 𝑟ꜛ, 𝑟ꜜ, 𝛥_𝑡𝑎𝑟𝑔𝑒𝑡 must all be the same                         !
    #          for the stuff to work.                                                 !




    𝑟ꜛ        = ( 𝛥_𝑡𝑎𝑟𝑔𝑒𝑡 ≥ 0/μs ? 𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤   : -𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤   )
    𝑟ꜜ        = ( 𝛥_𝑡𝑎𝑟𝑔𝑒𝑡 ≥ 0/μs ? 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 : -𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 )
    𝑡ᵒⁿ⁻ᵗᵃʳ   = 𝛥_𝑡𝑎𝑟𝑔𝑒𝑡/ 𝑟ꜛ               # time from "on" to reaching target value
    𝑡ᵖᵉᵃᵏ     = min(𝑡ᵒⁿ⁻ᵗᵃʳ , 𝑡ᵒᶠᶠ-𝑡ᵒⁿ)    # time from "on" to peak value
    𝛥ᵖᵉᵃᵏ     = 𝑡ᵖᵉᵃᵏ⋅𝑟ꜛ                   # peak value
    𝑡ᵖᵉᵃᵏ⁻⁰   = 𝛥ᵖᵉᵃᵏ / 𝑟ꜜ                 # time from peak value to zero


    𝑒𝑣 ::NTuple{5, μs_t{ℚ} } =
        (
            # wait before pulse
            𝑡ᵒⁿ,
            # ramp up
            𝑡ᵒⁿ + 𝑡ᵖᵉᵃᵏ,
            # plateau
            𝑡ᵒᶠᶠ,
            #   ramp down
            𝑡ᵒᶠᶠ + 𝑡ᵖᵉᵃᵏ⁻⁰,
            # wait after pulse
            𝑇
        )

    𝑒𝑣[4] ≤ 𝑒𝑣[5]    || throw(ArgumentError("𝛥_BangBang pulse shape doesn't fit: gap \
                                             between 𝑡ᵒᶠᶠ=$(𝑡ᵒᶠᶠ) and 𝑇=$(𝑇) too small \
                                             for 𝛥_𝑡𝑎𝑟𝑔𝑒𝑡=$(𝛥_𝑡𝑎𝑟𝑔𝑒𝑡) with 𝛥ᵖᵉᵃᵏ=$(𝛥ᵖᵉᵃᵏ), \
                                             and 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 ($(𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤))."))

    return Pulse__Δ_BangBang{ℚ}(𝑒𝑣, 𝑟ꜛ, 𝛥ᵖᵉᵃᵏ, 𝑟ꜜ)
end

function _check(Δ::Pulse__Δ_BangBang{ℚ}) where{ℚ}                                                   #(3.2) _check() Pulse__Δ_BangBang
    @assert 0μs ≤ Δ.𝑒𝑣[1]        "Pulse__Δ_BangBang: \
                                  𝑒𝑣=$(Δ.𝑒𝑣) has negative time. This is a bug."
    @assert issorted(Δ.𝑒𝑣)       "Pulse__Δ_BangBang: \
                                  𝑒𝑣=$(Δ.𝑒𝑣) not sorted. This is a bug."
    @assert (
        Δ.𝛥≥0/μs && Δ.𝑟ꜛ>0/μs^2 && Δ.𝑟ꜜ>0/μs^2
        ||
        Δ.𝛥<0/μs && Δ.𝑟ꜛ<0/μs^2 && Δ.𝑟ꜜ<0/μs^2
        )                         "Pulse__Δ_BangBang: \
                                   sign mismatch between 𝛥,𝑟ꜛ,𝑟ꜜ. This is a bug."
    return true
end

function phase(Δ::Pulse__Δ_BangBang{ℚ})      where{ℚ}                                               #(3.2) phase() Pulse__Δ_BangBang
    throw(ErrorException("DAU catch: Δ-pulses have no phase."))
end

#
# This function is to demonstrate the pulse shape data, and maybe for plotting or whatnot.
#
function (Δ::Pulse__Δ_BangBang{ℚ})(𝑡 ::μs_t{𝕂}) ::Rad_per_μs_t{𝕂}   where{ℚ,𝕂}                      #(3.2) callable Pulse__Δ_BangBang

    (; 𝑒𝑣, 𝑟ꜛ, 𝑟ꜜ, 𝛥) = Δ

    β = (2^30+1)//2^30
    if            𝑡 < 0μs            throw(DomainError(𝑡,"Time cannot be negative."))
    elseif  0μs   ≤ 𝑡 ≤ 𝑒𝑣[1]        return 𝕂(0)/μs
    elseif  𝑒𝑣[1] < 𝑡 < 𝑒𝑣[2]        return ( 𝑡-𝑒𝑣[1] )⋅𝑟ꜛ
    elseif  𝑒𝑣[2] ≤ 𝑡 ≤ 𝑒𝑣[3]        return 𝛥
    elseif  𝑒𝑣[3] < 𝑡 < 𝑒𝑣[4]        return ( 𝑒𝑣[4]-𝑡 )⋅𝑟ꜜ
    elseif  𝑒𝑣[4] ≤ 𝑡 ≤ 𝑒𝑣[5]⋅β      return 𝕂(0)/μs
    elseif  𝑒𝑣[5]⋅β < 𝑡              throw(DomainError(𝑡,"Time exceeds upper bound, \
                                                          𝑇=$(𝑒𝑣[5])."))
    else                             @assert false "It's the Unitful-comparison's bug!"
    end
end #^ callable Pulse__Δ_BangBang

function 𝑎𝑣𝑔(Δ ::Pulse__Δ_BangBang{ℚ},                                                              #(3.2) 𝑎𝑣𝑔() Pulse__Δ_BangBang
             𝑡 ::μs_t{𝕂}
             ;
             𝛥𝑡 ::μs_t{𝕂}               ) ::Rad_per_μs_t{𝕂}       where{ℚ,𝕂}

    (;𝑒𝑣) = Δ
    𝑡ᵉⁿᵈ  = 𝑡+𝛥𝑡
    sum   = 𝕂(0)
    for j = 1 : length(𝑒𝑣)-1
        if 𝑡 < 𝑒𝑣[j+1] && 𝑒𝑣[j] < 𝑡ᵉⁿᵈ
            𝑠ⱼ = max(𝑒𝑣[j], 𝑡)
            𝑡ⱼ = min(𝑡ᵉⁿᵈ, 𝑒𝑣[j+1])
            if 𝑠ⱼ < 𝑡ⱼ
                𝛿ₛ = Δ(𝑠ⱼ)
                𝛿ₜ = Δ(𝑡ⱼ)
                sum += (𝑡ⱼ-𝑠ⱼ)⋅( 𝛿ₛ + 𝛿ₜ )/2
            end
        end
    end
    return sum/𝛥𝑡
end #^ 𝑎𝑣𝑔()

function 𝑠𝑡𝑒𝑝(Δ::Pulse__Δ_BangBang{ℚ},                                                              #(3.2) 𝑠𝑡𝑒𝑝() Pulse__Δ_BangBang
              𝑡 ::μs_t{𝕂}
              ;
              ε ::𝕂                     ) ::μs_t{𝕂}   where{ℚ,𝕂}

    @assert ε > 0   "ε ≤ 0"

    (; 𝑒𝑣, 𝑟ꜛ, 𝑟ꜜ) = Δ

    𝑠𝑡𝑝(𝑟) = √( 4ε / abs(𝑟) )
#    𝑠𝑡𝑝(𝑟) =  4ε / abs(𝑟) / μs

    # β = (2^30+1)//2^30
    if            𝑡 < 0μs      throw(DomainError(𝑡,"Time 𝑡=$(BigFloat(𝑡)) cannot be \
                                                      negative."))
    elseif  0μs < 𝑒𝑣[1] - 𝑡    return                𝑒𝑣[1] - 𝑡
    elseif  0μs < 𝑒𝑣[2] - 𝑡    return min( 𝑠𝑡𝑝(𝑟ꜛ) ,  𝑒𝑣[2] - 𝑡 )
    elseif  0μs < 𝑒𝑣[3] - 𝑡    return                𝑒𝑣[3] - 𝑡
    elseif  0μs < 𝑒𝑣[4] - 𝑡    return min( 𝑠𝑡𝑝(𝑟ꜜ) ,  𝑒𝑣[4] - 𝑡 )
    elseif  0μs < 𝑒𝑣[5] - 𝑡    return                𝑒𝑣[5] - 𝑡
    elseif  𝑒𝑣[5] - 𝑡 < 0μs    throw(DomainError(𝑡,"Time 𝑡=$(BigFloat(𝑡)) exceeds upper \
                                                    bound 𝑇=$(BigFloat(𝑒𝑣[5]))."))
    else                         @assert false "It's the Unitful-comparison's bug!"
    end

end #^ 𝑠𝑡𝑒𝑝()


function plotpulse(Δ::Pulse__Δ_BangBang) ::NamedTuple                                               #(3.2) plotpulse() Pulse__Δ_BangBang
    (; 𝑒𝑣, 𝛥) = Δ
    x⃗ = [ (0//1)μs  ,  𝑒𝑣[1]     , 𝑒𝑣[2] , 𝑒𝑣[3] , 𝑒𝑣[4]     ,  𝑒𝑣[5]     ]
    y⃗ = [ (0//1)/μs ,  (0//1)/μs , 𝛥     , 𝛥     , (0//1)/μs ,  (0//1)/μs ]
    return (x⃗=x⃗, y⃗=y⃗)
    # 𝙇𝙖𝙯𝙮 𝙫𝙚𝙧𝙨𝙞𝙤𝙣 (𝙙𝙚𝙛𝙚𝙧𝙧𝙞𝙣𝙜 𝙩𝙤 𝙘𝙖𝙡𝙡𝙖𝙗𝙡𝙚):
    #    𝑋 = Iterators.flatten( [ [(0//1)μs], (𝑡 for 𝑡 ∈ Δ.𝑒𝑣) ] )
    #    return (  x⃗ = collect(𝑋),
    #              y⃗ = [ Δ(𝑥) for 𝑥 ∈ 𝑋 ]  )
end



# ——————————————————————————————————————————————————————————————————————————————————————————————————— 3.3. Ω_BangBang Pulse

@doc raw"""
Struct `Pulse__Ω_BangBang` `<:` `Pulse`

!!! note "Note!"
    The Rabi frequency technically cannot be negative, but we allow negative values in the
    constructors of the Ω-pulse shapes.  The sign is then hiden in the phase.  (Only reason for
    the phase, indeed.)  The additional type-parameter ℝ is needed to encode the phase.

## Constructor
```julia
Pulse__Ω_BangBang{ℚ,ℝ}( 𝑡ᵒⁿ      ::μs_t{ℚ},
                        𝑡ᵒᶠᶠ     ::μs_t{ℚ},
                        𝑇        ::μs_t{ℚ},
                        𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 ::Rad_per_μs_t{ℚ}
                        ;
                        𝛺ₘₐₓ           ::Rad_per_μs_t{ℚ},
                        𝛺ᵣₑₛ           ::Rad_per_μs_t{ℚ},
                        𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤   ::Radperμs_per_μs_t{ℚ},
                        𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 ::Radperμs_per_μs_t{ℚ} = 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤,
                        φᵣₑₛ           ::ℚ,
                        𝑡ₘₐₓ           ::μs_t{ℚ},
                        𝑡ᵣₑₛ           ::μs_t{ℚ},
                        𝛥𝑡ₘᵢₙ          ::μs_t{ℚ}               ) ::Pulse__Ω_BangBang{ℚ,ℝ}
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
struct Pulse__Ω_BangBang{ℚ,ℝ} <: Pulse                                                              #(3.3) struct Pulse__Ω_BangBang
    γ    ::Complex{ℝ}                   # phase
    𝑒𝑣   ::NTuple{5, μs_t{ℚ} }          # events
    𝑟ꜛ   ::Radperμs_per_μs_t{ℚ}         # up-ramp rate
    𝛺    ::Rad_per_μs_t{ℚ}              # top plateau value
    𝑟ꜜ   ::Radperμs_per_μs_t{ℚ}         # down-ramp rate
end

function Pulse__Ω_BangBang{ℚ,ℝ}(𝑡ᵒⁿ      ::μs_t{ℚ},                                                 #(3.3) constructor Pulse__Ω_BangBang
                                𝑡ᵒᶠᶠ     ::μs_t{ℚ},
                                𝑇        ::μs_t{ℚ},
                                𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 ::Rad_per_μs_t{ℚ}
                                ;
                                𝛺ₘₐₓ           ::Rad_per_μs_t{ℚ},
                                𝛺ᵣₑₛ           ::Rad_per_μs_t{ℚ},
                                𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤   ::Radperμs_per_μs_t{ℚ},
                                𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 ::Radperμs_per_μs_t{ℚ}
                                                 = 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤,
                                φᵣₑₛ           ::ℚ,                     # "\varphi"
                                𝑡ₘₐₓ           ::μs_t{ℚ},
                                𝑡ᵣₑₛ           ::μs_t{ℚ},
                                𝛥𝑡ₘᵢₙ          ::μs_t{ℚ}                ) ::
                                                          Pulse__Ω_BangBang{ℚ,ℝ}   where{ℚ,ℝ}

    ℂ = Complex{ℝ}

    0μs ≤ 𝑡ᵒⁿ                || throw(ArgumentError("Need  0μs ≤ 𝑡ᵒⁿ           ."))
          𝑡ᵒⁿ < 𝑡ᵒᶠᶠ         || throw(ArgumentError("Need        𝑡ᵒⁿ < 𝑡ᵒᶠᶠ    ."))
                𝑡ᵒᶠᶠ ≤ 𝑇     || throw(ArgumentError("Need              𝑡ᵒᶠᶠ ≤ 𝑇."))

    𝛺ₘₐₓ > 0/μs              || throw(ArgumentError("𝛺ₘₐₓ must be positive."))
    𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤 > 0/μs^2    || throw(ArgumentError("Max slew rate 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤 must \
                                                     be positive."))
    𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 > 0/μs^2  || throw(ArgumentError("Max slew rate 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 must \
                                                     be positive."))
    𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 % 𝛺ᵣₑₛ == 0/μs  || throw(ArgumentError("𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 ($(𝛺_𝑡𝑎𝑟𝑔𝑒𝑡)) is not integer \
                                                     multiple of 𝛺ᵣₑₛ ($(𝛺ᵣₑₛ))."))
    -𝛺ₘₐₓ ≤ 𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 ≤ +𝛺ₘₐₓ || throw(ArgumentError("𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 ($(𝛺_𝑡𝑎𝑟𝑔𝑒𝑡)) is not in \
                                                     range [-𝛺ₘₐₓ,+𝛺ₘₐₓ] ($(𝛺ₘₐₓ))."))
    𝑡ᵒⁿ + 𝛥𝑡ₘᵢₙ ≤ 𝑡ᵒᶠᶠ       || throw(ArgumentError("Gap 𝑡ᵒⁿ → 𝑡ᵒᶠᶠ ($(𝑡ᵒᶠᶠ-𝑡ᵒⁿ)) \
                                                     smaller than 𝛥𝑡ₘᵢₙ ($(𝛥𝑡ₘᵢₙ))."))
    𝑡ᵒⁿ ≤ 0μs || 𝑡ᵒⁿ ≥ 𝛥𝑡ₘᵢₙ || throw(ArgumentError("Gap 0μs → 𝑡ᵒⁿ ($(𝑡ᵒⁿ)) \
                                                     smaller than 𝛥𝑡ₘᵢₙ ($(𝛥𝑡ₘᵢₙ))."))
    𝑡ᵒⁿ %  𝑡ᵣₑₛ == 0μs       || throw(ArgumentError("𝑡ᵒⁿ ($(𝑡ᵒⁿ)) is not integer multiple \
                                                     of 𝑡ᵣₑₛ ($(𝑡ᵣₑₛ)."))
    𝑡ᵒᶠᶠ % 𝑡ᵣₑₛ == 0μs       || throw(ArgumentError("𝑡ᵒᶠᶠ ($(𝑡ᵒᶠᶠ)) is not integer multiple \
                                                     of 𝑡ᵣₑₛ ($(𝑡ᵣₑₛ))."))

    γ::ℂ =
        if 𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 < 0/μs
            𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 = -𝛺_𝑡𝑎𝑟𝑔𝑒𝑡           # Warning! Change sign of 𝛺_𝑡𝑎𝑟𝑔𝑒𝑡  𝗪𝗮𝗿𝗻𝗶𝗻𝗴!
            cis( δround(ℝ(π);δ=φᵣₑₛ) )     #          𝛺_𝑡𝑎𝑟𝑔𝑒𝑡 (and 𝑟ꜛ, 𝑟ꜜ)           !
        else                               #          must be positive, sign          !
            ℂ(1)                           #          is hidden in the phase.         !
        end


    𝑟ꜛ        = 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤
    𝑟ꜜ        = 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤
    𝑡ᵒⁿ⁻ᵗᵃʳ   = 𝛺_𝑡𝑎𝑟𝑔𝑒𝑡/ 𝑟ꜛ               # time from "on" to reaching target value
    𝑡ᵖᵉᵃᵏ     = min(𝑡ᵒⁿ⁻ᵗᵃʳ , 𝑡ᵒᶠᶠ-𝑡ᵒⁿ)    # time from "on" to peak value
    𝛺ᵖᵉᵃᵏ     = 𝑡ᵖᵉᵃᵏ⋅𝑟ꜛ                   # peak value
    𝑡ᵖᵉᵃᵏ⁻⁰   = 𝛺ᵖᵉᵃᵏ / 𝑟ꜜ                 # time from peak value to zero


    𝑒𝑣 ::NTuple{5, μs_t{ℚ} } =
        (
            # wait before pulse
            𝑡ᵒⁿ,
            # ramp up
            𝑡ᵒⁿ + 𝑡ᵖᵉᵃᵏ,
            # plateau
            𝑡ᵒᶠᶠ,
            #   ramp down
            𝑡ᵒᶠᶠ + 𝑡ᵖᵉᵃᵏ⁻⁰,
            # wait after pulse
            𝑇
        )

    𝑒𝑣[4] ≤ 𝑒𝑣[5]    || throw(ArgumentError("𝛺_BangBang pulse shape doesn't fit: gap \
                                             between 𝑡ᵒᶠᶠ=$(𝑡ᵒᶠᶠ) and 𝑇=$(𝑇) too small \
                                             for 𝛺_𝑡𝑎𝑟𝑔𝑒𝑡=$(𝛺_𝑡𝑎𝑟𝑔𝑒𝑡) with 𝛺ᵖᵉᵃᵏ=$(𝛺ᵖᵉᵃᵏ), \
                                             and 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 ($(𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤))."))

    return Pulse__Ω_BangBang{ℚ,ℝ}(γ, 𝑒𝑣, 𝑟ꜛ, 𝛺ᵖᵉᵃᵏ, 𝑟ꜜ)
end

function _check(Ω::Pulse__Ω_BangBang{ℚ,ℝ}) where{ℚ,ℝ}                                               #(3.3) _check() Pulse__Ω_BangBang
    @assert 0μs ≤ Ω.𝑒𝑣[1]          "Pulse__Ω_BangBang: \
                                    𝑒𝑣=$(Ω.𝑒𝑣) has negative time. This is a bug."
    @assert issorted(Ω.𝑒𝑣)         "Pulse__Ω_BangBang: \
                                    𝑒𝑣=$(Ω.𝑒𝑣) not sorted. This is a bug."
    @assert Ω.𝛺 ≥ 0/μs             "Pulse__Ω_BangBang: \
                                    negative 𝛺=$(Ω.𝛺). This is a bug."
    @assert Ω.𝑟ꜛ > 0/μs^2  &&
        Ω.𝑟ꜜ > 0/μs^2               "Pulse__Ω_BangBang: \
                                    negative slew rate (𝑟ꜛ=$(Ω.𝑟ꜛ), 𝑟ꜜ=$(Ω.𝑟ꜜ)). \
                                    This is a bug."
    return true
end

function phase(Ω::Pulse__Ω_BangBang{ℚ,ℝ}) ::Complex{ℝ}      where{ℚ,ℝ}                              #(3.3) phase() Pulse__Ω_BangBang
    # let's take the opportunity to run some checks:
    _check(Ω)

    return Ω.γ
end

#
# This function is to demonstrate the pulse shape data, and maybe for plotting or whatnot.
#
function (Ω::Pulse__Ω_BangBang{ℚ,ℝ})(𝑡 ::μs_t{𝕂}) ::Rad_per_μs_t{𝕂}   where{ℚ,ℝ,𝕂}                  #(3.3) callable Pulse__Ω_BangBang

    (; 𝑒𝑣, 𝑟ꜛ, 𝑟ꜜ, 𝛺) = Ω

    β = (2^30+1)//2^30
    if            𝑡 < 0μs            throw(DomainError(𝑡,"Time cannot be negative."))
    elseif  0μs   ≤ 𝑡 ≤ 𝑒𝑣[1]        return 𝕂(0)/μs
    elseif  𝑒𝑣[1] < 𝑡 < 𝑒𝑣[2]        return ( 𝑡-𝑒𝑣[1] )⋅𝑟ꜛ
    elseif  𝑒𝑣[2] ≤ 𝑡 ≤ 𝑒𝑣[3]        return 𝛺
    elseif  𝑒𝑣[3] < 𝑡 < 𝑒𝑣[4]        return ( 𝑒𝑣[4]-𝑡 )⋅𝑟ꜜ
    elseif  𝑒𝑣[4] ≤ 𝑡 ≤ 𝑒𝑣[5]⋅β      return 𝕂(0)/μs
    elseif  𝑒𝑣[5]⋅β < 𝑡              throw(DomainError(𝑡,"Time exceeds upper bound, \
                                                          𝑇=$(𝑒𝑣[5])."))
    else                             @assert false "It's the Unitful-comparison's bug!"
    end
end #^ callable Pulse__Ω_BangBang

function 𝑎𝑣𝑔(Ω ::Pulse__Ω_BangBang{ℚ,ℝ},                                                            #(3.3) 𝑎𝑣𝑔() Pulse__Ω_BangBang
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
                sum += (𝑡ⱼ-𝑠ⱼ)⋅( 𝜔ₛ + 𝜔ₜ )/2
            end
        end
    end
    return sum/𝛥𝑡
end #^ 𝑎𝑣𝑔()

function 𝑠𝑡𝑒𝑝(Ω::Pulse__Ω_BangBang{ℚ,ℝ},                                                            #(3.3) 𝑠𝑡𝑒𝑝() Pulse__Ω_BangBang
              𝑡 ::μs_t{𝕂}
              ;
              ε ::𝕂                     ) ::μs_t{𝕂}   where{ℚ,ℝ,𝕂}

    @assert ε > 0   "ε ≤ 0"

    (; 𝑒𝑣, 𝑟ꜛ, 𝑟ꜜ) = Ω

    𝑠𝑡𝑝(𝑟) = √( 4ε / abs(𝑟) )
#    𝑠𝑡𝑝(𝑟) =  4ε / abs(𝑟) / μs

    # β = (2^30+1)//2^30
    if            𝑡 < 0μs        throw(DomainError(𝑡,"Time 𝑡=$(BigFloat(𝑡)) cannot be \
                                                      negative."))
    elseif  0μs < 𝑒𝑣[1] - 𝑡    return                𝑒𝑣[1]-𝑡
    elseif  0μs < 𝑒𝑣[2] - 𝑡    return min( 𝑠𝑡𝑝(𝑟ꜛ) ,  𝑒𝑣[2]-𝑡 )
    elseif  0μs < 𝑒𝑣[3] - 𝑡    return                𝑒𝑣[3]-𝑡
    elseif  0μs < 𝑒𝑣[4] - 𝑡    return min( 𝑠𝑡𝑝(𝑟ꜜ) ,  𝑒𝑣[4]-𝑡 )
    elseif  0μs < 𝑒𝑣[5] - 𝑡    return                𝑒𝑣[5]-𝑡
    elseif  𝑒𝑣[5] - 𝑡 < 0μs    throw(DomainError(𝑡,"Time 𝑡=$(BigFloat(𝑡)) exceeds upper \
                                                      bound 𝑇=$(BigFloat(𝑒𝑣[5]))."))
    else                         @assert false "It's the Unitful-comparison's bug!"
    end

end #^ 𝑠𝑡𝑒𝑝()


function plotpulse(Ω::Pulse__Ω_BangBang) ::NamedTuple                                               #(3.3) plotpulse() Pulse__Ω_BangBang
    (; 𝑒𝑣, 𝛺) = Ω
    x⃗ = [ (0//1)μs  ,  𝑒𝑣[1]     , 𝑒𝑣[2] , 𝑒𝑣[3] , 𝑒𝑣[4]     ,  𝑒𝑣[5]     ]
    y⃗ = [ (0//1)/μs ,  (0//1)/μs , 𝛺     , 𝛺     , (0//1)/μs ,  (0//1)/μs ]
    return (x⃗=x⃗, y⃗=y⃗)
    # 𝙇𝙖𝙯𝙮 𝙫𝙚𝙧𝙨𝙞𝙤𝙣 (𝙙𝙚𝙛𝙚𝙧𝙧𝙞𝙣𝙜 𝙩𝙤 𝙘𝙖𝙡𝙡𝙖𝙗𝙡𝙚):
    #    𝑋 = Iterators.flatten( [ [(0//1)μs], (𝑡 for 𝑡 ∈ Ω.𝑒𝑣) ] )
    #    return (  x⃗ = collect(𝑋),
    #              y⃗ = [ Ω(𝑥) for 𝑥 ∈ 𝑋 ]  )
end

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 4. Sub-module `Schrödinger`

include("Schrödinger.mod.jl")
import .Schrödinger: schröd!
import .Schrödinger: 𝔑, 𝔛

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 5. Sub-module `HW_Descriptions`

include("HW_Descriptions.mod.jl")

end # module DOT_RydSim
# EOF
