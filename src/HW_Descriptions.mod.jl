########################################################################
#                                                                      #
# DOT_RydSim/src/HW_Descriptions.mod.jl                                #
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
#    0.  Header, Imports & Helpers                                                                            |
#                                                                                                             |
#        1.1.  Exports                                                                                        |
#        1.2.  Imports                                                                                        |
#                                                                                                             |
#                                                                                                             |
#    1.  The Types                                                                                            |
#                                                                                                             |
#    2.  File Input                                                                                           |
#                                                                                                             |
#        2.0.  Selection                                                                                      |
#        2.1.  read AWS-QuEra                                                                                 |
#                                                                                                             |
#                                                                                                             |
#    3.  Default                                                                                              |
#                                                                                                             |
#    4.  Cheatify                                                                                             |
#                                                                                                             |
#—————————————————————————————————————————————————————————————————————————————————————————————————————————————+


# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 0. Header, Imports & Helpers
"""
Module `HW_Descriptions`

Exports:
  * Type `HW_Descr{ℚ}`
  * Function
    `default_HW_Descr(options::Symbol...; ℤ=Int128) :: HW_Descr{Rational{ℤ}}`
  * Function
    `fileread_HW_Descr( `*FileType*` ; filename ::String, ℤ=Int128) ::HW_Descr{Rational{ℤ}}`
  * Functions
    ```julia
        cheatify_res( ::HW_Descr{ℚ} ; factor::ℚ =1//1000) ::HW_Descr{ℚ} where{ℚ}
        cheatify_𝑡ₘₐₓ(::HW_Descr{ℚ} ; factor::ℚ =  10//1) ::HW_Descr{ℚ} where{ℚ}
    ```
    which returns a copied & modified version of the HW-description with: Ω,Δ resolutions
    improved; 𝑡ₘₐₓ improved, resp.

The argument "*FileType*" is a type.  The type constants are not exported, and they are:

  * `HW_AWS_QuEra`
  * and that's it for now.

Both `default_HW_Descr()` and `fileread_HW_Descr(HW_AWS_QuEra,...)` accept the following
keyword arguments:

  * `Ω_downslew_factor`
  * `Δ_downslew_factor`
both of type rational.

Admissible `options` for `default_HW_Descr()` are, currently:
  * none — read from file "hw_default.json", which describes QuERA's device, and has been
    downloaded from AWS .
"""
module HW_Descriptions
export HW_Descr, default_HW_Descr, fileread_HW_Descr,  cheatify_res, cheatify_𝑡ₘₐₓ

import ..μs_t, ..Rad_per_μs_t, ..Radperμs_per_μs_t
import ..is_δrounded

using  ..DOT_NiceMath


using Unitful
using Unitful: Length, 𝐋, 𝐓, m, s
@derived_dimension Length⁶_per_Time 𝐋^6/𝐓

using JSON

import Base: rationalize
rationalize(ℤ,x::Integer) = ℤ(x)//ℤ(1)       # bulit-in `rationalize()` works only for float.

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1. The Types

@kwdef struct Lattice_Descr
	width          ::Length
	height         ::Length
	radialΔₘᵢₙ     ::Length
	verticalΔₘᵢₙ   ::Length
	posᵣₑₛ         ::Length
	numsitesₘₐₓ    ::Int32
end

"""
Struct `HW_Descr{ℚ}`

Holds the relevant data of the Rydberg atom array quantum device.
"""
@kwdef struct HW_Descr{ℚ}

    lattice        ::Lattice_Descr

    𝐶₆             ::Length⁶_per_Time

    𝛺ₘₐₓ           ::Rad_per_μs_t{ℚ}
    𝛺ᵣₑₛ           ::Rad_per_μs_t{ℚ}
    𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤   ::Radperμs_per_μs_t{ℚ}
    𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 ::Radperμs_per_μs_t{ℚ}

    𝛥ₘₐₓ           ::Rad_per_μs_t{ℚ}
    𝛥ᵣₑₛ           ::Rad_per_μs_t{ℚ}
    𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤   ::Radperμs_per_μs_t{ℚ}
    𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 ::Radperμs_per_μs_t{ℚ}

    φₘₐₓ           ::ℚ
    φᵣₑₛ           ::ℚ                     # "\varphi"

    𝑡ₘₐₓ           ::μs_t{ℚ}
    𝑡ᵣₑₛ           ::μs_t{ℚ}
    𝛥𝑡ₘᵢₙ          ::μs_t{ℚ}
end



# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2. File input

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.0. Selection

abstract type                 HW_Descr_Format end
struct        HW_AWS_QuEra <: HW_Descr_Format end

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.1. read AWS-QuEra

function fileread_HW_Descr(::Type{HW_AWS_QuEra}
                           ;
                           filename          ::String,
                           ℤ                 ::Type{<:Integer}   = Int128,
                           Ω_downslew_factor ::Rational          = 3//1,
                           Δ_downslew_factor ::Rational          = 1//3)

    ℚ       = Rational{ℤ}
    rat(x)  = rationalize(ℤ,x)

    j       = JSON.parse( read(filename,String) )
    la_area = j["lattice"]["area"]
    la_geo  = j["lattice"]["geometry"]
    R       = j["rydberg"]
    Rg      = R["rydbergGlobal"]

    lattice = Lattice_Descr(
	width         = u"μm"(  rat( la_area["width"]  )m   ),
	height        = u"μm"(  rat( la_area["height"] )m   ),
	radialΔₘᵢₙ    = u"μm"(  rat( la_geo["spacingRadialMin"]     )m   ),
	verticalΔₘᵢₙ  = u"μm"(  rat( la_geo["spacingVerticalMin"]   )m   ),
	posᵣₑₛ        = u"μm"(  rat( la_geo["positionResolution"]   )m   ),
	numsitesₘₐₓ   = la_geo["numberSitesMax"]
    )

    return 	HW_Descr{ℚ}(
        ;
	lattice,
	𝐶₆             = u"μm^6/μs"( R["c6Coefficient"]m^6/s ),

	𝛺ₘₐₓ           = u"μs^(-1)"(       rat( Rg["rabiFrequencyRange"][2]    )/s   )::Rad_per_μs_t{ℚ},
	𝛺ᵣₑₛ           = u"μs^(-1)"(       rat( Rg["rabiFrequencyResolution"]  )/s   )::Rad_per_μs_t{ℚ},
        𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤   = u"μs^(-2)"( rat( Rg["rabiFrequencySlewRateMax"] )/s^2 )      ::Radperμs_per_μs_t{ℚ},
        𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 = (
            Ω_downslew_factor ⋅ u"μs^(-2)"( rat( Rg["rabiFrequencySlewRateMax"] )/s^2 )
        )::Radperμs_per_μs_t{ℚ},

	𝛥ₘₐₓ           = u"μs^(-1)"(       rat( Rg["detuningRange"][2]    )/s   )     ::Rad_per_μs_t{ℚ},
	𝛥ᵣₑₛ           = u"μs^(-1)"(       rat( Rg["detuningResolution"]  )/s   )     ::Rad_per_μs_t{ℚ},
        𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤   = u"μs^(-2)"( rat( Rg["detuningSlewRateMax"] )/s^2 )           ::Radperμs_per_μs_t{ℚ},
        𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 = (
            Δ_downslew_factor ⋅ u"μs^(-2)"( rat( Rg["detuningSlewRateMax"] )/s^2 )
        )::Radperμs_per_μs_t{ℚ},

	φₘₐₓ           = rat( Rg["phaseRange"][2]   )                                 ::ℚ,
	φᵣₑₛ           = rat( Rg["phaseResolution"] )                                 ::ℚ,

	𝑡ₘₐₓ           = u"μs"( rat( Rg["timeMax"]        )s )                        ::μs_t{ℚ},
	𝑡ᵣₑₛ           = u"μs"( rat( Rg["timeResolution"] )s )                        ::μs_t{ℚ},
	𝛥𝑡ₘᵢₙ          = u"μs"( rat( Rg["timeDeltaMin"]   )s )                        ::μs_t{ℚ}
    )
end #^ input_HW_Descr()

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 3. Default

function default_HW_Descr(O :: Symbol...    ;
                          ℤ                 ::Type{<:Integer} = Int128,
                          Ω_downslew_factor ::Rational        = 1//1,
                          Δ_downslew_factor ::Rational        = 1//1)

    ALL_O = Set{Symbol}([])
    O ⊆ ALL_O || throw(ArgumentError("Unrecognized options: $(setdiff(O,ALL_O))"))

    filename = pkgdir(@__MODULE__,
                      "Resources",
                      "hw_default.json"
                      )

    fileread_HW_Descr(HW_AWS_QuEra ;
                      filename,
                      ℤ,
                      Ω_downslew_factor, Δ_downslew_factor)
end #^ default_HW_Descr()

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 4. Cheatify

function cheatify_res(hw::HW_Descr{𝐐}  ; factor ::𝐐 =1//1000) ::HW_Descr{𝐐} where{𝐐}

    factor > 0                       || throw(ArgumentError("`factor` must be positive"))

    (;
     lattice,
     𝐶₆,
     𝛺ₘₐₓ,
     𝛺ᵣₑₛ,
     𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤,
     𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
     𝛥ₘₐₓ,
     𝛥ᵣₑₛ,
     𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤,
     𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
     φₘₐₓ,
     φᵣₑₛ,
     𝑡ₘₐₓ,
     𝑡ᵣₑₛ,
     𝛥𝑡ₘᵢₙ            ) = hw

    new_𝛺ᵣₑₛ = factor⋅𝛺ᵣₑₛ
    new_𝛥ᵣₑₛ = factor⋅𝛥ᵣₑₛ

    is_δrounded(𝛺ₘₐₓ;𝛿=new_𝛺ᵣₑₛ)  || throw(ArgumentError("new 𝛺ᵣₑₛ |̸ 𝛺ₘₐₓ"))
    is_δrounded(𝛥ₘₐₓ;𝛿=new_𝛥ᵣₑₛ)  || throw(ArgumentError("new 𝛥ᵣₑₛ |̸ 𝛥ₘₐₓ"))

    return HW_Descr{𝐐}(;
                       lattice,
                       𝐶₆,
                       𝛺ₘₐₓ,
                       𝛺ᵣₑₛ              = new_𝛺ᵣₑₛ,
                       𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤,
                       𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
                       𝛥ₘₐₓ,
                       𝛥ᵣₑₛ              = new_𝛥ᵣₑₛ,
                       𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤,
                       𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
                       φₘₐₓ,
                       φᵣₑₛ,
                       𝑡ₘₐₓ,
                       𝑡ᵣₑₛ,
                       𝛥𝑡ₘᵢₙ)
end #^ cheatify_res()

function cheatify_𝑡ₘₐₓ(hw::HW_Descr{𝐐}  ; factor ::𝐐 =10//1) ::HW_Descr{𝐐} where{𝐐}

    factor > 0                       || throw(ArgumentError("`factor` must be positive"))

    (;
     lattice,
     𝐶₆,
     𝛺ₘₐₓ,
     𝛺ᵣₑₛ,
     𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤,
     𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
     𝛥ₘₐₓ,
     𝛥ᵣₑₛ,
     𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤,
     𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
     φₘₐₓ,
     φᵣₑₛ,
     𝑡ₘₐₓ,
     𝑡ᵣₑₛ,
     𝛥𝑡ₘᵢₙ            ) = hw

    new_𝑡ₘₐₓ = factor⋅𝑡ₘₐₓ

    is_δrounded(new_𝑡ₘₐₓ  ;𝛿=𝑡ᵣₑₛ)  || throw(ArgumentError("𝑡ᵣₑₛ |̸ new 𝑡ₘₐₓ"))
    new_𝑡ₘₐₓ > 𝛥𝑡ₘᵢₙ                || throw(ArgumentError("new 𝑡ₘₐₓ must be > 𝛥𝑡ₘᵢₙ"))

    return HW_Descr{𝐐}(;
                       lattice,
                       𝐶₆,
                       𝛺ₘₐₓ,
                       𝛺ᵣₑₛ,
                       𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤,
                       𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
                       𝛥ₘₐₓ,
                       𝛥ᵣₑₛ,
                       𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤,
                       𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
                       φₘₐₓ,
                       φᵣₑₛ,
                       𝑡ₘₐₓ              = new_𝑡ₘₐₓ,
                       𝑡ᵣₑₛ,
                       𝛥𝑡ₘᵢₙ)
end #^ cheatify_𝑡ₘₐₓ()



end #^ module HW_Descriptions
