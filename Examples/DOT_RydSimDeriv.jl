module DOT_RydSimDeriv
export load_hw
export get_hw_data, get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺
export evf_Ω, evf_Ω

using DOT_NiceMath
using DOT_NiceMath.NumbersF64

using DOT_RydSim
using DOT_RydSim:
    μs_t,
    Rad_per_μs_t


using DOT_RydSim.HW_Descriptions:
    HW_Descr,
	default_HW_Descr,
	fileread_HW_Descr,
	HW_AWS_QuEra


using Unitful: μs
using LinearAlgebra: Hermitian


"""
Function
```julia
load_hw(filename = :default
		; Ω_downslew_factor = 1//1,
		  Δ_downslew_factor = 1//1  ) ::HW_Descr
```

The `filename` can be either a string identifying a file name, or the symbol `:default`.

The only file type currently supported is AWS-QuEra's (`HW_AWS_QuEra`).
"""
function load_hw(filename ::String
				;
				Ω_downslew_factor = 1//1,
				Δ_downslew_factor = 1//1              ) ::HW_Descr{ℚ}

	return fileread_HW_Descr(HW_AWS_QuEra
							;   filename,
								ℤ,
								Ω_downslew_factor,
								Δ_downslew_factor )
end

function load_hw(select ::Symbol =:default
				;
				Ω_downslew_factor = 1//1,
				Δ_downslew_factor = 1//1              ) ::HW_Descr{ℚ}

	@assert select == :default  "What?!??"
	return default_HW_Descr(;
							ℤ,
							Ω_downslew_factor,
							Δ_downslew_factor)
end





_NT = @NamedTuple{
				_blah::Nothing,
				𝛺ₘₐₓ        ::Rad_per_μs_t{ℚ},
				𝛺ᵣₑₛ        ::Rad_per_μs_t{ℚ},
				𝛥ₘₐₓ        ::Rad_per_μs_t{ℚ}, 𝛥ᵣₑₛ::Rad_per_μs_t{ℚ},
				𝑡ᵒᶠᶠₘₐₓ     ::μs_t{ℚ},
				𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ ::μs_t{ℚ},
				𝑡ᵣₑₛ        ::μs_t{ℚ},
				𝑡ₘₐₓ        ::μs_t{ℚ}
		}

@doc raw"""
Function `get_hw_data(::HW_Descr) ::NamedTuple`

Returns a named tuple with the following fields, all of
unitful rational number types:
* `𝛺ₘₐₓ`, `𝛺ᵣₑₛ`; 
* `𝛥ₘₐₓ` `𝛥ᵣₑₛ`;
* `𝑡ᵣₑₛ`;
* `𝑡ₘₐₓ`           — max total evolution time
* `𝑡ᵒᶠᶠₘₐₓ`        — largest switch-off time which allows full range of 𝛺 and 𝛥
* `𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ`    — smallest duration ``t^{\text{off}}-t^{\text{on}}``
  which allows full range of 𝛺 and 𝛥
"""
function get_hw_data(hw ::HW_Descr{ℚ}) ::_NT

	( ; 𝛺ₘₐₓ, 𝛺ᵣₑₛ, 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤, φᵣₑₛ,
		𝛥ₘₐₓ, 𝛥ᵣₑₛ, 𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
		𝑡ₘₐₓ, 𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ                               ) = hw

	𝛺_𝑢𝑝𝑡𝑖𝑚𝑒 = 𝛺ₘₐₓ / 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤 ; 𝛺_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒 = 𝛺ₘₐₓ / 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤
	𝛥_𝑢𝑝𝑡𝑖𝑚𝑒 = 𝛥ₘₐₓ / 𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤 ; 𝛥_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒 = 𝛥ₘₐₓ / 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤

	return (_blah=nothing,
			𝛺ₘₐₓ, 𝛺ᵣₑₛ,
			𝛥ₘₐₓ, 𝛥ᵣₑₛ,
			𝑡ᵒᶠᶠₘₐₓ    = δround_down(
							𝑡ₘₐₓ - max(𝛺_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒,
									   𝛥_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒);
							𝛿=𝑡ᵣₑₛ),
			𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ = δround_up(
							max( 𝛺_𝑢𝑝𝑡𝑖𝑚𝑒,
								 𝛥_𝑢𝑝𝑡𝑖𝑚𝑒,
								 𝛥𝑡ₘᵢₙ);
							𝛿=𝑡ᵣₑₛ),
			𝑡ᵣₑₛ, 𝑡ₘₐₓ)
end


@doc raw"""
Function `get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺(hw::HW_Descr ;  𝛺 =hw.𝛺ₘₐₓ, 𝛥 =hw.𝛥ₘₐₓ) `

𝛺-pulse must end this quantity *later* than 𝛥-pulse in order not to break the RWA with max 𝛺,𝛥
"""
function
get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺( hw ::HW_Descr{ℚ}
					;
					𝛺 :: Rad_per_μs_t{ℚ} = hw.𝛺ₘₐₓ,
					𝛥 :: Rad_per_μs_t{ℚ} = hw.𝛥ₘₐₓ ) ::μs_t{ℚ}

	( ; 𝛺ₘₐₓ, 𝛺ᵣₑₛ, 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤, φᵣₑₛ,
		𝛥ₘₐₓ, 𝛥ᵣₑₛ, 𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
		𝑡ₘₐₓ, 𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ                               ) = hw

	𝛺 = abs(𝛺)
	𝛥 = abs(𝛥)
	@assert 𝛺 > 0/μs

	𝛺_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒 = 𝛺 / 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤
	𝛥_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒 = 𝛥 / 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤

	return δround_up(  abs(𝛥_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒 - 𝛺_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒)
					; 𝛿=𝑡ᵣₑₛ)
end

function ϕeⁱᴴψ(ϕ ::Vector{ℂ}, R ::Hermitian{ℂ,Matrix{ℂ}}, ψ ::Vector{ℂ}
               ;
               𝛺      ::Rad_per_μs_t{ℚ},
               𝛥      ::Rad_per_μs_t{ℚ},

			   Ω_𝑡ᵒⁿ  ::μs_t{ℚ},
               Ω_𝑡ᵒᶠᶠ ::μs_t{ℚ},
			   Δ_𝑡ᵒⁿ  ::μs_t{ℚ},
               Δ_𝑡ᵒᶠᶠ ::μs_t{ℚ},
			   𝑇      ::μs_t{ℚ},

               ε      ::ℝ,

			   hw     ::HW_Descr                                          ) ::ℂ

    @assert length(ϕ) == length(ψ)
    @assert ( length(ϕ) , length(ψ) ) == size(R)


	(; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(hw)
	# (; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ, 𝑡ᵒᶠᶠₘₐₓ,𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ) = get_hw_info(hw)

	pΩ = Pulse__Ω_BangBang{ℚ,ℝ}(Ω_𝑡ᵒⁿ, Ω_𝑡ᵒᶠᶠ, 𝑇, 𝛺
								;   hw.𝛺ₘₐₓ, hw.𝛺ᵣₑₛ,
									hw.𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, hw.𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
									hw.φᵣₑₛ,
									hw.𝑡ₘₐₓ, hw.𝑡ᵣₑₛ, hw.𝛥𝑡ₘᵢₙ)
	DOT_RydSim._check(pΩ)

	pΔ = Pulse__Δ_BangBang{ℚ}(Δ_𝑡ᵒⁿ, Δ_𝑡ᵒᶠᶠ, 𝑇, 𝛥
							  ;  hw.𝛥ₘₐₓ, hw.𝛥ᵣₑₛ,
							     hw.𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, hw.𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
								 hw.𝑡ₘₐₓ, hw.𝑡ᵣₑₛ, hw.𝛥𝑡ₘᵢₙ)
	DOT_RydSim._check(pΔ)

	ψᵤₛₑ = copy(ψ)

	schröd!(  ψᵤₛₑ, ℝ(𝑇)
			  ;
              Ω = pΩ,
			  Δ = pΔ,
			  R,
			  ε )

	return ϕ' ⋅ ψᵤₛₑ
end

evf(ϕ ::Vector{ℂ},
    R ::Hermitian{ℂ,Matrix{ℂ}},
    ψ ::Vector{ℂ}               ; kwargs... ) ::ℝ   = 1 - 2⋅abs²( ϕeⁱᴴψ(ϕ,R,ψ;kwargs...) )
    # ;
    # 𝛺 ::Rad_per_μs_t{ℚ},
    # 𝛥 ::Rad_per_μs_t{ℚ},
	# Ω_𝑡ᵒⁿ, Ω_𝑡ᵒᶠᶠ,
	# Δ_𝑡ᵒⁿ, Δ_𝑡ᵒᶠᶠ,
	# 𝑇,
	# hw                           )

function evf_Ω(𝛺  ::Rad_per_μs_t
               ;
               𝛥  ::Rad_per_μs_t,
               ϕ  ::Vector{ℂ},
               R  ::Hermitian{ℂ,Matrix{ℂ}},
               ψ  ::Vector{ℂ},
               ε  ::ℝ,
               hw ::HW_Descr              ) ::ℝ

	(; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ, 𝑡ᵒᶠᶠₘₐₓ,𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ, 𝑡ᵣₑₛ,𝑡ₘₐₓ) = get_hw_data(hw)
	(;𝛥𝑡ₘᵢₙ) = hw

	evf(ϕ,R,ψ ; 𝛺, 𝛥,
	 	Ω_𝑡ᵒⁿ=𝛥𝑡ₘᵢₙ, Ω_𝑡ᵒᶠᶠ=𝑡ᵒᶠᶠₘₐₓ,
	 	Δ_𝑡ᵒⁿ=𝛥𝑡ₘᵢₙ, Δ_𝑡ᵒᶠᶠ=𝑡ᵒᶠᶠₘₐₓ-get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺(hw;𝛺=𝛺ₘₐₓ, 𝛥),
	 	𝑇=𝑡ₘₐₓ,
        ε,
	 	hw)
end

function evf_Δ(𝛥  ::Rad_per_μs_t
               ;
               𝛺  ::Rad_per_μs_t,
               ϕ  ::Vector{ℂ},
               R  ::Hermitian{ℂ,Matrix{ℂ}},
               ψ  ::Vector{ℂ},
               ε  ::ℝ,
               hw ::HW_Descr              ) ::ℝ
    (; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ, 𝑡ᵒᶠᶠₘₐₓ,𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ, 𝑡ᵣₑₛ,𝑡ₘₐₓ) = get_hw_data(hw)
	(;𝛥𝑡ₘᵢₙ) = hw

    evf(ϕ,R,ψ ; 𝛺, 𝛥,
	 	Ω_𝑡ᵒⁿ=𝛥𝑡ₘᵢₙ, Ω_𝑡ᵒᶠᶠ=𝑡ᵒᶠᶠₘₐₓ+get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺(hw;𝛺, 𝛥=𝛥ₘₐₓ),
	 	Δ_𝑡ᵒⁿ=𝛥𝑡ₘᵢₙ, Δ_𝑡ᵒᶠᶠ=𝑡ᵒᶠᶠₘₐₓ,
	 	𝑇=𝑡ₘₐₓ,
        ε,
	 	hw)
end


end #^ module DOT_RydSimDeriv
# EOF
