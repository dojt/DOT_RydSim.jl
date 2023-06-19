using .DOT_RydSimDeriv

using DOT_NiceMath
using DOT_NiceMath.NumbersF64


using Unitful
using Unitful: 	Time, Frequency,
				μs

using Plots
plotly();


using LinearAlgebra: eigvals, Hermitian, normalize


hw = load_hw(;  Ω_downslew_factor = 1//3,
				Δ_downslew_factor = 1//2)


N_ATOMS  = 1
R_STDDEV = 64
LOGε     = -1.0

ψ = randn(ℂ,2^N_ATOMS) |> normalize

R = let A = randn(ℂ,2^N_ATOMS,2^N_ATOMS) ; Hermitian( (A+A')⋅R_STDDEV/2 ) end

println("λ⃗ = ", eigvals(R) )

let
    global hw
	(; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(hw)
	𝛥 = 𝛥ᵣₑₛ #-𝛥ₘₐₓ/2
	scatter( 𝛺 -> evf_Ω(𝛺;𝛥)|>ℜ , -𝛺ₘₐₓ: 7𝛺ᵣₑₛ :+𝛺ₘₐₓ
			; label="",
			markersize=0.5, markerstrokewidth=0,
			xaxis="𝛺")
	scatter!(𝛺 -> evf_Ω(𝛺;𝛥)|>ℑ , -𝛺ₘₐₓ: 7𝛺ᵣₑₛ :+𝛺ₘₐₓ
			; label="",
			markersize=0.5, markerstrokewidth=0,
			xaxis="𝛺")
end

let
	(; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(hw)
	𝛺 = -𝛺ₘₐₓ/100
	scatter( 𝛥 -> evf_Δ(𝛥; 𝛺)|>ℜ , -𝛥ₘₐₓ: 100001𝛥ᵣₑₛ :+𝛥ₘₐₓ
			; label="",
			markersize=0.5, markerstrokewidth=0,
			xaxis="𝛥")
	scatter!(𝛥 -> evf_Δ(𝛥; 𝛺)|>ℑ , -𝛥ₘₐₓ: 100001𝛥ᵣₑₛ :+𝛥ₘₐₓ
			; label="",
			markersize=0.5, markerstrokewidth=0,
			xaxis="𝛥")
end
