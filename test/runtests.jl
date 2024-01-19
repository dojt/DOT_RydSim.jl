########################################################################
#                                                                      #
# DOT_RydSim/test/runtests.jl                                          #
#                                                                      #
# (c) Dirk Oliver Theis 2023                                           #
#                                                                      #
# License:                                                             #
#                                                                      #
#             Apache 2.0                                               #
#                                                                      #
########################################################################

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 0. Packages & Helpers

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 0.1. Packages

using Test
using JET

using Logging

using Unitful
using Unitful: μs


using DOT_RydSim


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 0.2. ......


# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1. Tests


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.0: dummy()
module Test__Dummy
export test__dummy
using Test
using JET
using Logging
using DOT_RydSim

function test__dummy(option_set::Symbol...)
    ALL_OPTS = [ ]
    option_set ⊆ ALL_OPTS  || throw(ArgumentError("Options not recognized: $(setdiff(option_set,ALL_OPTS))"))

	@testset verbose=true "Dummy" begin
        @testset "......." begin
        end
    end
end #^ test__dummy()
end #^ module Test__Dummy
using .Test__Dummy

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.1: Units
function test__units(option_set::Symbol...)
    ALL_OPTS = [ ]
    option_set ⊆ ALL_OPTS  || throw(ArgumentError("Options not recognized: $(setdiff(option_set,ALL_OPTS))"))

	@testset verbose=true "Testing the unit types" begin
        @testset verbose=true "μs_t" begin
            @test                  1.0μs    isa DOT_RydSim.μs_t{Float64}
            @test                  (1//2)μs isa DOT_RydSim.μs_t{Rational{Int}}
            @test_throws TypeError (1//2)μs isa DOT_RydSim.μs_t{Char}
        end
        @testset verbose=true "Rad_per_μs_t" begin
            @test                  0.5/μs isa DOT_RydSim.Rad_per_μs_t{Float64}
            @test                  1/2μs  isa DOT_RydSim.Rad_per_μs_t{Float64}
            @test                  1//2μs isa DOT_RydSim.Rad_per_μs_t{Rational{Int}}
            @test_throws TypeError 1/2μs  isa DOT_RydSim.Rad_per_μs_t{String}
        end
        @testset "......." begin
        end
    end
end #^ test__units()



# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.2: Pulses
module Test__Pulses
export test__pulses

using Test
using JET
using Logging
using Unitful: μs
using LinearAlgebra: Hermitian
using GenericLinearAlgebra
using Unitful
using QuadGK

using  DOT_NiceMath
import DOT_NiceMath.NumbersF64
import DOT_NiceMath.NumbersBig

using DOT_RydSim
using DOT_RydSim: δround,
                  Pulse, phase, 𝑎𝑣𝑔, 𝑠𝑡𝑒𝑝
using DOT_RydSim.HW_Descriptions

function _test_pulse(p::Pulse, 𝑇 ::DOT_RydSim.μs_t{𝐑}) ::Nothing    where{𝐑}

    test_opt(   𝑎𝑣𝑔, (typeof(p),typeof(𝑇))  )
    @test_call  𝑎𝑣𝑔(p,𝑇;𝛥𝑡=𝑇)

    for _iter = 1:100
        ≈(a,b) = isapprox(a,b;rtol=1e-2)

        𝑡  = rand()⋅𝑇
        𝛥𝑡 = rand()⋅(𝑇-𝑡)
        𝑎  = 𝑎𝑣𝑔(p,𝑡;𝛥𝑡)
        ∫  = quadgk( t->p(t), 𝑡,𝑡+𝛥𝑡 )[1] / 𝛥𝑡
        @test 1/μs+𝑎 ≈ 1/μs+∫
    end

    for _iter = 1:10
        𝑡  = rand()⋅𝑇
        ε  = rand()⋅ustrip(μs,𝑇-𝑡)⋅1e-9
        𝛥𝑡 = 𝑠𝑡𝑒𝑝(p,𝑡;ε)
        @test 𝛥𝑡 > 0μs
        𝑎  = 𝑎𝑣𝑔(p,𝑡;𝛥𝑡)
        ∫  = quadgk( t -> p(t)-𝑎, 𝑡,𝑡+𝛥𝑡 )[1]
        @test abs( ∫ ) ≤ ε
    end

    nothing;
end

function test__pulses(Opts::Symbol...)
    ALL_OPTS = [ :Big ]
    Opts ⊆ ALL_OPTS  ||  throw(
        ArgumentError("Options not recognized: $(setdiff(Opts,ALL_OPTS))")
    )

    My_Numbers = ( :Big ∈ Opts ? NumbersBig : NumbersF64 )

    ℝ  = My_Numbers.ℝ
    ℂ  = My_Numbers.ℂ
    𝒊  = My_Numbers.𝒊
    𝒊π = My_Numbers.𝒊
    ℤ  = My_Numbers.ℤ
    ℚ  = My_Numbers.ℚ

    @testset verbose=true """Test pulses $(:Big∈Opts ? "(w/ BigFloat)" : "")""" begin
        @testset verbose=true "Helpers" begin
            @testset "δround{,_up,_down,_to0}()" begin

                @test_call δround(0.123;δ=1//31)
                @test_call δround(0.123μs;𝛿=(1//31)μs)
                @test_call δround((123//100)μs;𝛿=(1//31)μs)
                for i = 1:100
                    δ = ℚ(  rationalize(Int16,abs(randn()))  )
                    x = rand(-3:+3)⋅δ
                    x̃ = x + (rand()-0.5)⋅(1-1e-5)⋅δ
                    @test δround(x̃;δ) == x
                end
                for i = 1:100
                    𝛿 = ℚ(  rationalize(Int16,abs(randn()))  )⋅u"kg*100m/s"
                    𝑥 = rand(-3:+3)⋅𝛿
                    𝑥̃ = 𝑥 + (rand()-0.5)⋅(1-1e-5)⋅𝛿
                    @test δround(𝑥̃;𝛿) == 𝑥
                end

                @test_call δround_down(0.123μs;𝛿=(1//31)μs)
                @test_call δround_down((123//100)μs;𝛿=(1//31)μs)
                for i = 1:100
                    𝛿 = ℚ(  rationalize(Int16,abs(randn()))  )⋅u"kg*100m/s"
                    𝑥 = rand(-3:+3)⋅𝛿
                    𝑥̃ = 𝑥 + rand()⋅(1-1e-5)⋅𝛿
                    @test δround_down(𝑥̃;𝛿) == 𝑥
                end

                @test_call δround_down(0.123μs;𝛿=(1//31)μs)
                @test_call δround_down((123//100)μs;𝛿=(1//31)μs)
                for i = 1:100
                    𝛿 = ℚ(  rationalize(Int16,abs(randn()))  )⋅u"kg*100m/s"
                    𝑥 = rand(-3:+3)⋅𝛿
                    𝑥̃ = ( 𝑥>0u"kg*100m/s" ?
                              𝑥 + rand()⋅(1-1e-5)⋅𝛿
                            : 𝑥 - rand()⋅(1-1e-5)⋅𝛿 )
                    @test δround_to0(𝑥̃;𝛿) == 𝑥
                end

                @test_call δround_up(0.123μs;𝛿=(1//31)μs)
                @test_call δround_up((123//100)μs;𝛿=(1//31)μs)
                for i = 1:100
                    𝛿 = ℚ(  rationalize(Int16,abs(randn()))  )⋅u"kg*100m/s"
                    𝑥 = rand(-3:+3)⋅𝛿
                    𝑥̃ = 𝑥 - rand()⋅(1-1e-5)⋅𝛿
                    @test δround_up(𝑥̃;𝛿) == 𝑥
                end
            end
        end #^ testset "Helpers"

        hw = nothing
        @testset "Load default HW-description" begin
            # @test_call default_HW_Descr(;ℤ)
            hw = default_HW_Descr(;ℤ)
            @test hw.𝑡ₘₐₓ isa DOT_RydSim.μs_t{ℚ}
            let hw_hr = cheatify_res(hw ; factor=ℚ(1//17))
                @test hw_hr.𝛺ᵣₑₛ           == hw.𝛺ᵣₑₛ/17
                @test hw_hr.𝛥ᵣₑₛ           == hw.𝛥ᵣₑₛ/17
                @test hw_hr.𝛺ₘₐₓ           == hw.𝛺ₘₐₓ
                @test hw_hr.𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 == hw.𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤
                @test hw_hr.𝑡ₘₐₓ           == hw.𝑡ₘₐₓ
                @test hw_hr.𝑡ᵣₑₛ           == hw.𝑡ᵣₑₛ
            end
            let hw_lg = cheatify_𝑡ₘₐₓ(hw ; factor=ℚ(7//1))
                @test hw_lg.𝛺ᵣₑₛ           == hw.𝛺ᵣₑₛ
                @test hw_lg.𝛺ₘₐₓ           == hw.𝛺ₘₐₓ
                @test hw_lg.𝛥ᵣₑₛ           == hw.𝛥ᵣₑₛ
                @test hw_lg.𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤 == hw.𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤
                @test hw_lg.𝑡ₘₐₓ           == 7hw.𝑡ₘₐₓ
                @test hw_lg.𝑡ᵣₑₛ           == hw.𝑡ᵣₑₛ
            end
        end
        𝑇 = hw.𝑡ₘₐₓ

        (;𝛺ₘₐₓ, 𝛺ᵣₑₛ, 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
         𝛥ₘₐₓ, 𝛥ᵣₑₛ,  𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
         φᵣₑₛ,  𝑡ₘₐₓ, 𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ                    ) = hw

        @testset "Δ-BangBang" begin
            p = Pulse__Δ_BangBang{ℚ}(𝑇/10, 9⋅𝑇/10, 𝑇 , -𝛥ₘₐₓ/2;
                                     𝛥ₘₐₓ, 𝛥ᵣₑₛ, 𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
                                     𝑡ₘₐₓ, 𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ)
            @test      DOT_RydSim._check(p)
            @test_call DOT_RydSim._check(p)
            _test_pulse(p,𝑇)
            p = Pulse__Δ_BangBang{ℚ}(𝑇/10, 9⋅𝑇/10, 𝑇 , ℚ(0//1)/μs;
                                     𝛥ₘₐₓ, 𝛥ᵣₑₛ, 𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
                                     𝑡ₘₐₓ, 𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ)
            @test      DOT_RydSim._check(p)
            @test_call DOT_RydSim._check(p)
            _test_pulse(p,𝑇)
        end

        @testset "Ω-BangBang" begin
            p = Pulse__Ω_BangBang{ℚ,ℝ}(𝑇/10, 9⋅𝑇/10, 𝑇 , -𝛺ₘₐₓ/2;
                                       𝛺ₘₐₓ, 𝛺ᵣₑₛ, 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
                                       φᵣₑₛ,
                                       𝑡ₘₐₓ, 𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ)
            @test      DOT_RydSim._check(p)
            @test_call DOT_RydSim._check(p)
            _test_pulse(p,𝑇)
            p = Pulse__Ω_BangBang{ℚ,ℝ}(𝑇/10, 9⋅𝑇/10, 𝑇 , ℚ(0//1)/μs;
                                       𝛺ₘₐₓ, 𝛺ᵣₑₛ, 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
                                       φᵣₑₛ,
                                       𝑡ₘₐₓ, 𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ)
            @test      DOT_RydSim._check(p)
            @test_call DOT_RydSim._check(p)
            _test_pulse(p,𝑇)
        end

    end #^ testset
end #^ test__pulses()
end #^ module Test__Pulses
using .Test__Pulses: test__pulses

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.3: Schrödinger
module Test__Schrödinger
export test__secrets, test__schröd!

using Test
using JET
using Logging
using Unitful: μs
using LinearAlgebra: Hermitian
using GenericLinearAlgebra

using DOT_NiceMath
import DOT_NiceMath.NumbersF64
import DOT_NiceMath.NumbersBig

using DOT_RydSim
using DOT_RydSim: μs_t, Rad_per_μs_t, Radperμs_per_μs_t
using DOT_RydSim.Schrödinger
using DOT_RydSim.HW_Descriptions

expi(A::Hermitian) = cis(A)

δ(x,y) = ( x==y    ? 1 : 0 )
N(x,y) = ( x==2==y ? 1 : 0 )
X(x,y;γ) = if     x==1 && y==2    γ
           elseif x==2 && y==1    γ'
           else                   zero(typeof(γ))  end

function test__secrets(Opts::Symbol...)
    ALL_OPTS = [ :Big ]
    Opts ⊆ ALL_OPTS  ||  throw(
        ArgumentError("Options not recognized: $(setdiff(Opts,ALL_OPTS))")
    )

    My_Numbers = ( :Big ∈ Opts ? NumbersBig : NumbersF64 )

    ℝ  = My_Numbers.ℝ
    ℂ  = My_Numbers.ℂ
    𝒊  = My_Numbers.𝒊
    𝒊π = My_Numbers.𝒊
    ℤ  = My_Numbers.ℤ
    ℚ  = My_Numbers.ℚ

	@testset verbose=true """Sub-module `Schrödinger`: Secrets $(:Big∈Opts ? "(w/ BigFloat)" : "")""" begin
        @testset "Helpers" begin
            let N1 = Schrödinger.N₁(1,ℂ)
                @test N1 isa Hermitian{ℂ,Matrix{ℂ}}
#                @test_opt  Schrödinger.N₁(1,ℂ)
#                @test_call Schrödinger.N₁(1,ℂ)
                for k=1:2
                    for ℓ=1:2
                        @test N1[k,ℓ] isa ℂ
                        @test N1[k,ℓ] == N(k,ℓ)
                    end
                end
            end #^ N(1)
            let N2 = Schrödinger.N₁(2,ℂ)
                @test N2 isa Hermitian{ℂ,Matrix{ℂ}}
#                @test_opt  Schrödinger.N₁(2,ℂ)
#                @test_call Schrödinger.N₁(2,ℂ)
                for k₁=1:2
                    for k₂=1:2
                        k = 1+ 2(k₁-1)+(k₂-1)
                        for ℓ₁=1:2
                            for ℓ₂=1:2
                                ℓ = 1+ 2(ℓ₁-1)+(ℓ₂-1)
                                @test N2[k,ℓ] isa ℂ
                                @test N2[k,ℓ] == N(k₁,ℓ₁)⋅δ(k₂,ℓ₂)
                            end
                        end
                    end
                end
            end #^ N(2)
            let N3 = Schrödinger.N₁(3,ℂ)
                @test N3 isa Hermitian{ℂ,Matrix{ℂ}}
#                @test_opt  Schrödinger.N₁(3,ℂ)
#                @test_call Schrödinger.N₁(3,ℂ)
                for k₁=1:2
                    for k₂=1:2
                        for k₃=1:2
                            k = 1+ 4(k₁-1)+2(k₂-1)+(k₃-1)
                            for ℓ₁=1:2
                                for ℓ₂=1:2
                                    for ℓ₃=1:2
                                        ℓ = 1+ 4(ℓ₁-1)+2(ℓ₂-1)+(ℓ₃-1)
                                        @test N3[k,ℓ] isa ℂ
                                        @test N3[k,ℓ] == N(k₁,ℓ₁)⋅δ(k₂,ℓ₂)⋅δ(k₃,ℓ₃)
                                    end
                                end
                            end
                        end
                    end
                end
            end #^ N(3)

            let γ::ℂ  = randn(ComplexF64),
                X1    = Schrödinger.X₁(1;γ)
                @test X1 isa Hermitian{ℂ,Matrix{ℂ}}
#                @test_opt  Schrödinger.X₁(1;γ)
#                @test_call Schrödinger.X₁(1;γ)
                for k=1:2
                    for ℓ=1:2
                        @test X1[k,ℓ] isa ℂ
                        @test X1[k,ℓ] == X(k,ℓ;γ)
                    end
                end
            end #^ X(1)=#
            let γ::ℂ  = randn(ComplexF64),
                X2    = Schrödinger.X₁(2;γ)
                @test X2 isa Hermitian{ℂ,Matrix{ℂ}}
#                @test_opt  Schrödinger.X₁(2;γ)
#                @test_call Schrödinger.X₁(2;γ)
                for k₁=1:2
                    for k₂=1:2
                        k = 1+ 2(k₁-1)+(k₂-1)
                        for ℓ₁=1:2
                            for ℓ₂=1:2
                                ℓ = 1+ 2(ℓ₁-1)+(ℓ₂-1)
                                @test X2[k,ℓ] isa ℂ
                                @test X2[k,ℓ] == X(k₁,ℓ₁;γ)⋅δ(k₂,ℓ₂)
                            end
                        end
                    end
                end
            end #^ X(2)
            let γ::ℂ  = randn(ComplexF64),
                X3 = Schrödinger.X₁(3;γ)
                @test X3 isa Hermitian{ℂ,Matrix{ℂ}}
#                @test_opt  Schrödinger.X₁(3;γ)
#                @test_call Schrödinger.X₁(3;γ)
                for k₁=1:2
                    for k₂=1:2
                        for k₃=1:2
                            k = 1+ 4(k₁-1)+2(k₂-1)+(k₃-1)
                            for ℓ₁=1:2
                                for ℓ₂=1:2
                                    for ℓ₃=1:2
                                        ℓ = 1+ 4(ℓ₁-1)+2(ℓ₂-1)+(ℓ₃-1)
                                        @test X3[k,ℓ] isa ℂ
                                        @test X3[k,ℓ] == X(k₁,ℓ₁;γ)⋅δ(k₂,ℓ₂)⋅δ(k₃,ℓ₃)
                                    end
                                end
                            end
                        end
                    end
                end
            end #^ X(3)
        end #^ testset "Helpers"

        @testset "timestep!()" begin
            for 𝟐ᴬ in [2,4,8]
                ψ₀ = Vector{ℂ}( randn(ComplexF64,𝟐ᴬ) )
                X  = let M::Matrix{ℂ}=randn(ComplexF64,𝟐ᴬ,𝟐ᴬ) ; (M+M')/2 |> Hermitian end
                N  = let M::Matrix{ℂ}=randn(ComplexF64,𝟐ᴬ,𝟐ᴬ) ; (M+M')/2 |> Hermitian end
                R  = let M::Matrix{ℂ}=randn(ComplexF64,𝟐ᴬ,𝟐ᴬ) ; (M+M')/2 |> Hermitian end
                for _iter = 1:3
                    ψ     = copy(ψ₀)
                    Δt::ℝ = rand()
                    ω ::ℝ = rand()
                    δ ::ℝ = rand()
                    Schrödinger.timestep!(ψ, Δt⋅μs
                                          ;
                                          𝜔=ω/μs, 𝛿=δ/μs,
                                          X_2=X/2,N,R)
                    @test ψ ≈ expi(-Δt⋅(ω⋅X/2 - δ⋅N + R))⋅ψ₀
#                    @test_call Schrödinger.timestep!(ψ, Δt⋅μs
#                                                     ;
#                                                     𝜔=ω/μs, 𝛿=δ/μs,
#                                                     X_2=X/2,N,R)
                end
            end
        end #^ tstset "timestep!()"
    end #^ function-testset
end #^ test__secrets()

function test__schröd!(Opts::Symbol...)
    ALL_OPTS = [ :Big ]
    Opts ⊆ ALL_OPTS  ||  throw(
        ArgumentError("Options not recognized: $(setdiff(Opts,ALL_OPTS))")
    )

    My_Numbers = ( :Big ∈ Opts ? NumbersBig : NumbersF64 )

    ℝ  = My_Numbers.ℝ
    ℂ  = My_Numbers.ℂ
    𝒊  = My_Numbers.𝒊
    𝒊π = My_Numbers.𝒊
    ℤ  = My_Numbers.ℤ
    ℚ  = My_Numbers.ℚ

    @testset verbose=true """Sub-module `Schrödinger`: schröd!() $(:Big∈Opts ? "(w/ BigFloat)" : "")""" begin

        @testset verbose=true "schröd!()" begin
            @testset "Let's just run it!" begin
                𝟐ᴬ = 8
                ψ ::Vector{ℂ}              = randn(𝟐ᴬ)
                R ::Hermitian{ℂ,Matrix{ℂ}} = let A=randn(𝟐ᴬ,𝟐ᴬ) ; Hermitian((A'+A)/2) end

                (;Ω,Δ,𝑇) = let hw = default_HW_Descr(;ℤ)
                    𝑇 = hw.𝑡ₘₐₓ
                    (;𝛺ₘₐₓ, 𝛺ᵣₑₛ, 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
                     𝛥ₘₐₓ, 𝛥ᵣₑₛ,  𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
                     φᵣₑₛ,  𝑡ₘₐₓ, 𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ                    ) = hw

                    Δ = Pulse__Δ_BangBang{ℚ}(𝑇/10, 9⋅𝑇/10, 𝑇 , -𝛥ₘₐₓ/2;
                                             𝛥ₘₐₓ, 𝛥ᵣₑₛ, 𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
                                             𝑡ₘₐₓ, 𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ)
                    Ω = Pulse__Ω_BangBang{ℚ,ℝ}(𝑇/10, 9⋅𝑇/10, 𝑇 , -𝛺ₘₐₓ/2;
                                               𝛺ₘₐₓ, 𝛺ᵣₑₛ, 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
                                               φᵣₑₛ,
                                               𝑡ₘₐₓ, 𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ)
                    (Ω=Ω,Δ,𝑇)
                end

                # @test      schröd!(ψ,𝑇 ; Ω, Δ, R) === nothing       skip=true
#                @test_call schröd!(ψ,𝑇 ; Ω, Δ, R)
            end
        end #^ testset "schröd!()"
    end #^ function-testset
end #^ test__schröd!

end #^ module Test__Schrödinger
import .Test__Schrödinger

# ——————————————————————————————————————————————————————————————————————————————————————————————————— X. Main



using JSON # Only for ignoring by JET

@testset verbose=true "Testing DOT_RydSim.jl" begin
    test__units()
    test__pulses()
    test__pulses(:Big)
    Test__Schrödinger.test__secrets()
    Test__Schrödinger.test__secrets(:Big)
    Test__Schrödinger.test__schröd!()
    Test__Schrödinger.test__schröd!(:Big)

    #
    # Basic JET-based package test:

#   Stupid shit doesn't work.  Why can't these moronic millennials declare variables?!
#   test_package(DOT_RydSim, ignored_modules=(AnyFrameModule(JSON.Parser),) )

end

#  @testset "A broken test:" begin
#      @test DOODELDIDOO skip=true
#  end

#runtests.jl
#EOF
