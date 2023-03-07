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
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 0. Packages

using Test
using Logging

using Unitful
using Unitful: μs

#
using DOT_RydSim

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1. Tests


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.0: dummy()
module Test__Dummy
export test__dummy
using Test
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


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.1: Schrödinger
module Test__Schrödinger
export test__Schrödinger

using Test
using Logging
using Unitful: μs
using LinearAlgebra: Hermitian
using GenericLinearAlgebra

using DOT_NiceMath
import DOT_NiceMath.NumbersF64
import DOT_NiceMath.NumbersBig

using DOT_RydSim.Schrödinger
using DOT_RydSim: μs_t, Rad_per_μs_t, Radperμs_per_μs_t

expi(A::Hermitian) = cis(A)

δ(x,y) = ( x==y    ? 1 : 0 )
N(x,y) = ( x==2==y ? 1 : 0 )
X(x,y;γ) = if     x==1 && y==2    γ
           elseif x==2 && y==1    γ'
           else                   zero(typeof(γ))  end

function test__Schrödinger(Opts::Symbol...)
    ALL_OPTS = [ :Big ]
    Opts ⊆ ALL_OPTS  ||  throw(
        ArgumentError("Options not recognized: $(setdiff(Opts,ALL_OPTS))")
    )

    if :Big ∈ Opts
        ℝ  = NumbersBig.ℝ
        ℂ  = NumbersBig.ℂ
        𝒊  = NumbersBig.𝒊
        𝒊π = NumbersBig.𝒊
    else
        ℝ  = NumbersF64.ℝ
        ℂ  = NumbersF64.ℂ
        𝒊  = NumbersBig.𝒊
        𝒊π = NumbersBig.𝒊
    end

	@testset verbose=true """Sub-module `Schrödinger` $(:Big∈Opts ? "(w/ BigFloat)" : "")""" begin

        @testset "Helpers" begin
            let N1 = Schrödinger.N₁(1,ℂ)
                @test N1 isa Hermitian{ℂ,Matrix{ℂ}}
                for k=1:2
                    for ℓ=1:2
                        @test N1[k,ℓ] isa ℂ
                        @test N1[k,ℓ] == N(k,ℓ)
                    end
                end
            end #^ N(1)
            let N2 = Schrödinger.N₁(2,ℂ)
                @test N2 isa Hermitian{ℂ,Matrix{ℂ}}
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
        end

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
                                          X,N,R)
                    @test ψ ≈ expi(-Δt⋅(ω⋅X - δ⋅N + R))⋅ψ₀
                end
            end
        end #^ for N

        @testset verbose=true "schröd!()" begin
            @testset "Let's just run it!" begin
                𝟐ᴬ = 8
                ψ ::Vector{ℂ} = randn(𝟐ᴬ)
                𝑇 ::μs_t{ℝ}   = 1μs
                γ ::ℂ         = exp( 𝒊π ⋅ 1.01 )
                R             = Matrix{ℂ}( randn(𝟐ᴬ,𝟐ᴬ) ) |> Hermitian

                function 𝜔( 𝑡 ::μs_t{REAL} ) ::Rad_per_μs_t{REAL} where{REAL}
                    0.9876/μs
                end
                function 𝛿( 𝑡 ::μs_t{REAL} ) ::Rad_per_μs_t{REAL} where{REAL}
                    0.6789/μs
                end
                @testset "Self-test:" begin
                    @test 𝜔(1.1μs) == 0.9876/μs
                    @test 𝛿(1.1μs) == 0.6789/μs
                end

                @test schröd!(ψ,𝑇,γ ; 𝜔, 𝛿, R) === nothing
            end
        end
    end
end #^ test__Schrödinger()
end #^ module Test__Schrödinger
using .Test__Schrödinger

@testset verbose=true "Testing DOT_RydSim.jl" begin
    test__units()
    test__Schrödinger()
    test__Schrödinger(:Big)
    @testset "A broken test:" begin
        @test fasle skip=true
    end
end

#runtests.jl
#EOF
