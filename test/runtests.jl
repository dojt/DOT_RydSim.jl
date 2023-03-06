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
        ℝ = NumbersBig.ℝ
        ℂ = NumbersBig.ℂ
    else
        ℝ = NumbersF64.ℝ
        ℂ = NumbersF64.ℂ
    end

	@testset verbose=true """Sub-module `Schrödinger` $(:Big∈Opts ? "(w/ BigFloat)" : "")""" begin

        @testset "Helpers" begin
            let N₁ = Schrödinger.N(1,ℂ)
                @test N₁ isa Hermitian{ℂ,Matrix{ℂ}}
                for k=1:2
                    for ℓ=1:2
                        @test N₁[k,ℓ] isa ℂ
                        @test N₁[k,ℓ] == N(k,ℓ)
                    end
                end
            end #^ N(1)
            let N₂ = Schrödinger.N(2,ℂ)
                @test N₂ isa Hermitian{ℂ,Matrix{ℂ}}
                for k₁=1:2
                    for k₂=1:2
                        k = 1+ 2(k₁-1)+(k₂-1)
                        for ℓ₁=1:2
                            for ℓ₂=1:2
                                ℓ = 1+ 2(ℓ₁-1)+(ℓ₂-1)
                                @test N₂[k,ℓ] isa ℂ
                                @test N₂[k,ℓ] == N(k₁,ℓ₁)⋅δ(k₂,ℓ₂)
                            end
                        end
                    end
                end
            end #^ N(2)
            let N₃ = Schrödinger.N(3,ℂ)
                @test N₃ isa Hermitian{ℂ,Matrix{ℂ}}
                for k₁=1:2
                    for k₂=1:2
                        for k₃=1:2
                            k = 1+ 4(k₁-1)+2(k₂-1)+(k₃-1)
                            for ℓ₁=1:2
                                for ℓ₂=1:2
                                    for ℓ₃=1:2
                                        ℓ = 1+ 4(ℓ₁-1)+2(ℓ₂-1)+(ℓ₃-1)
                                        @test N₃[k,ℓ] isa ℂ
                                        @test N₃[k,ℓ] == N(k₁,ℓ₁)⋅δ(k₂,ℓ₂)⋅δ(k₃,ℓ₃)
                                    end
                                end
                            end
                        end
                    end
                end
            end #^ N(3)

            let γ::ℂ  = randn(ComplexF64),
                X₁    = Schrödinger.X(1,γ)
                @test X₁ isa Hermitian{ℂ,Matrix{ℂ}}
                for k=1:2
                    for ℓ=1:2
                        @test X₁[k,ℓ] isa ℂ
                        @test X₁[k,ℓ] == X(k,ℓ;γ)
                    end
                end
            end #^ X(1)=#
            let γ::ℂ  = randn(ComplexF64),
                X₂    = Schrödinger.X(2,γ)
                @test X₂ isa Hermitian{ℂ,Matrix{ℂ}}
                for k₁=1:2
                    for k₂=1:2
                        k = 1+ 2(k₁-1)+(k₂-1)
                        for ℓ₁=1:2
                            for ℓ₂=1:2
                                ℓ = 1+ 2(ℓ₁-1)+(ℓ₂-1)
                                @test X₂[k,ℓ] isa ℂ
                                @test X₂[k,ℓ] == X(k₁,ℓ₁;γ)⋅δ(k₂,ℓ₂)
                            end
                        end
                    end
                end
            end #^ X(2)
            let γ::ℂ  = randn(ComplexF64),
                X₃ = Schrödinger.X(3,γ)
                @test X₃ isa Hermitian{ℂ,Matrix{ℂ}}
                for k₁=1:2
                    for k₂=1:2
                        for k₃=1:2
                            k = 1+ 4(k₁-1)+2(k₂-1)+(k₃-1)
                            for ℓ₁=1:2
                                for ℓ₂=1:2
                                    for ℓ₃=1:2
                                        ℓ = 1+ 4(ℓ₁-1)+2(ℓ₂-1)+(ℓ₃-1)
                                        @test X₃[k,ℓ] isa ℂ
                                        @test X₃[k,ℓ] == X(k₁,ℓ₁;γ)⋅δ(k₂,ℓ₂)⋅δ(k₃,ℓ₃)
                                    end
                                end
                            end
                        end
                    end
                end
            end #^ X(3)
        end

        @testset "timestep!()" begin
            for N in [2,4,8]
                ψ₀ = Vector{ℂ}( randn(ComplexF64,N) )
                X  = let A::Matrix{ℂ}=randn(ComplexF64,N,N) ; (A+A')/2 |> Hermitian end
                Z  = let A::Matrix{ℂ}=randn(ComplexF64,N,N) ; (A+A')/2 |> Hermitian end
                R  = let A::Matrix{ℂ}=randn(ComplexF64,N,N) ; (A+A')/2 |> Hermitian end
                for _iter = 1:3
                    ψ     = copy(ψ₀)
                    Δt::ℝ = rand()
                    ω ::ℝ = rand()
                    δ ::ℝ = rand()
                    Schrödinger.timestep!(ψ, Δt⋅μs
                                          ;
                                          𝜔=ω/μs, 𝛿=δ/μs,
                                          X,Z,R)
                    @test ψ ≈ expi(-Δt⋅(ω⋅X - δ⋅Z + R))⋅ψ₀
                end
            end
        end #^ for N

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
