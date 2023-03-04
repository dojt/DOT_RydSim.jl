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
using DOT_NiceMath.NumbersBig

using DOT_RydSim.Schrödinger

expi(A::Hermitian) = cis(A)

function test__Schrödinger(option_set::Symbol...)
    ALL_OPTS = [ ]
    option_set ⊆ ALL_OPTS  || throw(ArgumentError("Options not recognized: $(setdiff(option_set,ALL_OPTS))"))

	@testset verbose=true "Sub-module `Schrödinger`" begin
        @testset "Sub-module loaded :)" begin end

        @testset "timestep!()" begin
            ψ₀ = randn(ComplexF64,8)
            X = let A=randn(ComplexF64,8,8) ; (A+A')/2 |> Hermitian end
            Z = let A=randn(ComplexF64,8,8) ; (A+A')/2 |> Hermitian end
            R = let A=randn(ComplexF64,8,8) ; (A+A')/2 |> Hermitian end
            for _iter = 1:3
                ψ     = copy(ψ₀)
                Δt::ℝ = rand()
                ω ::ℝ = rand()
                δ ::ℝ = rand()
                Schrödinger.timestep!(ψ, Δt⋅μs
                                      ;
                                      𝜔=ω/μs, 𝛿=δ/μs,
                                      X,Z,R)
                @test ψ ≈ expi(-Δt⋅(ω⋅X + δ⋅Z + R))⋅ψ₀
            end
        end

        @testset "......." begin
            @test false skip=true
        end
    end
end #^ test__Schrödinger()
end #^ module Test__Schrödinger
using .Test__Schrödinger

@testset verbose=true "Testing DOT_RydSim.jl" begin
    test__units()
    test__Schrödinger()
end

#runtests.jl
#EOF
