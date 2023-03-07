########################################################################
#                                                                      #
# DOT_RydSim/src/Schrödinger.mod.jl                                    #
#                                                                      #
# (c) Dirk Oliver Theis 2023                                           #
#                                                                      #
# License:                                                             #
#                                                                      #
#             Apache 2.0                                               #
#                                                                      #
########################################################################


"""
Module `Schrödinger`

Simulation of multi-qubit quantum evolution under time-dependent Hamiltonian.

# Exports

## General
* (nothing yet)

## Number definitions (in sub-modules `Numbers`𝑥𝑦𝑧)
* (also nothing)

# Sub-modules

Sub-module names are not exported by `DOT_RydSim`.

* `DOT_NiceMath.` — nope...
"""
module Schrödinger
export schröd!

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 0. TOC  +
#                                                                                                             |
#  Table of Contents                                                                                          |
#  —————————————————                                                                                          |
#                                                                                                             |
#                                                                                                             |
#  1. Imports & Helpers                                                                                       |
#                                                                                                             |
#     1.1. X,N matrices                                                                                       |
#     1.1. Shit                                                                                               |
#                                                                                                             |
#                                                                                                             |
#  2. Work horses                                                                                             |
#                                                                                                             |
#     2.1. timestep!()                                                                                        |
#     2.2. schröd!()                                                                                          |
#                                                                                                             |
#—————————————————————————————————————————————————————————————————————————————————————————————————————————————+


# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1. Imports & Helpers

import ..μs_t, ..Rad_per_μs_t, ..Radperμs_per_μs_t
using  ..DOT_NiceMath

using LinearAlgebra: Hermitian, I as Id

using Unitful
using Unitful: μs

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.1: X,Z matrices
#

@doc raw"""
Function `N₁(A ::Integer, ℂ ::Type{<Complex})` — Rydberg |r⟩⟨r| operator for the 1st atom

Letter is customarily `N` as it's kinda a number operator, except of course, here, |g⟩=|1⟩
and |r⟩=|2⟩ — so it ain't exactly no *number* operator.

## Arguments:
1. `A` — number of atoms
2. `ℂ` — complex number type to use as `eltype`.

### Returns:
* 2ᴬ×2ᴬ Hermitian matrix; type `Hermitian{ℂ,Matrix{ℂ}}`
"""
N₁(A ::Integer, ℂ::Type{<:Complex}) ::Hermitian{ℂ,Matrix{ℂ}} = begin
    @assert A ≥ 1
    let 𝙽 = [  ℂ(0)      0
                0       +1  ]

        𝙽 ⊗ Id(2^(A-1)) |> Hermitian
    end
end

(
    X₁(A ::Integer ; γ::ℂ)   ::Hermitian{ℂ,Matrix{ℂ}}
) where{ℂ} =
    begin
        @assert A ≥ 1
        let 𝚇 = [   0     γ
                    γ'    0   ]

            𝚇 ⊗ Id(2^(A-1)) |> Hermitian
        end
    end

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.2: Shit
#

log_of_pow2(x ::Integer) = begin @assert ispow2(x) ; trailing_zeros(x) end


# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2. Work horses

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.1: timestep!()
function  timestep!(ψ  ::Vector{ℂ},
                    𝛥𝑡 ::μs_t{ℝ}
                    ;
                    𝜔  ::Rad_per_μs_t{ℝ},
                    X  ::Hermitian{ℂ,𝕄_t},
                    𝛿  ::Rad_per_μs_t{ℝ},
                    N  ::Hermitian{ℂ,𝕄_t},
                    R  ::Hermitian{ℂ,𝕄_t}            ) ::Nothing      where{ℝ,ℂ,𝕄_t}

    let ωΔt ::ℝ              = ustrip(NoUnits, 𝜔⋅𝛥𝑡),
        δΔt ::ℝ              = ustrip(NoUnits, 𝛿⋅𝛥𝑡),
        Δt  ::ℝ              = ustrip(μs, 𝛥𝑡),
        A   ::Hermitian{ℂ,𝕄_t} = ωΔt⋅X - δΔt⋅N + Δt⋅R

        ψ .= cis(A)'ψ
    end
    nothing
end #^ timestep!()

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.2: schröd!()
@doc raw"""
Function `schröd!(ψ  ::Vector{ℂ},  𝑇 ::μs_t{ℝ},  γ ::ℂ ; ...)   where{ℝ,ℂ,𝕄_t}`

Simulates time evolution under time-dependent Hamiltonian
```math
    H(t)/\hbar = \omega(t) Xᵧ - \delta(t) N + R
```

where `Xᵧ` is defined as |1⟩⟨2|⋅γ + |2⟩⟨1|⋅γ̄.  The phase, γ, is not time-dependent.

The number of atoms, `A`, is derived from the dimensions of ``R``, which is assumed to be
a 2ᴬ×2ᴬ matrix.  The matrices ``X,N`` are constructed to act on the "first" atom, i.e.,
X⊗1⊗...⊗1 and N⊗1⊗...⊗1; see the function `⊗` above.

## Arguments

Note the convention: *Italic variables are unitful*. That also applies to functions: The
returned value of an italic-type function must be unitful.

### Position variables
1. `ψ` — state vector at time 𝑡=0μs; will be modified in place.
2. `𝑇` — end-time of evolution, in μs.
3. `γ` — phase for the X-term; time independent, dimensionless

### Returns
...nothing.

### Keyword arguments

#### Mandatory keyword arguments: Defining the Hamiltonian

* `𝜔` — function which, essentially, gives the Rabi frequency as a function of time.  (The
  Rabi frequency must always be non-negative.)   See § *The Functions* below.
* `𝛿` — function which, essentially, gives the Rabi frequency as a function of time.  See
  § *The Functions* below.
* `R` — Hermitian operators, of type `Hermitian{ℂ,𝕄_t}`

#### Optional keyword arguments
* `𝑚𝑎𝑥_𝜔_𝑠𝑙𝑒𝑤` — (type `::Radperμs_per_μs_t{ℝ}`) defaults to infinity; is checked in every
  time step.
* `𝑚𝑎𝑥_𝛿_𝑠𝑙𝑒𝑤` — (type `::Radperμs_per_μs_t{ℝ}`) defaults to infinity; is checked in every
  time step.

## The Functions: 𝜔, 𝛿
The signatures must be
```julia
    𝜔( 𝑡 ::μs_t{ℝ}, ε ::ℝ )  ::TimeStep_t{ℝ}
    𝛿( 𝑡 ::μs_t{ℝ}, ε ::ℝ )  ::TimeStep_t{ℝ}
```

With ``f`` the Rabi frequency (case `𝜔`) or detuning (case `𝛿`), respectively, on input
``t \in \left[0,T\right[`` and ``\varepsilon > 0``, the functions return, in the struct,
the information:
* `𝛥𝑡` — largest ``\delta \in \left]0,T-t\right]`` (in μs) for which
```math
    \int_t^{t+\delta} | f(s) - μ_\delta|\,ds \le \varepsilon,
```
with
```math
    \mu_δ := \tfrac{1}{\delta} \int_t^{t+\delta} f(s) \,ds.
```
* `𝑎𝑣𝑔` the value of ``\mu_{\Delta\!t}``:
```math
    \tfrac{1}{\Delta\!t} \int_t^{t+\Delta\!t} f(s) \,ds,
```
in rad/μs.
"""
function schröd!(ψ  ::Vector{ℂ},
                 𝑇  ::μs_t{ℝ},
                 γ  ::ℂ
                 ;
                 𝜔  ::Function,
                 𝛿  ::Function,
                 R  ::Hermitian{ℂ,𝕄_t},
                 ε          ::ℝ                    = ℝ(1e-3),
                 𝑚𝑎𝑥_𝜔_𝑠𝑙𝑒𝑤 ::Radperμs_per_μs_t{ℝ} = ℝ(1e50)/μs^2,
                 𝑚𝑎𝑥_𝛿_𝑠𝑙𝑒𝑤 ::Radperμs_per_μs_t{ℝ} = ℝ(1e50)/μs^2  ) ::Nothing  where{ℝ,ℂ,𝕄_t}


    A    = log_of_pow2( length(ψ) )       ; @assert A ≥ 1               "Need at least one atom, i.e., length ψ ≥ 2."
    𝟐ᴬ   = length(ψ)                      ; @assert 2^A == 𝟐ᴬ           "Crazy bug #1"
    N    = N₁(A,ℂ)                        ; @assert size(N) == size(R)  "Sizes of `ψ` and `R` don't match."
    X    = X₁(A;γ)                        ; @assert size(X) == size(N)  "Crazy bug #2"

    𝑡 ::μs_t{ℝ} = 0μs

    while 𝑡  <  𝑇 - 1e-50μs

        Ω_𝛥𝑡 ::μs_t{ℝ} = min(𝑇-𝑡, 𝜔(STEP, 𝑡;ε) )
        Δ_𝛥𝑡 ::μs_t{ℝ} = min(𝑇-𝑡, 𝛿(STEP, 𝑡;ε) )

        let 𝛺𝑠𝑙𝑒𝑤, 𝛥𝑠𝑙𝑒𝑤
            Ω_𝛥𝑡 > 1e-50μs ||
                throw(Ctrl_Exception("Time-step for Ω is non-positive: $(Ω_𝛥𝑡) ≤ 0μs"))

            Δ_𝛥𝑡 > 1e-50μs ||
                throw(Ctrl_Exception("Time-step for Δ is non-positive: $(Δ_𝛥𝑡) ≤ 0μs"))

            𝛺𝑠𝑙𝑒𝑤 = 4ε/Ω_𝛥𝑡^2
            𝛥𝑠𝑙𝑒𝑤 = 4ε/Δ_𝛥𝑡^2

            if 𝛺𝑠𝑙𝑒𝑤 > 𝑚𝑎𝑥_𝜔_𝑠𝑙𝑒𝑤
                throw(Ctrl_Exception("Slew rate for Ω exceeded: $(𝛺𝑠𝑙𝑒𝑤) > $(𝑚𝑎𝑥_𝜔_𝑠𝑙𝑒𝑤))"))
            end
            if 𝛥𝑠𝑙𝑒𝑤 > 𝑚𝑎𝑥_𝛿_𝑠𝑙𝑒𝑤
                throw(Ctrl_Exception("Slew rate for Ω exceeded: $(𝛥𝑠𝑙𝑒𝑤) > $(𝑚𝑎𝑥_𝛿_𝑠𝑙𝑒𝑤))"))
            end
        end

        𝛥𝑡 = min( Ω.𝛥𝑡, Δ.𝛥𝑡 )

        Ω_𝜇 ::Rad_per_μs_t{ℝ} = 𝜔(AVG, 𝑡; 𝛥𝑡)
        Δ_𝜇 ::Rad_per_μs_t{ℝ} = 𝛿(AVG, 𝑡; 𝛥𝑡)


        timestep!(ψ, 𝛥𝑡 ; 𝜔=Ω_𝜇, 𝛿=Δ_𝜇,
                          X, N, R)

        𝑡 += 𝛥𝑡

    end


    nothing
end #^ schröd!()


end # module Schrödinger

#
# Thoughts about the algorithm
#
# Matrix exponential
#
# For matrices up to 4×4:       cis( StaticMatrix )           ~1 μs  (vs ~2μs)
#
# For matrices from 16×16       cis( Hermitian(Matrix) )     ~26 μs (vs ∞)
#
#
#EOF
