# DOT_RydSim.jl
Quantum simulation of arrays of Rydberg atoms &amp; related stuff

## Version History

Time goes up!


####  **v0.3.0**

* Added functions`𝔛()` and `𝔑()` (for X- and N- operators)


####  **v0.2.0**

* Added function `δround_to0()`

####  **v0.1.29**

* Fixed bug in `HW_Descriptions.default_HW_Descr` options handling

####  **v0.1.28**

* Added options to fn `HW_Descriptions.default_HW_Descr()`
* ... currently the only one: `:hires`, which loads 1/1000 resolution of Ω, Δ from new file "Resources/hw_hires.json"

####  **v0.1.27**

* Fixed bug in `is_δrounded()` -- ARGH!!!

####  **v0.1.26**

* Added exported fn `is_δrounded()`

####  **v0.1.25**

* Fixed missing factor 1/2 on Rabi frequency

####  **v0.1.24**

* `schröd!()` now has 𝑡₀ optional keyword arg

####  **v0.1.23**

* More detailed error checks in pulse constructor

####  **v0.1.22**

* Fixed bug w/ γ=0 in `Pulse__Ω_BangBang`

####  **v0.1.21**

* Fixed bug w/ 𝛥𝑡ₘᵢₙ

####  **v0.1.20**

* Fixed bug w/ `rationalize( x::Rational )`

####  **v0.1.19**

* Added `δround_down()`, `δround_up()`
####  **v0.1.18**

* Rotating-Wave-Approx breakdown is now thrown as `Ctrl_Exception`
* Improved docs in `HW_Descriptions`

####  **v0.1.17**

* Added `JET.jl` in tests
* Bug fix

####  **v0.1.16**

* Attempt to speed up timestep!() by adding memory workspace

####  **v0.1.15**

* Attempt to speed up comparisons between μs_t{} for float/rational

####  **v0.1.14**

* Fixed 💩

####  **v0.1.13**

* Improved (I hope) warning for Δ≠0 & Ω=0

####  **v0.1.12**

* Added simple check for Δ≠0 & Ω=0
* Improved file-docs for DOT_RydSim.jl

####  **v0.1.11**

**Done**

* 𝑠𝑡𝑒𝑝() fns seem to work now

####  **v0.1.10**

**WIP**

* Try to fix bugs in shapes 𝑠𝑡𝑒𝑝() fn.

####  **v0.1.9**

**WIP**

* Bug fixes in pulse shapes 𝑠𝑡𝑒𝑝() fn.

####  **v0.1.8**

**WIP**

* Bug fix in `Schrödinger.schröd!()` -- again!

####  **v0.1.7**

**Done**

* Bug fix in `Schrödinger.schröd!()`

####  **v0.1.6**

**Done**

* Reorganized `Ctrl_Exception` throwing in module `Schrödinger`

####  **v0.1.5**

**Done**

* Fixed bug (forgot to define `Ctrl_Exception`)

####  **v0.1.4**

**Done**

* Fixed bug with import-export of `schröd!()`
* Fixed docs of `schröd!()`

####  **v0.1.3**

**Done**

* Fixed bug `Pulse__Δ_BangBang` constructor, added test

####  **v0.1.2**

**Done**

* Improved docs, exceptions
* Fixed bug in `Pulse__Δ_BangBang` constructor

####  **v0.1.1**

**Done**

* Drafted some pulse shapes in module `DOT_RydSim`
* Drafted `Schrödinger` submodule
* Drafted `HW_Descriptions` submodule and added default HW-descr file in `Resources`

**Todo**

* Combination of pulses with time-evolution is untested!!


####  **v0.1.0** Initial version
* Mostly empty
