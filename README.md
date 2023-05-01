# DOT_RydSim.jl
Quantum simulation of arrays of Rydberg atoms &amp; related stuff

## Version History

Time goes up!

####  **v0.1.12**

**Done**

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
