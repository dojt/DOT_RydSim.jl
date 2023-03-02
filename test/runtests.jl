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

using Test

#
using DOT_RydSim



function test__dummy(option_set::Symbol...)                                                  # dummy()
    ALL_OPTS = [ ]
    option_set ⊆ ALL_OPTS  || throw(ArgumentError("Options not recognized: $(setdiff(option_set,ALL_OPTS))"))

	@testset verbose=true "......." begin
        @testset "......." begin
        end
    end
end #^ test__dummy()




@testset verbose=true "Testing DOT_RydSim.jl" begin
    test__dummy()
end

#runtests.jl
#EOF
