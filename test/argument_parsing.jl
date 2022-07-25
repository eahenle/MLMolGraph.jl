module Arg_test

using Test, MLMolGraph

@testset "target constraints" begin
    args = make_arg_dict()
    @test args[:target] == ""
end

end