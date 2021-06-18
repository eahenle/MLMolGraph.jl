module Cache_tools_test

using Test, MLMolGraph

@testset "setup/clear cache" begin
    old_cache = rc[:paths][:cache] # save current cache location
    mkdir("temp")
    cd("temp")
    setup_cache(pwd())
    # test moving cache
    @test isdirpath.(rc[:cache][2])
    open(joinpath(pwd(), "temp.txt"), "w") do io
        print(io, "foo")
    end
    # test clearing cache
    @test isfile(joinpath(rc[:paths][:cache], "temp.txt"))
    clear_cache()
    @test isdirpath.(rc[:cache][2])
    @test !isfile(joinpath(rc[:paths][:cache], "temp.txt"))
    setup_cache(old_cache) # restore cache location
end

@testset "cached" begin
    clear_cache()
    bar = :baz
    foo = cached("foo.jld2") do 
        return bar
    end
    @test foo == :bar
    @test isfile(joinpath(rc[:data][:cache], "foo.jld2"))
    bar = :norf
    friz = cached("foo.jld2") do 
        return bar
    end
    @test friz == :bar
    clear_cache()
    friz = cached("foo.jld2") do
        return bar
    end
    @test friz == :norf
end

end
