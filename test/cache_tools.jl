module Cache_tools_test

using Test, MLMolGraph

@testset "setup/clear cache" begin
    old_cache = rc[:paths][:cache] # save current cache location
    temp_cache = joinpath(rc[:paths][:data], "temp")
    if !isdir(temp_cache)
        mkdir(temp_cache)
    end
    old_pwd = pwd()
    cd(temp_cache)
    setup_cache(pwd())
    # test moving cache
    @test all(isdir.(path[2] for path ∈ rc[:cache]))
    open(joinpath(pwd(), "temp.txt"), "w") do io
        print(io, "foo")
    end
    # test clearing cache
    @test isfile(joinpath(rc[:paths][:cache], "temp.txt"))
    clear_cache()
    @test all(isdir.(path[2] for path ∈ rc[:cache]))
    @test !isfile(joinpath(rc[:paths][:cache], "temp.txt"))
    setup_cache(old_cache) # restore cache location
    cd(old_pwd) # restore working directory path
end

@testset "cached" begin
    clear_cache()
    bar = :baz
    foo = cached("foo.jld2") do 
        return bar
    end
    @test foo == :baz
    @test isfile(joinpath(rc[:paths][:cache], "foo.jld2"))
    bar = :norf
    friz = cached("foo.jld2") do 
        return bar
    end
    @test friz == :baz
    clear_cache()
    friz = cached("foo.jld2") do
        return bar
    end
    @test friz == :norf
end

end
