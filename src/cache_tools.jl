function setup_cache()
    rc[:cache][:primitive] = joinpath(rc[:paths][:data], "cache", "primitive")
    rc[:cache][:bonded_xtals] = joinpath(rc[:paths][:data], "cache", "bonded_xtals")
    for dir ∈ rc[:cache]
        if !isdirpath(dir[2]) 
            mkpath(dir[2])
        end
    end
end


function clear_cache()
    try
        rm(joinpath(rc[:paths][:data], "cache"), recursive=true)
    catch
    end
    setup_cache()
end


function cached(f::Function, cachefile::String)
    cachefile = joinpath(rc[:paths][:data], "cache", cachefile)
    if isfile(cachefile)
        @load cachefile obj
    else
        obj = f()
        @save cachefile obj
    end
    return obj
end
