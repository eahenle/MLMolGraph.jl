function setup_cache(cache)
    rc[:paths][:cache] = cache
    rc[:cache][:primitive] = joinpath(cache, "primitive")
    rc[:cache][:bonded_xtals] = joinpath(cache, "bonded_xtals")
    for dir ∈ rc[:cache]
        if !isdirpath(dir[2]) 
            mkpath(dir[2])
        end
    end
end


function clear_cache()
    try
        rm(rc[:paths][:cache], recursive=true)
    catch
    end
    setup_cache(rc[:paths][:cache])
end


function cached(f::Function, cachefile::String)
    cachefile = joinpath(rc[:paths][:cache], cachefile)
    if isfile(cachefile)
        @debug "Loading cache" cachefile
        @load cachefile obj
    else
        obj = f()
        @save cachefile obj
    end
    return obj
end
