export TimeChunk, splittime

struct TimeChunk
    starttime::Any
    idx0::Int64
    nsamples::Int64
end

"""
    function splittime(timelen, nchunks, sampfreq, time0)

Given a time span lasting `timelen` seconds, if the sampling
frequency is `sampfreq` Hz, split the span into `nchunks`
sections. Return an array of `TimeChunk` instances. The
parameter `time0` is used as a constant offset for each
member `TimeChunk.starttime`.

# Examples
```julia-repl
julia> splittime(10, 3, 2)
3-element Array{TimeChunk,1}:
 TimeChunk(0.0, 1, 7)     
 TimeChunk(3.33333, 8, 6) 
 TimeChunk(6.66667, 14, 7)
julia> splittime(10, 3, 2, 100) # Same as above, but time0=100
 3-element Array{TimeChunk,1}:
 TimeChunk(100.0, 1, 7)     
 TimeChunk(103.333, 8, 6) 
 TimeChunk(106.667, 14, 7)
```
"""
function splittime(timelen, nchunks, sampfreq, time0 = 0)
    result = Array{TimeChunk}(nchunks)

    nrem = convert(Int, sampfreq * timelen)
    idx0 = 1
    for i = 1:nchunks
        starttime = time0 + (i - 1) * timelen / nchunks

        if i < nchunks
            n = convert(Int, round(nrem // (nchunks - i + 1)))
        else
            n = nrem
        end

        result[i] = TimeChunk(starttime, idx0, n)

        idx0 += n
        nrem -= n
    end

    result
end
