__precompile__()
module Stripeline
using PyCall

astropy_time = PyNULL()
astropy_coordinates = PyNULL()
astropy_units = PyNULL()

function __init__()
    copy!(astropy_time, pyimport_conda("astropy.time", "astropy"))
    copy!(astropy_coordinates, pyimport_conda("astropy.coordinates", "astropy"))
    copy!(astropy_units, pyimport_conda("astropy.units", "astropy"))
end


include("instrumentdb.jl")
include("scanning.jl")
include("timesplit.jl")
include("mapmaker.jl")
include("tod_splitter.jl")

end
