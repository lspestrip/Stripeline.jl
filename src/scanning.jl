using Quaternions
import Healpix

export timetorotang

function timetorotang(time, rpm)
    if rpm == 0
        0.0
    else
        2 * π * time * (rpm / 60)
    end
end

function genpointings(wheelanglesfn, dir, timerange_s; latitude_deg=0.0)
    
    dirs = Array{Float64}(length(timerange_s), 2)
    ψ = Array{Float64}(length(timerange_s))

    for (idx, time_s) = enumerate(timerange_s)
        (wheel1ang, wheel2ang, wheel3ang) = wheelanglesfn(time_s)
        
        qwheel1 = qrotation(dir, wheel1ang)
        qwheel2 = qrotation([1, 0, 0], wheel2ang)
        qwheel3 = qrotation([0, 0, 1], wheel3ang)
        
        # This is in the ground reference frame
        groundq = qwheel3 * (qwheel2 * qwheel1)
        
        # Now from the ground reference frame to the Earth reference frame
        locq = qrotation([1, 0, 0], deg2rad(90 - latitude_deg))
        earthq = qrotation([0, 0, 1], 2 * π * time_s / 86400)
        quat = earthq * (locq * groundq)
        rotmatr = rotationmatrix(quat)
        
        vector = rotmatr * [0; 0; 1]
        poldir = rotmatr * [1; 0; 0]
        
        # The North for a vector v is just dv/dθ, as θ is the
        # colatitude and moves along the meridian
        (θ, ϕ) = Healpix.vec2ang(vector[1], vector[2], vector[3])
        dirs[idx, :] = [θ, ϕ]
        northdir = [cos(θ) * cos(ϕ); cos(θ) * sin(ϕ); sin(θ)]
        
        cosψ = clamp(dot(northdir, poldir), -1, 1)
        crosspr = northdir × poldir
        sinψ = clamp(dot(crosspr, crosspr), -1, 1)
        ψ[idx] = atan2(cosψ, sinψ) * sign(dot(crosspr, vector))
    end
    
    (dirs, ψ)
end