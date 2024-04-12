export TodNoiseProperties, build_noise_properties, tod2map_mpi, baseline2map_mpi, destripe, baselines_covmat, condnumber_mpi, get_QU_weights

using LinearAlgebra

try
    import MPI
catch
end

@doc raw"""
This structure holds a number of parameters relative to the noise level and
the 1/f baselines measured or simulated for a certain polarimeter.

Field              | Type           | Meaning
:----------------- |:-------------- |:----------------------------------------------
`pol_number`       | Int            | ID number of the polarimeter
`sigma`            | Float          | sigma of white noise
`num_of_samples`   | Int            | total number of samples
`baseline_lengths` | Array{Int}     | Array containing the length of each baseline

The structure provides a constructor that accepts the following
three keyword parameters:

- `pol`: number of the polarimeter
- `rms`: the sigma level
- `baseline_lengths`: array containing the elements of each baseline

# Example

The followng code creates a `TodNoiseProperties` object referring to polarimeter #1
(int this case, STRIP31) and specifying 0.002463 as the noise level. The data sample
contains 500 samples, and it is split into 5 baselines of equal length.

```julia
prop = TodNoiseProperties(
    pol = 1,
    rms = 0.002463, 
    baselines = Int[100, 100, 100, 100, 100, 100],
)
```

"""
struct TodNoiseProperties
    polarimeter::Int
    sigma::Array{Float64}
    number_of_samples::Int
    baseline_lengths::Array{Int}

    TodNoiseProperties(;
        pol::Int,
        rms::Array{Float64},
        baselines::Array{Int}) = length(rms) != sum(baselines) ? error(" length(rms) /= sum(baselines) ") : new(pol, rms, sum(baselines), baselines)
end


@doc raw"""
    function build_noise_properties(detector_list, rms_list, num_of_baselines, num_of_samples) -> data_properties

This function builds a list of `TodNoiseProperties` object. These are
needed for destriping and calculation of `baselines_covmat`.

It requires in input 4 arrays containing:
- the ID number of the polarimeters that the current rank will simulate
- the white noise sigma for each polarimeter that the current rank will simulate
- the number of 1/f baselines for each polarimeter that the current rank will simulate
- the total number of samples for each polarimeter that the current rank will simulate

N.B. the ID numbers, the number of 1/f baselines and the total number of samples can be obtained using the function `get_chunk_properties`.

It returns an array of `TodNoiseProperties`, of length equal to the number of polarimeters simulated by current rank.

"""
function build_noise_properties(detector_list, rms_list, num_of_baselines, num_of_samples)
    @assert length(detector_list) == length(rms_list)
    @assert length(detector_list) == length(num_of_baselines)
    @assert length(detector_list) == length(num_of_samples)

    data_properties  = Array{TodNoiseProperties}(undef, length(detector_list))
    for i in 1:length(detector_list)
        # We use baselines of equal lengths
        baseline_lengths  = repeat([Int(num_of_samples[i] / num_of_baselines[i])], num_of_baselines[i])
        data_properties[i] = TodNoiseProperties(
            pol = detector_list[i],
            rms = rms_list[detector_list[i]],
            baselines = baseline_lengths,
        )
    end
    data_properties
end


@doc raw"""
    get_QU_weights(twopsi,mode="sum") -> wq ,wu

Given twice the detector polarization angle, the function returns the weights mapping the sky polarization signal (Q_sky, U_sky) into the detector output

y = wq Q_sky + wu Q_sky +noise

Mode is one of

- "Q" : demodulated Q output
- "U" : demodulated U output
- "sum": demodulated Q + U 
- "diff": demodulated Q - U 

"""
function get_QU_weights(twopsi;tod_mode="sum")
    lmode = lowercase(tod_mode)
    
    if (lmode == lowercase("Q"))
        # y = Q cos(2psi) + U sin(2psi) +noise
        wq = cos(twopsi)
        wu = sin(twopsi)
    elseif (lmode == lowercase("U"))
        # y = - Q sin(2psi) + U cos(2psi) +noise
        wq = -sin(twopsi)
        wu =  cos(twopsi)
    elseif (lmode == lowercase("sum"))
        # y = Q (cos(2psi) - sin(2psi)) + U (cos(2psi) + sin(2psi)) +noise
        wq = sqrt(2.)*cos(twopsi +π/4)
        wu = sqrt(2.)*cos(twopsi -π/4)
    elseif (lmode == lowercase("diff"))
        # y = Q (cos(2psi) + sin(2psi)) - U (cos(2psi) - sin(2psi)) +noise
        wq =  sqrt(2.)*cos(twopsi -π/4)
        wu = -sqrt(2.)*cos(twopsi +π/4)
    else
        wq = NaN
        wu = NaN
    end

    wq ,wu
end

function get_QU_weights(twopsi,σ_k;tod_mode="sum")
    
    wq ,wu = get_QU_weights(twopsi,tod_mode=tod_mode)

    wq /= σ_k
    wu /= σ_k

    wq ,wu
end


@doc raw"""
    tod2map(pix_idx, tod, num_of_pixels, twopsi; comm = nothing) -> binned_map
    tod2map(pix_idx, tod, num_of_pixels, twopsi, data_properties; comm = nothing) -> binned_map

This function creates a binned map from the time-ordered data kept in
the array `tod`, assuming that each sample is observing the pixel
whose index is in `pix_idx` at an angle half of `twopsi`. The parameter `num_of_pixels` contains
the number of pixels in the Healpix map, and it is used as an upper
bound for the values in `pix_idx`. The parameter `comm` must be a MPI
communicator, or `nothing` if you are not using MPI.

If `comm` is not `nothing`, the function is parallelized using MPI, and
each process computes a map from its available data.  All partial maps
are then combined together with `MPI.Allreduce`. The function returns
an array containing the binned map.

If Array of structures `TodNoiseProperties`is passed to the function, the
output binned map will be a weighted binned map. Each sample is weighted
according to the inverse white noise variance sigma^2 of the corrispondingù
polarimeter. In this way, the less noisy polarimeters will count more in the
estimation of the map.

# Requirements
- The length of the arrays `twopsi`, `pix_idx`, and `tod` must be the same

"""
function tod2map_mpi(pix_idx, tod, num_of_pixels, twopsi;
                       comm = nothing, unseen = NaN, tod_mode = "sum" )
    @assert length(pix_idx) == length(tod)
    @assert length(twopsi)  == length(tod)

    #lpl
    #This solves for x = (R^t N^-1 R)^-1 (R^t N^-1 y) for demodulated polarization tod
    # x 2n_pix array, ordering x : x^t = (mQ_1,mU_1,mQ_2,mU_2,...,mQ_npix,mU_npix)
    # y n_sample tod
    # N n_sample white noise covariance in time domain, here N = 1 
    # R (n_sample x 2n_pix) pointing matrix so that
    # R_i,2j-1 = Wq_i  if sample i falls into pix j, 0 otherwise
    # R_i,2j   = Wu_i  if sample i falls into pix j, 0 otherwise
    # Wq, and Wu depend on the specific tod considered, e.g. for demodulated Q tod
    # Wq_i = cos (2 psi_i)
    # Wu_i = sin (2 psi_i)    
    # (R^tN^-1R) is a (2x2)-block diagonal matrix  
        
    T = eltype(tod)
    N = eltype(pix_idx)

    #(R^t N^1 R) on-diagonal blocks are symmetric, so only 3 numbers/block
    #are needed
    #Allocating a num_of_pixels array rather than just observed pixels is
    #wasteful, but likely small compared to the tod, and simplyfies life 

    binned_map  = zeros(T, (2,num_of_pixels))
    binned_rnr  = zeros(T,(3,num_of_pixels))

    #don't need these anymore, since we check on (R^t N^-1 R) 
    #partial_hits = zeros(N, num_of_pixels)
    #hits = zeros(N, num_of_pixels)

    @inbounds for i in eachindex(pix_idx)
        wq , wu = get_QU_weights(twopsi[i],tod_mode=tod_mode) 

        binned_map[1,pix_idx[i]] += tod[i]*wq
        binned_map[2,pix_idx[i]] += tod[i]*wu

        binned_rnr[1,pix_idx[i]] += wq*wq
        binned_rnr[2,pix_idx[i]] += wq*wu
        binned_rnr[3,pix_idx[i]] += wu*wu

    end

    if comm != nothing
        MPI.Allreduce!(MPI.IN_PLACE, binned_map, MPI.SUM, comm)
        MPI.Allreduce!(MPI.IN_PLACE, binned_rnr, MPI.SUM, comm)
    end

    @inbounds for i in eachindex(binned_rnr[1,:])
        if (binned_rnr[1,i] > 0.)
            #M =  |a b|    M^-1  = | c -b| / det(M) 
            #     |b c|            |-b  a| 
            det = binned_rnr[1,i]*binned_rnr[3,i] -binned_rnr[2,i]^2
            if(det >0.) 
                mq  = (+binned_rnr[3,i]*binned_map[1,i]
                       -binned_rnr[2,i]*binned_map[2,i])/det
                mu  = (-binned_rnr[2,i]*binned_map[1,i]
                       +binned_rnr[1,i]*binned_map[2,i])/det
                
                binned_map[1,i] = mq
                binned_map[2,i] = mu
            else
                binned_map[:,i] .= unseen
            end
        else
            binned_map[:,i] .= unseen
        end
    end
    binned_rnr = nothing
    binned_map
end

function tod2map_mpi(pix_idx, tod, num_of_pixels, twopsi,
                       data_properties; comm = nothing, unseen = NaN,
                       tod_mode = "sum")

    @assert length(pix_idx) == length(tod)
    @assert length(twopsi)  == length(tod)

    T = eltype(tod)
    
    binned_map   = zeros(T, (2,num_of_pixels))
    binned_rnr  = zeros(T,(3,num_of_pixels))

    offset = 0

    #loop on detectors
    for j in eachindex(data_properties)                 
        #end_idx = start_idx + data_properties[j].number_of_samples - 1

        #loop on samples
        #for i in start_idx:end_idx             
        for i in 1:data_properties[j].number_of_samples             
            ioff = i+offset

            wq , wu = get_QU_weights(twopsi[ioff],data_properties[j].sigma[i],tod_mode=tod_mode)
            #sigma2 = (data_properties[j].sigma)^2       
            
            tos = tod[ioff]/data_properties[j].sigma[i]

            #partial_map[1,pix_idx[i]] += wq*tod[i]/sigma2
            #partial_map[2,pix_idx[i]] += wu*tod[i]/sigma2

            binned_map[1,pix_idx[ioff]] += wq*tos#*tod[i]/sigma2
            binned_map[2,pix_idx[ioff]] += wu*tos#*tod[i]/sigma2

            binned_rnr[1,pix_idx[ioff]] += wq*wq#/sigma2
            binned_rnr[2,pix_idx[ioff]] += wq*wu#/sigma2
            binned_rnr[3,pix_idx[ioff]] += wu*wu#/sigma2

        end
        #start_idx += data_properties[j].number_of_samples
        offset += data_properties[j].number_of_samples
    end

    if comm != nothing
        MPI.Allreduce!(MPI.IN_PLACE, binned_map, MPI.SUM, comm)
        MPI.Allreduce!(MPI.IN_PLACE, binned_rnr, MPI.SUM, comm)
    end

    @inbounds for i in eachindex(binned_rnr[1,:])
        if (binned_rnr[1,i] > 0.)
            det = binned_rnr[1,i]*binned_rnr[3,i] -binned_rnr[2,i]^2
            if (det > 0.) 
                mq  = (+binned_rnr[3,i]*binned_map[1,i]
                       -binned_rnr[2,i]*binned_map[2,i])/det
                mu  = (-binned_rnr[2,i]*binned_map[1,i]
                       +binned_rnr[1,i]*binned_map[2,i])/det

                binned_map[1,i] = mq
                binned_map[2,i] = mu
            else
                binned_map[:,i] .= unseen
            end                
        else
            binned_map[:,i] .= unseen
        end
    end
    binned_rnr = nothing
    binned_map
end

#preallocate map and rnr buffers to see if this fixes MPI memory leaks 
function tod2map_mpi!(pix_idx, tod, num_of_pixels, twopsi,
    data_properties, binned_map, binned_rnr; comm = nothing, unseen = NaN,
    tod_mode = "sum")

    @assert length(pix_idx) == length(tod)
    @assert length(twopsi)  == length(tod)

    T = eltype(tod)

    binned_map .= zero(T)
    #binned_rnr .= zero(T)

    offset = 0

    #loop on detectors
    for j in eachindex(data_properties)                 
        #end_idx = start_idx + data_properties[j].number_of_samples - 1

        #loop on samples
        #for i in start_idx:end_idx             
        for i in 1:data_properties[j].number_of_samples             
            ioff = i+offset

            wq , wu = get_QU_weights(twopsi[ioff],data_properties[j].sigma[i],tod_mode=tod_mode)
            #sigma2 = (data_properties[j].sigma)^2       

            tos = tod[ioff]/data_properties[j].sigma[i]

            #partial_map[1,pix_idx[i]] += wq*tod[i]/sigma2
            #partial_map[2,pix_idx[i]] += wu*tod[i]/sigma2

            binned_map[1,pix_idx[ioff]] += wq*tos#*tod[i]/sigma2
            binned_map[2,pix_idx[ioff]] += wu*tos#*tod[i]/sigma2

            #binned_rnr[1,pix_idx[ioff]] += wq*wq#/sigma2
            #binned_rnr[2,pix_idx[ioff]] += wq*wu#/sigma2
            #binned_rnr[3,pix_idx[ioff]] += wu*wu#/sigma2

        end
        #start_idx += data_properties[j].number_of_samples
        offset += data_properties[j].number_of_samples
    end

    if comm != nothing
        MPI.Allreduce!(MPI.IN_PLACE, binned_map, MPI.SUM, comm)
        #MPI.Allreduce!(MPI.IN_PLACE, binned_rnr, MPI.SUM, comm)
    end

    @inbounds for i in eachindex(binned_rnr[1,:])
        if (binned_rnr[1,i] > 0.)
            det = binned_rnr[1,i]*binned_rnr[3,i] -binned_rnr[2,i]^2
            if (det > 0.) 
                mq  = (+binned_rnr[3,i]*binned_map[1,i]
                    -binned_rnr[2,i]*binned_map[2,i])/det
                mu  = (-binned_rnr[2,i]*binned_map[1,i]
                    +binned_rnr[1,i]*binned_map[2,i])/det

                binned_map[1,i] = mq
                binned_map[2,i] = mu
            else
                binned_map[:,i] .= unseen
            end                
        else
            binned_map[:,i] .= unseen
        end
    end
end

function condnumber_mpi(pix_idx, num_of_pixels, twopsi,
                          data_properties; comm = nothing, unseen = NaN,
                          tod_mode = "sum" )
    @assert length(pix_idx) == length(twopsi)

    #lpl
    #computes several diagnostics of mapmaking stability:
    #* caustics maps :
    #  (R^t N^-1 R)^-1 R^t N^-1 E
    #  E^t = (1, 1, ..., 1) 
    #* determinant
    #* condition number
    
    T = eltype(twopsi)
    
    binned_map  = zeros(T,(4,num_of_pixels))
    binned_rnr  = zeros(T,(3,num_of_pixels))

    offset = 0

    #loop on detectors
    for j in eachindex(data_properties)                 
        #end_idx = start_idx + data_properties[j].number_of_samples - 1

        #loop on samples
        for i in 1:data_properties[j].number_of_samples
            ioff = i +offset             

            wq , wu = get_QU_weights(twopsi[ioff],tod_mode=tod_mode)
            sigma2 = (data_properties[j].sigma[i])^2       

            binned_map[1,pix_idx[ioff]] += wq*(wq +wu)/sigma2
            binned_map[2,pix_idx[ioff]] += wu*(wq +wu)/sigma2

            binned_rnr[1,pix_idx[ioff]] += wq*wq/sigma2
            binned_rnr[2,pix_idx[ioff]] += wq*wu/sigma2
            binned_rnr[3,pix_idx[ioff]] += wu*wu/sigma2

        end
        #start_idx += data_properties[j].number_of_samples
        offset += data_properties[j].number_of_samples

    end

    if comm != nothing
        MPI.Allreduce!(MPI.IN_PLACE,binned_map, MPI.SUM, comm)
        MPI.Allreduce!(MPI.IN_PLACE,binned_rnr, MPI.SUM, comm)
    end

    @inbounds for i in eachindex(binned_rnr[1,:])
        if(binned_rnr[1,i] > 0) 
            det = binned_rnr[1,i]*binned_rnr[3,i] -binned_rnr[2,i]^2
            if (det > 0) 
                mq  = (+binned_rnr[3,i]*binned_map[1,i]
                       -binned_rnr[2,i]*binned_map[2,i])/det
                mu  = (-binned_rnr[2,i]*binned_map[1,i]
                       +binned_rnr[1,i]*binned_map[2,i])/det

                m = (binned_rnr[1,i] +binned_rnr[3,i])/2
                λ_1 = m + sqrt(m^2 - det)
                λ_2 = m - sqrt(m^2 - det)
                
                binned_map[1,i] = mq
                binned_map[2,i] = mu
                binned_map[3,i] = det
                binned_map[4,i] = λ_1/λ_2
            else
                binned_map[:,i] .= unseen
            end
        else
            binned_map[:,i] .= unseen
        end
    end

    binned_rnr = nothing
    binned_map
end


function binned_noise_variance_mpi(pix_idx, num_of_pixels, twopsi, data_properties; comm = nothing, unseen = NaN,
    tod_mode = "sum" )
    @assert length(pix_idx) == length(twopsi)

    #lpl
    #computes the pixel (inverse) noise variance of binned map (R N^-1 R)

    T = eltype(twopsi)

    binned_rnr  = zeros(T,(3,num_of_pixels))
    offset = 0
    #loop on detectors
    for j in eachindex(data_properties)                 
        #end_idx = start_idx + data_properties[j].number_of_samples - 1

        #loop on samples
        for i in 1:data_properties[j].number_of_samples
            ioff = i +offset             

            wq , wu = get_QU_weights(twopsi[ioff],tod_mode=tod_mode)
            sigma2 = (data_properties[j].sigma[i])^2       

            binned_rnr[1,pix_idx[ioff]] += wq*wq/sigma2
            binned_rnr[2,pix_idx[ioff]] += wq*wu/sigma2
            binned_rnr[3,pix_idx[ioff]] += wu*wu/sigma2

        end
        #start_idx += data_properties[j].number_of_samples
        offset += data_properties[j].number_of_samples

    end

    if comm != nothing
        MPI.Allreduce!(MPI.IN_PLACE,binned_rnr, MPI.SUM, comm)
    end

    @inbounds for i in eachindex(binned_rnr[1,:])
        if(binned_rnr[1,i] > 0) 
            det = binned_rnr[1,i]*binned_rnr[3,i] -binned_rnr[2,i]^2
            if (det > 0) 
            else
                #binned_map[:,i] .= unseen
            end
        else
            binned_rnr[:,i] .= unseen
        end
    end

    binned_rnr 
end


@doc raw"""
    baseline2map_mpi(pix_idx, baselines, baseline_lengths, num_of_pixels; comm = nothing)-> noise_map
    baseline2map_mpi(pix_idx, baselines, num_of_pixels, data_properties; comm = nothing) -> noise_map

This is a MPI based function: each MPI process computes a map from its
available data.  All partial maps are then combined together with
MPI.Allreduce.
The function returns an array containing the binned map.

If Array of structures `TodNoiseProperties`is passed to the function (instead of `baseline_lengths`) the output binned map will be a weighted binned map.
Each sample will be weighted according to the inverse white noise variance sigma^2 of the corrisponding polarimeter.
In this way, the less noisy polarimeters will count more in the estimation of the map.

# Requirements
- The length of `baselines` and `baseline_lengths` must be the same;
- The value `sum(baseline_lengths)` must be the same as the length of `pix_idx`.
"""
baseline2map_mpi 

function baseline2map_mpi(pix_idx, baselines, baseline_lengths, num_of_pixels,
                            twopsi;
                            comm = nothing, unseen = NaN, tod_mode = "sum")

    #lpl                                                                       
    #This solves for x = (R^t N^-1 R)^-1 (R^t N^-1 Fa) 
    #It's basically the same function as tod2map but reimplemented directly
    #to avoid allocating an entire new tod, I guess?
    
    T = eltype(baselines)
    N = eltype(pix_idx)

    binned_map   = zeros(T, (2,num_of_pixels))
    binned_rnr   = zeros(T, (3,num_of_pixels))

#    partial_hits = zeros(N, num_of_pixels)
#    hits = zeros(N, num_of_pixels)

    startidx = 1

    for i in eachindex(baseline_lengths)
        endidx = baseline_lengths[i] + startidx - 1

        for j in startidx:endidx

            wq , wu = get_QU_weights(twopsi[j],tod_mode=tod_mode)

            binned_map[1,pix_idx[j]] += baselines[i]*wq
            binned_map[2,pix_idx[j]] += baselines[i]*wu

            binned_rnr[1,pix_idx[j]] += wq*wq
            binned_rnr[2,pix_idx[j]] += wq*wu
            binned_rnr[3,pix_idx[j]] += wu*wu

 #           partial_hits[pix_idx[j]] += 1
        end
        startidx += baseline_lengths[i]
    end

    if comm != nothing
        MPI.Allreduce!(MPI.IN_PLACE,binned_map, MPI.SUM, comm)
        MPI.Allreduce!(MPI.IN_PLACE,binned_rnr, MPI.SUM, comm)
    end
 
    @inbounds for i in eachindex(binned_rnr[1,:])
        if (binned_rnr[1,i] > 0.)
            det = binned_rnr[1,i]*binned_rnr[3,i] -binned_rnr[2,i]^2
            if (det >0)
                mq  = (+binned_rnr[3,i]*binned_map[1,i]
                       -binned_rnr[2,i]*binned_map[2,i])/det
                mu  = (-binned_rnr[2,i]*binned_map[1,i]
                       +binned_rnr[1,i]*binned_map[2,i])/det
                
                binned_map[1,i] = mq
                binned_map[2,i] = mu
            else
                binned_map[:,i] .= unseen
            end
        else
            binned_map[:,i] .= unseen
        end
    end

    binned_rnr = nothing
    binned_map
end

function baseline2map_mpi(pix_idx, baselines, num_of_pixels, twopsi,
                            data_properties, num_of_baselines;
                            comm = nothing, unseen = NaN, tod_mode = "sum")

    T = eltype(baselines)

    binned_map = zeros(T, (2,num_of_pixels))
    binned_rnr  = zeros(T, (3,num_of_pixels))
    
    #startidx = 1
    offset = 0
    baseline_idx = 1

    for det_idx in eachindex(data_properties)  #loop on detectors
        sigma_offset = 0
        for sample_idx in eachindex(data_properties[det_idx].baseline_lengths)
            #endidx = data_properties[det_idx].baseline_lengths[sample_idx] + startidx - 1

            @inbounds for j in 1:data_properties[det_idx].baseline_lengths[sample_idx]
                joff   = j+offset
                jsoff  = j+sigma_offset
                #partial_map[pix_idx[j]] += baselines[baseline_idx] / (data_properties[det_idx].sigma)^2
                #sigma2 = data_properties[det_idx].sigma^2

                wq , wu = get_QU_weights(twopsi[joff],data_properties[det_idx].sigma[jsoff],tod_mode=tod_mode)
                
                bos = baselines[baseline_idx]/data_properties[det_idx].sigma[jsoff]
                #partial_map[1,pix_idx[j]] +=
                #    wq*(baselines[baseline_idx]/sigma2)
                #partial_map[2,pix_idx[j]] +=
                #    wu*(baselines[baseline_idx]/sigma2)

                binned_map[1,pix_idx[joff]] += wq*bos
                binned_map[2,pix_idx[joff]] += wu*bos

                binned_rnr[1,pix_idx[joff]] += wq*wq#/sigma2
                binned_rnr[2,pix_idx[joff]] += wq*wu#/sigma2
                binned_rnr[3,pix_idx[joff]] += wu*wu#/sigma2
                
            end
            #startidx += data_properties[det_idx].baseline_lengths[sample_idx]
            offset       += data_properties[det_idx].baseline_lengths[sample_idx]
            sigma_offset += data_properties[det_idx].baseline_lengths[sample_idx]
            baseline_idx += 1
        end
    end

    if comm != nothing
        MPI.Allreduce!(MPI.IN_PLACE,binned_map, MPI.SUM, comm)
        MPI.Allreduce!(MPI.IN_PLACE,binned_rnr, MPI.SUM, comm)
    end

    @inbounds for i in eachindex(binned_rnr[1,:])
        if (binned_rnr[1,i] > 0)

            det = binned_rnr[1,i]*binned_rnr[3,i] -binned_rnr[2,i]^2
            if (det > 0) 
                mq  = (+binned_rnr[3,i]*binned_map[1,i]
                       -binned_rnr[2,i]*binned_map[2,i])/det
                mu  = (-binned_rnr[2,i]*binned_map[1,i]
                       +binned_rnr[1,i]*binned_map[2,i])/det

                binned_map[1,i] = mq
                binned_map[2,i] = mu
            else
                binned_map[:,i] .= unseen
            end
        else
            binned_map[:,i] .= unseen
        end
    end

    binned_rnr = nothing
    binned_map
end

function baseline2map_mpi!(pix_idx, baselines, num_of_pixels, twopsi,
        data_properties, num_of_baselines, binned_map, binned_rnr;
        comm = nothing, unseen = NaN, tod_mode = "sum")

    T = eltype(baselines)

    #startidx = 1
    offset = 0
    baseline_idx = 1

    binned_map .= zero(T)
    #binned_rnr .= zero(T)

    for det_idx in eachindex(data_properties)  #loop on detectors
        sigma_offset = 0
        for sample_idx in eachindex(data_properties[det_idx].baseline_lengths)
        #endidx = data_properties[det_idx].baseline_lengths[sample_idx] + startidx - 1

        @inbounds for j in 1:data_properties[det_idx].baseline_lengths[sample_idx]
            joff   = j+offset
            jsoff  = j+sigma_offset
            #partial_map[pix_idx[j]] += baselines[baseline_idx] / (data_properties[det_idx].sigma)^2
            #sigma2 = data_properties[det_idx].sigma^2

            wq , wu = get_QU_weights(twopsi[joff],data_properties[det_idx].sigma[jsoff],tod_mode=tod_mode)

            bos = baselines[baseline_idx]/data_properties[det_idx].sigma[jsoff]
            #partial_map[1,pix_idx[j]] +=
            #    wq*(baselines[baseline_idx]/sigma2)
            #partial_map[2,pix_idx[j]] +=
            #    wu*(baselines[baseline_idx]/sigma2)

            binned_map[1,pix_idx[joff]] += wq*bos
            binned_map[2,pix_idx[joff]] += wu*bos

            #binned_rnr[1,pix_idx[joff]] += wq*wq#/sigma2
            #binned_rnr[2,pix_idx[joff]] += wq*wu#/sigma2
            #binned_rnr[3,pix_idx[joff]] += wu*wu#/sigma2

        end
        #startidx += data_properties[det_idx].baseline_lengths[sample_idx]
        offset       += data_properties[det_idx].baseline_lengths[sample_idx]
        sigma_offset += data_properties[det_idx].baseline_lengths[sample_idx]
        baseline_idx += 1
        end
    end

    if comm != nothing
        MPI.Allreduce!(MPI.IN_PLACE,binned_map, MPI.SUM, comm)
        #MPI.Allreduce!(MPI.IN_PLACE,binned_rnr, MPI.SUM, comm)
    end

    @inbounds for i in eachindex(binned_rnr[1,:])
        if (binned_rnr[1,i] > 0)

            det = binned_rnr[1,i]*binned_rnr[3,i] -binned_rnr[2,i]^2
            if (det > 0) 
                mq  = (+binned_rnr[3,i]*binned_map[1,i]
                    -binned_rnr[2,i]*binned_map[2,i])/det
                mu  = (-binned_rnr[2,i]*binned_map[1,i]
                    +binned_rnr[1,i]*binned_map[2,i])/det

                binned_map[1,i] = mq
                binned_map[2,i] = mu
            else
                binned_map[:,i] .= unseen
            end
        else
            binned_map[:,i] .= unseen
        end
    end

end

function applyz_and_sum(pix_idx, tod, num_of_pixels, twopsi, data_properties,
                          num_of_baselines;
                          comm = nothing, unseen = NaN, tod_mode = "sum")
    
    @assert length(tod) == length(pix_idx)

    #lpl
    #maps a tod -> baselines, according to
    #F^t N^-1 [ 1 - R(R^t N^-1 R) R^t N^-1 ] y = F^t N^-1 [y - R tod2map(y)]
    #so we can mostly use what's in place just accounting for the fact that
    #tod2map gives a (2 x n_pix) array  describing (Q,U) maps, and R encodes 
    #for polarization reference frame rotation
    
    baselines_sum = zeros(eltype(tod), num_of_baselines)

    binned_map = tod2map_mpi(pix_idx,
        tod,
        num_of_pixels,
        twopsi,                       
        data_properties,
        comm = comm,
        unseen = unseen,
        tod_mode = tod_mode)

    baseline_idx = 1
    #startidx = 1

    offset       = 0  
    for det_idx in eachindex(data_properties)  #loop on detectors
        sigma_offset = 0
        for baseline_idx in eachindex(data_properties[det_idx].baseline_lengths)
            #endidx = data_properties[det_idx].baseline_lengths[baseline_idx] + startidx - 1

            for j in 1:data_properties[det_idx].baseline_lengths[baseline_idx]
                joff  = j +offset
                jsoff = j +sigma_offset 
                wq , wu = get_QU_weights(twopsi[joff],tod_mode=tod_mode)
                baselines_sum[baseline_idx] +=
                    (tod[joff]
                     -wq*binned_map[1,pix_idx[joff]]
                     -wu*binned_map[2,pix_idx[joff]]
                    ) / (data_properties[det_idx].sigma[jsoff])^2
            end
            sigma_offset += data_properties[det_idx].baseline_lengths[baseline_idx]
            offset       += data_properties[det_idx].baseline_lengths[baseline_idx]
            #startidx += data_properties[det_idx].baseline_lengths[baseline_idx]
            baseline_idx += 1
        end
    end

    baselines_sum
end

function applyz_and_sum!(pix_idx, tod, num_of_pixels, twopsi, data_properties,
        num_of_baselines, baselines_sum, binned_map, binned_rnr;
        comm = nothing, unseen = NaN, tod_mode = "sum")

    @assert length(tod) == length(pix_idx)

    #lpl
    #maps a tod -> baselines, according to
    #F^t N^-1 [ 1 - R(R^t N^-1 R) R^t N^-1 ] y = F^t N^-1 [y - R tod2map(y)]
    #so we can mostly use what's in place just accounting for the fact that
    #tod2map gives a (2 x n_pix) array  describing (Q,U) maps, and R encodes 
    #for polarization reference frame rotation

    #baselines_sum = zeros(eltype(tod), num_of_baselines)
    baselines_sum .= zero(eltype(tod))

    tod2map_mpi!(pix_idx,
        tod,
        num_of_pixels,
        twopsi,                       
        data_properties,
        binned_map,
        binned_rnr,
        comm = comm,
        unseen = unseen,
        tod_mode = tod_mode)

    baseline_idx = 1
    #startidx = 1

    offset       = 0  
    for det_idx in eachindex(data_properties)  #loop on detectors
        sigma_offset = 0
        for baseline_idx in eachindex(data_properties[det_idx].baseline_lengths)
            #endidx = data_properties[det_idx].baseline_lengths[baseline_idx] + startidx - 1

            for j in 1:data_properties[det_idx].baseline_lengths[baseline_idx]
                joff  = j +offset
                jsoff = j +sigma_offset 
                wq , wu = get_QU_weights(twopsi[joff],tod_mode=tod_mode)
                baselines_sum[baseline_idx] +=
                    (tod[joff]
                    -wq*binned_map[1,pix_idx[joff]]
                    -wu*binned_map[2,pix_idx[joff]]
                    ) / (data_properties[det_idx].sigma[jsoff])^2
            end
            sigma_offset += data_properties[det_idx].baseline_lengths[baseline_idx]
            offset       += data_properties[det_idx].baseline_lengths[baseline_idx]
            #startidx += data_properties[det_idx].baseline_lengths[baseline_idx]
            baseline_idx += 1
        end
    end

end

function applya(baselines, pix_idx, num_of_baselines, num_of_pixels, twopsi,
                  data_properties;
                  comm = nothing, unseen = NaN, tod_mode = "sum")

    #lpl
    #maps baselines -> baselines, according to 
    #F^t N^-1 [1 -R(R^t N^-1 R)^-1 R^t N^-1] Fa

    @assert length(baselines) == num_of_baselines

    baselines_sum = zeros(eltype(baselines), num_of_baselines)
    total_sum = zero(eltype(baselines))

    binned_map = baseline2map_mpi(pix_idx,
        baselines,
        num_of_pixels,
        twopsi,                                                 
        data_properties,
        num_of_baselines,
        comm = comm,
        unseen = unseen,
        tod_mode = tod_mode)

    #startidx = 1
    offset = 0
    baseline_idx = 1

    for det_idx in eachindex(data_properties)
        sigma_offset = 0
        for baseline_idx in eachindex(data_properties[det_idx].baseline_lengths)
            #endidx = data_properties[det_idx].baseline_lengths[baseline_idx] + startidx - 1

            for j in 1:data_properties[det_idx].baseline_lengths[baseline_idx]
                joff  = j +offset
                jsoff = j +sigma_offset
                wq , wu = get_QU_weights(twopsi[joff],tod_mode=tod_mode)
                baselines_sum[baseline_idx] +=
                    (baselines[baseline_idx]
                     -wq*binned_map[1,pix_idx[joff]]
                     -wu*binned_map[2,pix_idx[joff]]
                     ) / (data_properties[det_idx].sigma[jsoff])^2
            end

            #startidx += data_properties[det_idx].baseline_lengths[baseline_idx]
            offset       += data_properties[det_idx].baseline_lengths[baseline_idx]
            sigma_offset += data_properties[det_idx].baseline_lengths[baseline_idx]

            baseline_idx += 1 #lpl: why ?
        end
    end

    #needed to assure that sum(baselines)==0

    if comm != nothing
        total_sum = MPI.Allreduce([sum(baselines)], MPI.SUM, comm)[1]
    else
        total_sum = sum(baselines)
    end

    baselines_sum .+= total_sum
end

function applya!(baselines, pix_idx, num_of_baselines, num_of_pixels, twopsi,
        data_properties, baselines_sum, binned_map, binned_rnr;
        comm = nothing, unseen = NaN, tod_mode = "sum")

        #lpl
        #maps baselines -> baselines, according to 
        #F^t N^-1 [1 -R(R^t N^-1 R)^-1 R^t N^-1] Fa

        @assert length(baselines) == num_of_baselines

        #baselines_sum = zeros(eltype(baselines), num_of_baselines)
        baselines_sum .= zero(eltype(baselines))
        total_sum = zero(eltype(baselines))

        baseline2map_mpi!(pix_idx,
            baselines,
            num_of_pixels,
            twopsi,                                                 
            data_properties,
            num_of_baselines,
            binned_map,
            binned_rnr,
            comm = comm,
            unseen = unseen,
            tod_mode = tod_mode)

        #startidx = 1
        offset = 0
        baseline_idx = 1

        for det_idx in eachindex(data_properties)
            sigma_offset = 0
            for baseline_idx in eachindex(data_properties[det_idx].baseline_lengths)
                #endidx = data_properties[det_idx].baseline_lengths[baseline_idx] + startidx - 1

                for j in 1:data_properties[det_idx].baseline_lengths[baseline_idx]
                    joff  = j +offset
                    jsoff = j +sigma_offset
                    wq , wu = get_QU_weights(twopsi[joff],tod_mode=tod_mode)
                    baselines_sum[baseline_idx] +=
                        (baselines[baseline_idx]
                           -wq*binned_map[1,pix_idx[joff]]
                           -wu*binned_map[2,pix_idx[joff]]
                           ) / (data_properties[det_idx].sigma[jsoff])^2
                end

                #startidx += data_properties[det_idx].baseline_lengths[baseline_idx]
                offset       += data_properties[det_idx].baseline_lengths[baseline_idx]
                sigma_offset += data_properties[det_idx].baseline_lengths[baseline_idx]

                baseline_idx += 1 #lpl: why ?
            end
        end

        #needed to assure that sum(baselines)==0

        if comm != nothing
            total_sum = MPI.Allreduce([sum(baselines)], MPI.SUM, comm)[1]
        else
            total_sum = sum(baselines)
        end

        baselines_sum .+= total_sum
end


function mpi_dot_prod(x, y; comm = nothing)

    @assert eltype(x) == eltype(y)
    result = zero(eltype(x))
    local_sum::eltype(x) = dot(x, y)

    if comm != nothing
        result = MPI.Allreduce([local_sum], MPI.SUM, comm)[1]
    else
        result = local_sum
    end

    return result
end

mutable struct DestripingResults
    threshold::Float64
    max_iter::Int
    convergence_param_list::Array{Float64}
    best_iteration::Int
    best_baselines::Array{Float64}
    best_sky_map::Any
    baseline_history::Any

    DestripingResults() = new(
        1e-9,
        1000,
        Float64[],
        -1,
        Float64[],
        nothing,
        nothing,
    )
end

function conj_grad(
    results::DestripingResults,
    baselines_sum,
    pix_idx,
    tod,
    num_of_pixels,
    twopsi,
    data_properties,
    num_of_baselines,
    rank;
    comm = nothing,
    save_baseline_history = false,
    callback = nothing,
    unseen = NaN,
    tod_mode = "sum",
    baselines_guess = nothing)

    T = eltype(tod)
    N = eltype(pix_idx)

    baselines = Array{T}(undef, num_of_baselines)
    r = Array{T}(undef, num_of_baselines)
    r_next = Array{T}(undef, num_of_baselines)
    p = Array{T}(undef, num_of_baselines)
    Ap = Array{T}(undef, num_of_baselines)
    
    if (baselines_guess != nothing)
        baselines .= baselines_guess
    else
        baselines .= 0    #starting baselines
    end
    convergence_parameter = zero(T)
    rdotr = zero(T)
    rdotr_next = zero(T)

    best_convergence_parameter = zero(T)
    best_baselines = zeros(T, num_of_baselines)
    results.best_iteration = 0

    r = baselines_sum - applya(baselines, pix_idx, num_of_baselines, num_of_pixels, twopsi, data_properties, comm = comm ,unseen = unseen ,tod_mode = tod_mode)

    #residual
    p .= r
    rdotr = mpi_dot_prod(r, r, comm = comm)
    best_convergence_parameter = sqrt(rdotr)

    results.baseline_history = save_baseline_history ? [] : nothing
    results.best_baselines = Array{Float64}(undef, num_of_baselines)
    results.best_baselines .= baselines

    rdotr == 0 && return results

    results.convergence_param_list = Float64[]
    iter_idx = 0
    while true
        Ap = applya(p, pix_idx,  num_of_baselines, num_of_pixels, twopsi,
                      data_properties, comm = comm ,unseen = unseen ,tod_mode = tod_mode)
            
        rdotr = mpi_dot_prod(r, r, comm = comm)
        pdotAp = mpi_dot_prod(p, Ap, comm = comm)

        alpha = rdotr / pdotAp
        @. baselines += alpha * p
        save_baseline_history && push!(results.baseline_history, copy(baselines))

        
        @. r_next = r - alpha * Ap
        rdotr_next = mpi_dot_prod(r_next, r_next, comm = comm)
        convergence_parameter = sqrt(rdotr_next)
        push!(results.convergence_param_list, convergence_parameter)

        if (convergence_parameter < best_convergence_parameter)
            best_convergence_parameter = convergence_parameter
            results.best_baselines .= baselines
            results.best_iteration = iter_idx
        end

        isnothing(callback) || callback(
            iter_idx = iter_idx,
            max_iter = results.max_iter,
            convergence_parameter = convergence_parameter,
            convergence_threshold = results.threshold,
        )
        
        ((convergence_parameter < results.threshold) || (iter_idx > results.max_iter)) && break

        beta = rdotr_next / rdotr
        #if (rank ==0)
        #    println("In ",rank," at ",iter_idx, "  Ap ",sum(Ap), "  rdotr ", rdotr,"  rdotr_next ",rdotr_next, "  pdotAp ",pdotAp, "  alpha ", alpha, "  beta ",beta )
        #end 
        @. p = r_next + beta * p
        r .= r_next
        iter_idx += 1
    end
end

#preallocates all buffers needed for internal calls
function conj_grad_prealloc(
    results::DestripingResults,
    baselines_sum,
    pix_idx,
    tod,
    num_of_pixels,
    twopsi,
    data_properties,
    num_of_baselines,
    rank;
    comm = nothing,
    save_baseline_history = false,
    callback = nothing,
    unseen = NaN,
    tod_mode = "sum",
    baselines_guess = nothing)

    T = eltype(tod)
    N = eltype(pix_idx)

    baselines = Array{T}(undef, num_of_baselines)
    r = Array{T}(undef, num_of_baselines)
    r_next = Array{T}(undef, num_of_baselines)
    p = Array{T}(undef, num_of_baselines)
    Ap = Array{T}(undef, num_of_baselines)
    
    map_buffer = zeros(T,(2,num_of_pixels))
    #rnr_buffer = zeros(T,(3,num_of_pixels))

    rnr_buffer = binned_noise_variance_mpi(pix_idx, num_of_pixels, twopsi, data_properties; comm = comm, unseen = unseen,
        tod_mode = tod_mode )
     

    #starting baselines
    if (baselines_guess != nothing)
        baselines .= baselines_guess
    else
        baselines .= 0
    end

    convergence_parameter = zero(T)
    rdotr = zero(T)
    rdotr_next = zero(T)

    best_convergence_parameter = zero(T)
    best_baselines = zeros(T, num_of_baselines)
    results.best_iteration = 0

    #r = baselines_sum - applya(baselines, pix_idx, num_of_baselines, num_of_pixels, twopsi, data_properties, comm = comm ,unseen = unseen ,tod_mode = tod_mode)
    applya!(baselines, pix_idx, num_of_baselines, num_of_pixels, twopsi, data_properties, Ap, map_buffer, rnr_buffer, comm = comm ,unseen = unseen ,tod_mode = tod_mode)
    r = baselines_sum - Ap

    #residual
    p .= r
    rdotr = mpi_dot_prod(r, r, comm = comm)
    best_convergence_parameter = sqrt(rdotr)

    results.baseline_history = save_baseline_history ? [] : nothing
    results.best_baselines = Array{Float64}(undef, num_of_baselines)
    results.best_baselines .= baselines

    rdotr == 0 && return results

    results.convergence_param_list = Float64[]
    iter_idx = 0
    while true
        #Ap = applya(p, pix_idx,  num_of_baselines, num_of_pixels, twopsi,
        #              data_properties, comm = comm ,unseen = unseen ,tod_mode = tod_mode)
            
        applya!(p, pix_idx,  num_of_baselines, num_of_pixels, twopsi,
            data_properties, Ap, map_buffer, rnr_buffer, comm = comm ,unseen = unseen ,tod_mode = tod_mode)

        rdotr = mpi_dot_prod(r, r, comm = comm)
        pdotAp = mpi_dot_prod(p, Ap, comm = comm)

        alpha = rdotr / pdotAp
        @. baselines += alpha * p
        save_baseline_history && push!(results.baseline_history, copy(baselines))

        
        @. r_next = r - alpha * Ap
        rdotr_next = mpi_dot_prod(r_next, r_next, comm = comm)
        convergence_parameter = sqrt(rdotr_next)
        push!(results.convergence_param_list, convergence_parameter)

        if (convergence_parameter < best_convergence_parameter)
            best_convergence_parameter = convergence_parameter
            results.best_baselines .= baselines
            results.best_iteration = iter_idx
        end

        isnothing(callback) || callback(
            iter_idx = iter_idx,
            max_iter = results.max_iter,
            convergence_parameter = convergence_parameter,
            convergence_threshold = results.threshold,
        )
        
        ((convergence_parameter < results.threshold) || (iter_idx > results.max_iter)) && break

        beta = rdotr_next / rdotr
        #if (rank ==0)
        #    println("In ",rank," at ",iter_idx, "  Ap ",sum(Ap), "  rdotr ", rdotr,"  rdotr_next ",rdotr_next, "  pdotAp ",pdotAp, "  alpha ", alpha, "  beta ",beta )
        #end 
        @. p = r_next + beta * p
        r .= r_next
        iter_idx += 1
    end
end


function destriped_map(baselines, pix_idx, tod, data_properties, num_of_pixels, twopsi,  num_of_baselines; comm = nothing, unseen = NaN, tod_mode = tod_mode)
    @assert length(tod) == length(pix_idx)
    @assert length(tod) == length(twopsi)

    tod2map_mpi(pix_idx, tod, num_of_pixels, twopsi, data_properties, comm = comm, tod_mode = tod_mode) - baseline2map_mpi(pix_idx, baselines, num_of_pixels, twopsi, data_properties, num_of_baselines, comm = comm, tod_mode = tod_mode)
end


@doc raw"""
    destripe(pix_idx, tod, num_of_pixels, data_properties, rank; comm = nothing, threshold = 1e-9, max_iter = 10000, save_baseline_history = false) -> DestripingResults

This MPI based function creates a map from a TOD and removes both 1/f
and white noise, using the destriping technique.

The parameters passed to the function have the following meaning:

- `pix_idx`: array containing the indices of the pixels visited by the instrument

- `tod`: the values measured by the polarimeters for each pixel
  (either I, Q, or U)

- `num_of_pixels`: the number of pixels in the map to be
   produced. This is used as an upper limit for the values in `pix_idx`

- `data_properties`: an array of structures `TodNoiseProperties`,
   holding information on each simulated polarimeter noise level,
   number and length of 1/f baselines and total number of samples.
   It can be obtained by using function `build_noise_properties`.

- `rank`: the rank of the current MPI process

- `comm`: the MPI communicator object.

The following arguments are optional:

- `threshold` is used by the conjugated-gradient algorithm. When the
   residual error of the iteration goes below this value, the iteration
   stops. The smaller the value, the more accurate the solution.

- `max_iter` is the maximum number of iterations to be executed in the
   conjugated-gradient algorithm. If the algorithm does not converge
   within this number of iterations, the process will quit without
   having reached the convergence threshold (see the `threshold`
   keyword above).

- If `save_baseline_history` is `true`, the return value will contain the
  sequence of baselines tested by the CG algorithm. Each MPI process will
  hold its own baselines.

The function returns a `DestripingResults` object containings the destriped
map, the sequence of baselines, and other information describing the
convergence of the CG algorithm.

# Remarks
- The length of the arrays `twopsi`, `pix_idx`, and `tod` must be the same;
- If you do not specify `comm`, no MPI will be used
"""
function destripe(
    pix_idx,
    tod,
    num_of_pixels,
    twopsi,
    data_properties,
    rank;
    comm = nothing,
    threshold = 1e-9,
    max_iter = 1_000,
    save_baseline_history = false,
    unseen = NaN,
    callback = nothing,
    tod_mode = "sum",
)

    @assert length(pix_idx) == length(tod)
    @assert length(twopsi)  == length(tod)

    num_of_baselines = 0
    for i in 1:length(data_properties)
        num_of_baselines += length(data_properties[i].baseline_lengths)
    end

    baselines_sum = applyz_and_sum(pix_idx,
        tod,
        num_of_pixels,
        twopsi,                             
        data_properties,
        num_of_baselines,
        comm = comm,
        unseen = unseen,
        tod_mode = tod_mode,
    )

    results = DestripingResults()
    results.threshold = threshold
    results.max_iter = max_iter

    conj_grad_prealloc(
        results,
        baselines_sum,
        pix_idx,
        tod,
        num_of_pixels,
        twopsi,
        data_properties,
        num_of_baselines,
        rank,
        comm = comm,
        save_baseline_history = save_baseline_history,
        callback = callback,
        unseen = unseen,
        tod_mode = tod_mode,
    )

    # once we have an estimate of the baselines, we can build the destriped map
    results.best_sky_map = destriped_map(results.best_baselines,
        pix_idx,
        tod,
        data_properties,
        num_of_pixels,
        twopsi,                                   
        num_of_baselines,
        comm = comm,
        unseen = unseen,
        tod_mode = tod_mode,                                   
    )

    results
end

#if we already have good guess for baselines, we input it to the destriper
function destripe(
    pix_idx,
    tod,
    num_of_pixels,
    twopsi,
    data_properties,
    rank;
    comm = nothing,
    threshold = 1e-9,
    max_iter = 1_000,
    save_baseline_history = false,
    unseen = NaN,
    callback = nothing,
    tod_mode = "sum",
    baselines_guess = nothing,
    )

    @assert length(pix_idx) == length(tod)
    @assert length(twopsi)  == length(tod)

    num_of_baselines = 0
    for i in 1:length(data_properties)
        num_of_baselines += length(data_properties[i].baseline_lengths)
    end

    baselines_sum = applyz_and_sum(pix_idx,
        tod,
        num_of_pixels,
        twopsi,                             
        data_properties,
        num_of_baselines,
        comm = comm,
        unseen = unseen,
        tod_mode = tod_mode,
    )

    results = DestripingResults()
    results.threshold = threshold
    results.max_iter = max_iter

    conj_grad_prealloc(
        results,
        baselines_sum,
        pix_idx,
        tod,
        num_of_pixels,
        twopsi,
        data_properties,
        num_of_baselines,
        rank,
        comm = comm,
        save_baseline_history = save_baseline_history,
        callback = callback,
        unseen = unseen,
        tod_mode = tod_mode,
        baselines_guess = baselines_guess,
    )

    # once we have an estimate of the baselines, we can build the destriped map
    results.best_sky_map = destriped_map(results.best_baselines,
        pix_idx,
        tod,
        data_properties,
        num_of_pixels,
        twopsi,                                   
        num_of_baselines,
        comm = comm,
        unseen = unseen,
        tod_mode = tod_mode,                                   
    )

    results
end

@doc raw"""
    baselines_covmat(polarimeters, sigma_k, baseline_length_s, fsamp_hz, total_time) -> covariance_matrix

This function produces the covariance matrix of the 1/f baselines computed by the destriper.
As an approximation, we ignore nondiagonal terms (i.e. the correlations between different baselines).
The function output is therefore an array corresponding to the covariance matrix diagonal.
Another assumption is that the white noise variance stays constant over a given baseline.

The baseline error is computed according to equation 29 in https://arxiv.org/abs/0904.3623

The parameters passed to the function have the following meaning:
    -`polarimeters`: the array of polarimeters ID numbers
    -`sigma_k`: the array of corresponding white noise sigma (in K)
    -`baseline_length_s` = the length (in s) of each 1/f baseline
    -`fsamp_hz` = the sampling frequency (in Hz)
    -`total_time` = the duration (in s) of the observation

The function firstly computes an Array of structures `TodNoiseProperties`,
which matches the white noise variance sigma with the number of 1/f baselines and their lengths, for each polarimeter.
Then it computes the covariance matrix according to these matches.
"""
function baselines_covmat(polarimeters, sigma_k, baseline_length_s, fsamp_hz, total_time)

    #needs to be rewritten to account for polarization
    ALL_data_properties = Array{TodNoiseProperties_Q}(undef, length(polarimeters))
    baselines_per_pol =  Int64(total_time / baseline_length_s)
    for i in 1:length(polarimeters)
        baseline_lengths  = repeat([baseline_length_s * fsamp_hz], baselines_per_pol)
        ALL_data_properties[i]  = TodNoiseProperties_Q(pol = polarimeters[i], 
            rms = sigma_k[i],
            baselines = baseline_lengths,
        )
    end

    covariance_matrix = []
    for i in 1:length(ALL_data_properties)
        partial_covmat = Array{Float64}(undef, length(ALL_data_properties[i].baseline_lengths))

        for j in 1:length(ALL_data_properties[i].baseline_lengths)
            sigma = ALL_data_properties[i].sigma
            this_baseline_length = ALL_data_properties[i].baseline_lengths[j]

            partial_covmat[j] = sigma^2 / this_baseline_length
        end
        covariance_matrix = append!(covariance_matrix, partial_covmat)
    end
    return covariance_matrix
end
