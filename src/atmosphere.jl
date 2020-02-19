export initialize_kolmogorov, atm_update, atm_observe

using FFTW
using Random

@doc raw"""
        Some documentation...
"""
function initialize_kolmogorov(N::Int64,
                               seed::Int64,
                               scale_factor::Float64,
                               corr_length::Float64
                              )

        rng = MersenneTwister(seed)
        random_realization = randn(rng, N, N, N)

        inds = ((0:(N-1)) .+ 0.5 .- (N/2))
        X = reshape(repeat(inds, N*N), N, N, N)
        Y = reshape(repeat(reshape(X[:,:,1]', N*N), N), N, N, N)
        Z = reshape(sort(reshape(repeat(inds, N*N)', N*N*N)), N, N, N)

        R = ((X .^2 .+ Y .^2 .+ Z .^2).^0.5) * scale_factor
        a = 1 / corr_length
        mag_k = ((2.0 * π) ./ (R .+ 0.01 )).^(5.0 / 3.0) .* exp.(- ((R.^2) ./ (2*a^2)) )

        return rng, random_realization, mag_k

end


function atm_update(random_realization::Array{Float64, 3},
                    rng::Any,
                    N::Int64,
                    w_speed::Float64,
                    w_dir::Float64,
                    scale_factor::Float64,
                    f_samp::Float64)

    # Angle from the north eastwards
    vx = w_speed * cos(deg2rad(w_dir))
    vy = w_speed * sin(deg2rad(w_dir))

    nx = Int(round((scale_factor / vx) * f_samp))
    ny = Int(round((scale_factor / vy) * f_samp))


    random_realization[2:end, :, :] .= random_realization[1:end-1, :, :]
    random_realization[1, :, :] .= randn(rng, N, N)

    return random_realization

end

function atm_observe(random_realization::Array{Float64, 3},
                     kol_spec::Array{Float64, 3},
                     az::Float64,
                     el::Float64)

        fft_rand = fft(random_realization)
        atm_map = ifft(fft_rand .* fftshift(kol_spec))
        atm_map_norm = real.(atm_map) ./ maximum(abs.(atm_map))

        # Integration along the line of sight of each detector
	 	# and evaluate the covariance matrix

        return atm_map_norm
end
