# HEAD

- (**Breaking change**) Turn `rank` and `comm` parameters into
  keywords for `generate_noise_mpi`
  [#39](https://github.com/lspestrip/Stripeline.jl/pull/39)
- (**Breaking change**) Use `nothing` instead of `missing` for MPI
  communicators
  [#39](https://github.com/lspestrip/Stripeline.jl/pull/39)
- Ensure that `generate_noise_mpi` produces different seeds for each
  polarimeter
  [#39](https://github.com/lspestrip/Stripeline.jl/pull/39)
- (**Breaking change**) Make `destripe` return a `DestripingResults`
  structure instead of a tuple
  [#39](https://github.com/lspestrip/Stripeline.jl/pull/39)
- Remove the dependency on `Quaternions.jl` [#38](https://github.com/lspestrip/Stripeline.jl/pull/38)
- Update the instrument database with the number of the polarizer for each horn [#37](https://github.com/lspestrip/Stripeline.jl/pull/37)
- Make the 1/f slope always in the range 0..2 [#34](https://github.com/lspestrip/Stripeline.jl/pull/34)

# Version 0.4.3

- Add functions to simulate ADCs [#24](https://github.com/lspestrip/Stripeline.jl/pull/24)
- Add a high-level API for accessing the instrument database [#30](https://github.com/lspestrip/Stripeline.jl/pull/30)
- Properly use the slope of the 1/f noise when generating maps from TODs [#31](https://github.com/lspestrip/Stripeline.jl/pull/31)
- Fix bandshapes for W-band polarimeters [#32](https://github.com/lspestrip/Stripeline.jl/pull/32)
- Plot `BandshapeInfo` and `SpectrumInfo` objects [#33](https://github.com/lspestrip/Stripeline.jl/pull/33)

# Version 0.4

- Conform coordinate systems to LSPE-STRIP-SP-017
  ([#30a0f](https://github.com/lspestrip/Stripeline.jl/commit/30a0fbdb5fe45fa20cd7a2fef08bc114ad3d7956),
  [#0418a](https://github.com/lspestrip/Stripeline.jl/commit/0418a40a489cd2dfd7607effe661c55af1ca649e))
- Add `day_duration_s` keyword to `genpointings`
- Add error bars to bandshape and always measure them as positive quantities ([#70dc6](https://github.com/lspestrip/Stripeline.jl/commit/70dc6612e3784e4b3cfded55540e01cccec0bbf3), [PR#22](https://github.com/lspestrip/Stripeline.jl/pull/22))
- Print the instrument database in Markdown format [#f2c2b](https://github.com/lspestrip/Stripeline.jl/commit/f2c2b11b317b149131ee4ab447a4ffe680148f2d)
- Add functions `sensitivity_tant`, `t_to_trj`, `trj_to_t`,
  `deltat_to_deltatrj, `deltatrj_to_deltat` [#23](https://github.com/lspestrip/Stripeline.jl/pull/23)

# Version 0.3.2

- Fix bugs in mapmaker.jl
  ([#f4349](https://github.com/lspestrip/Stripeline.jl/commit/f434916605201fd3e3daa81497248270b6378d76))


# Version 0.3.1

- Fix bug in generate_alldetectors_map.Jl


# Version 0.3.0

- New noise generation employs MPI to reduce wall clock times
- The seed for noise generation can now be specified as an input in `generate_noise_mpi`
