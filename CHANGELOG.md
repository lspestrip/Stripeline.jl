# HEAD

- Conform coordinate systems to LSPE-STRIP-SP-017
  ((https://github.com/lspestrip/Stripeline.jl/commit/30a0fbdb5fe45fa20cd7a2fef08bc114ad3d7956)[#30a0f],
  (https://github.com/lspestrip/Stripeline.jl/commit/0418a40a489cd2dfd7607effe661c55af1ca649e)[#0418a])


# Version 0.3.2

- Fix bugs in mapmaker.jl
  ((https://github.com/lspestrip/Stripeline.jl/commit/f434916605201fd3e3daa81497248270b6378d76)[#f4349])


# Version 0.3.1

- Fix bug in generate_alldetectors_map.Jl


# Version 0.3.0

- New noise generation employs MPI to reduce wall clock times
- The seed for noise generation can now be specified as an input in `generate_noise_mpi`
