# NavigationSystems.jl

A Julia package supporting the design of multi-sensor integrated navigation systems


## Set-up ##

This package is in development and currently not working.

This package is intended to work with Julia v1.0 and later.
To use with Julia prior to v1.0, please see the branch v0.5. 

To add the package, in the package manager (press ] in the REPL)
```julia
pkg> add "https://github.com/cipherlab-poly/NavigationSystems.jl.git"
```

To use the package functionalities
```julia
using NavigationSystems
```

## Run Tests ##

```julia
pkg> test NavigationSystems
```


## Dependencies

- [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) for faster computations with small matrices. This implies that only Julia versions ≥ 0.5 will be supported.
- [Rotations.jl](https://github.com/FugroRoames/Rotations.jl) for various representations of rotations. 
- [CoordinateTransformations.jl](https://github.com/FugroRoames/CoordinateTransformations.jl).
- [Geodesy.jl](https://github.com/JuliaGeo/Geodesy.jl) for the conversion between standard coordinate systems (LLA, ECEF, etc.), which leverages the CoordinateTransformations.jl package.
- [Plots.jl](https://github.com/JuliaPlots/Plots.jl) at least during initial development.

The following dependencies have been removed for now
- [PolarFact.jl](https://github.com/weijianzhang/PolarFact.jl) for matrix polar decomposition, which provides the solution to certain "best" rotation approximations. Does not seem to work with Julia >= 1.0
- [NMEA.jl](https://github.com/furface/NMEA.jl) to parse GPS messages.
