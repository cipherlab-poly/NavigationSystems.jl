# Navigation.jl
A Julia package supporting the design of multi-sensor navigation systems


## Set-up ##

```julia
Pkg.clone("https://github.com/cipherlab-poly/Navigation.jl.git")
using Navigation
```


## Run Tests ##

```julia
Pkg.test("Navigation")
```


## Dependencies

- [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) for faster computations with small matrices. This implies that only Julia versions ≥ 0.5 will be supported.
- [Rotations.jl](https://github.com/FugroRoames/Rotations.jl) for various representations of rotations. 
- [CoordinateTransformations.jl](https://github.com/FugroRoames/CoordinateTransformations.jl).
- [Geodesy.jl](https://github.com/JuliaGeo/Geodesy.jl) for the conversion between standard coordinate systems (LLA, ECEF, etc.), which leverages the CoordinateTransformations.jl package.
- [PolarFact.jl](https://github.com/weijianzhang/PolarFact.jl) for matrix polar decomposition, which provides the solution to certain "best" rotation approximations.
- [NMEA.jl](https://github.com/furface/NMEA.jl) to parse GPS messages.
- [Plots.jl](https://github.com/JuliaPlots/Plots.jl) at least during initial development.