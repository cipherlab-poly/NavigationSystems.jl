__precompile__()

module Navigation

using StaticArrays
using CoordinateTransformations, Rotations
using PolarFact
using Geodesy
using NMEA

export
  # coordinates, complementing Geodesy.jl
  NED,

  # rotation utilities
  closestOrthMatrix,

  # rotational kinematics
  crossmat, crossmat1,
  poissonUpdate,

  # sensors
  simulateBrownianMotion

include("coordinates.jl")
include("rotations.jl")
include("kinematics.jl")
include("sensors.jl")

end
