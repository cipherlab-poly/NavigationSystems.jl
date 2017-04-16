__precompile__()

module NavigationSystems

using StaticArrays
using CoordinateTransformations, Rotations
using PolarFact
using Geodesy
using NMEA

export
  # coordinates, complementing Geodesy.jl
  NED,
  NEDfromECEF, ECEFfromNED, NEDfromLLA, LLAfromNED, NEDfromUTMZ, UTMZfromNED,
  NEDfromUTM, UTMfromNED, NEDfromENU, ENUfromNED

  # rotation utilities
  closestOrthMatrix,

  # rotational kinematics
  crossmat, crossmat1,
  poissonUpdate,

  # sensors
  simulateBrownianMotion

include("ned.jl")
include("rotations.jl")
include("kinematics.jl")
include("sensors.jl")

end
