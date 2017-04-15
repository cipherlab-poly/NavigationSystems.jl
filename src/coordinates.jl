######################################
# Frames of Reference and coordinates
######################################

"""
    NED(n, e, d = 0.0)
North-East-Down (NED) coordinates. A local Cartesian coordinate system.
Defined similarly to Geodesy.ENU
"""
immutable NED{T <: Number} <: FieldVector{T}
    n::T
    e::T
    d::T
end
NED{T}(x :: T, y :: T) = NED(x, y, zero(T))
@inline function NED(x,y,z)
    T = promote_type(promote_type(typeof(x),typeof(y)), typeof(z))
    NED{T}(x,y,z)
end
Base.show(io::IO, ::MIME"text/plain", ned::NED) =
  print(io, "NED($(ned.e), $(ned.n), $(ned.d))")
