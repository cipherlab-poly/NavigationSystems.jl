##############################################################
# NED Frames of Reference, complementing Geodesy.jl ENU frames
##############################################################

# NED frame related definitions
# Based on Geodesy.ENU
# Currently we essentially copy the functions from Geodesy.jl
# Should probably be done more efficiently using metaprogramming

# See Geodesy.jl --> points.jl
const FieldVector{T} = StaticArrays.FieldVector{3, T}

"""
    NED(n, e, d = 0.0)
North-East-Down (NED) coordinates. A local Cartesian coordinate system.
Defined similarly to Geodesy.ENU
"""
struct NED{T <: Number} <: FieldVector{T}
    n::T
    e::T
    d::T
end
NED(x :: T, y :: T) where {T} = NED(x, y, zero(T))
@inline function NED(x,y,z)
    T = promote_type(promote_type(typeof(x),typeof(y)), typeof(z))
    NED{T}(x,y,z)
end
Base.show(io::IO, ::MIME"text/plain", ned::NED) = print(io, "NED($(ned.n), $(ned.e), $(ned.d))")


# See Geodesy.jl --> conversion.jl
NED(ned::NED, datum) = ned
NED(ecef::ECEF, origin, datum) = NEDfromECEF(origin, datum)(ecef)
ECEF(ned::NED, origin, datum) = ECEFfromNED(origin, datum)(ned)
NED(lla::LLA, origin, datum) =  NEDfromLLA(origin, datum)(lla)
LLA(ned::NED, origin, datum) = LLAfromNED(origin, datum)(ned)
NED(utm::UTMZ, origin, datum) = NEDfromUTMZ(origin, datum)(utm)
UTMZ(ned::NED, origin, datum) = UTMZfromNED(origin, datum)(ned)
NED(utm::UTM, zone::Integer, hemisphere::Bool, origin, datum) =
    NEDfromUTM(origin, zone, hemisphere, datum)(utm)
UTM(ned::NED, zone::Integer, hemisphere::Bool, origin, datum) =
    UTMfromNED(origin, zone, hemisphere, datum)(ned)
# Conversion between NED and ENU frames assumed to be centered at the same origin
NED(enu::ENU) = NEDfromENU(enu)
ENU(ned::NED) = ENUfromNED(ned)

# See Geodesy.jl --> transformations.jl

# Conversion between NED and ENU frames assumed to be centered at the same origin
struct NEDfromENU <: Transformation
end
NEDfromENU(enu::ENU) = NED(enu.n, enu.e, -enu.u)
struct ENUfromNED <: Transformation
end
ENUfromNED(ned::NED) = ENU(ned.e, ned.n, -ned.d)

"""
    NEDfromECEF(origin, datum)
    NEDfromECEF(origin::UTM, zone, isnorth, datum)
    NEDfromECEF(origin::ECEF, lat, lon)
Construct a `Transformation` object to convert from global `ECEF` coordinates
to local `NED` coordinates centered at the `origin`. This object pre-caches both the
ECEF coordinates and latitude and longitude of the origin for maximal efficiency.
"""
struct NEDfromECEF{T} <: Transformation
    origin::ECEF{T}
    lat::T
    lon::T
end

NEDfromECEF(origin::LLA, datum) = NEDfromECEF(ECEFfromLLA(datum)(origin), origin.lat, origin.lon)
function NEDfromECEF(origin::ECEF, datum)
    origin_lla = LLAfromECEF(datum)(origin)
    NEDfromECEF(origin, origin_lla.lat, origin_lla.lon)
end
Base.show(io::IO, trans::NEDfromECEF) = print(io, "NEDfromECEF($(trans.origin), lat=$(trans.lat)°, lon=$(trans.lon)°)")
Base.isapprox(t1::NEDfromECEF, t2::NEDfromECEF; kwargs...) = isapprox(t1.origin, t2.origin; kwargs...) && isapprox(t1.lat, t2.lat; kwargs...) && isapprox(t1.lon, t2.lon; kwargs...)

function (trans::NEDfromECEF)(ecef::ECEF)
    ϕdeg, λdeg = trans.lat, trans.lon

    ∂x = ecef.x - trans.origin.x
    ∂y = ecef.y - trans.origin.y
    ∂z = ecef.z - trans.origin.z

    # Compute rotation matrix
    sinλ, cosλ = sind(λdeg), cosd(λdeg)
    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)

    # R = [     -sinλ       cosλ  0.0
    #      -cosλ*sinϕ -sinλ*sinϕ cosϕ
    #       cosλ*cosϕ  sinλ*cosϕ sinϕ]
    #
    # east, north, up = R * [∂x, ∂y, ∂z]
    east  = ∂x * -sinλ      + ∂y * cosλ       + ∂z * 0.0
    north = ∂x * -cosλ*sinϕ + ∂y * -sinλ*sinϕ + ∂z * cosϕ
    up    = ∂x * cosλ*cosϕ  + ∂y * sinλ*cosϕ  + ∂z * sinϕ

    return NED(north, east, -up)
end

"""
    ECEFfromNED(origin, datum)
    ECEFfromNED(origin::UTM, zone, isnorth, datum)
    ECEFfromNED(origin::ECEF, lat, lon)
Construct a `Transformation` object to convert from local `NED` coordinates
centred at `origin` to global `ECEF` coodinates. This object pre-caches both the
ECEF coordinates and latitude and longitude of the origin for maximal efficiency.
"""
struct ECEFfromNED{T} <: Transformation
    origin::ECEF{T}
    lat::T
    lon::T
end
ECEFfromNED(origin::LLA, datum) = ECEFfromNED(ECEFfromLLA(datum)(origin), origin.lat, origin.lon)
function ECEFfromNED(origin::ECEF, datum)
    origin_lla = LLAfromECEF(datum)(origin)
    ECEFfromNED(origin, origin_lla.lat, origin_lla.lon)
end
Base.show(io::IO, trans::ECEFfromNED) = print(io, "ECEFfromNED($(trans.origin), lat=$(trans.lat)°, lon=$(trans.lon)°)")
Base.isapprox(t1::ECEFfromNED, t2::ECEFfromNED; kwargs...) = isapprox(t1.origin, t2.origin; kwargs...) && isapprox(t1.lat, t2.lat; kwargs...) && isapprox(t1.lon, t2.lon; kwargs...)

function (trans::ECEFfromNED)(ned::NED)
    ϕdeg, λdeg = trans.lat, trans.lon

    # Compute rotation matrix
    sinλ, cosλ = sind(λdeg), cosd(λdeg)
    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)

    # Rᵀ = [-sinλ -cosλ*sinϕ cosλ*cosϕ
    #        cosλ -sinλ*sinϕ sinλ*cosϕ
    #         0.0       cosϕ      sinϕ]
    # Δx, Δy, Δz = Rᵀ * [east, north, up]
    Δx = -sinλ * ned.e + -cosλ*sinϕ * ned.n - cosλ*cosϕ * ned.d
    Δy =  cosλ * ned.e + -sinλ*sinϕ * ned.n - sinλ*cosϕ * ned.d
    Δz =   0.0 * ned.e +       cosϕ * ned.n -      sinϕ * ned.d

    X = trans.origin.x + Δx
    Y = trans.origin.y + Δy
    Z = trans.origin.z + Δz

    return ECEF(X,Y,Z)
end

Base.inv(trans::ECEFfromNED) = NEDfromECEF(trans.origin, trans.lat, trans.lon)
Base.inv(trans::NEDfromECEF) = ECEFfromNED(trans.origin, trans.lat, trans.lon)

"""
    NEDfromLLA(origin, datum)
Creates composite transformation `NEDfromECEF(origin, datum) ∘ ECEFfromLLA(datum)`.
"""
NEDfromLLA(origin, datum) = NEDfromECEF(origin, datum) ∘ ECEFfromLLA(datum)

"""
    LLAfromNED(origin, datum)
Creates composite transformation `LLAfromECEF(datum) ∘ ECEFfromNED(origin, datum)`.
"""
LLAfromNED(origin, datum) = LLAfromECEF(datum) ∘ ECEFfromNED(origin, datum)

# UTMZ conversions

NEDfromECEF(origin::UTMZ, datum) = NEDfromECEF(LLAfromUTMZ(datum)(origin), datum)
ECEFfromNED(origin::UTMZ, datum) = ECEFfromNED(LLAfromUTMZ(datum)(origin), datum)

"""
    NEDfromUTMZ(origin, datum)
Creates composite transformation `NEDfromLLA(origin, datum) ∘ LLAfromUTMZ(datum)`.
"""
NEDfromUTMZ(origin, datum) = NEDfromLLA(origin, datum) ∘ LLAfromUTMZ(datum)

"""
    UTMZfromNED(origin, datum)
Creates composite transformation `UTMZfromLLA(datum) ∘ LLAfromNED(origin, datum)`.
"""
UTMZfromNED(origin, datum) = UTMZfromLLA(datum) ∘ LLAfromNED(origin, datum)

NEDfromECEF(origin::UTM, zone::Integer, isnorth::Bool, datum) = NEDfromECEF(LLAfromUTM(zone, isnorth, datum)(origin), datum)
ECFfromNED(origin::UTM, zone::Integer, isnorth::Bool, datum) = ECEFfromNED(LLAfromUTM(zone, isnorth, datum)(origin), datum)

# Assume origin and utm point share the same zone and hemisphere
UTMfromNED(origin::UTM, zone::Integer, isnorth::Bool, datum) = UTMfromLLA(zone, isnorth, datum) ∘ LLAfromNED(UTMZ(origin, zone, isnorth), datum)
NEDfromUTM(origin::UTM, zone::Integer, isnorth::Bool, datum) = NEDfromLLA(UTMZ(origin, zone, isnorth), datum) ∘ LLAfromUTM(zone, isnorth, datum)

"""
    UTMfromNED(origin, zone, isnorth, datum)
Creates composite transformation `UTMfromLLA(zone, isnorth, datum) ∘ LLAfromNED(origin, datum)`.
If `origin` is a `UTM` point, then it is assumed it is in the given specified zone and hemisphere.
"""
UTMfromNED(origin, zone::Integer, isnorth::Bool, datum) = UTMfromLLA(zone, isnorth, datum) ∘ LLAfromNED(origin, datum)

"""
    NEDfromUTM(origin, zone, isnorth, datum)
Creates composite transformation `UTMfromLLA(zone, isnorth, datum) ∘ LLAfromNED(origin, datum)`.
If `origin` is a `UTM` point, then it is assumed it is in the given specified zone and hemisphere.
"""
NEDfromUTM(origin, zone::Integer, isnorth::Bool, datum) = NEDfromLLA(origin, datum) ∘ LLAfromUTM(zone, isnorth, datum)
