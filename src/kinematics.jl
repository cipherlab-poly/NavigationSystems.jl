######################################
# Rotational kinematics
######################################

"""
    m = crossmat(a)
Input: a is a 3 dimensional vector

Output: `m`, a static array, so that `typeof(m) = SMatrix{3, 3, eltype(a))`.
`m` is the 3x3 matrix representing the cross product operator a x .

Use `m = crossmat1(a)` instead if you want an output in the form of a standard matrix,
without the overhead of conversion from a static array.
"""
function crossmat(a::AbstractVector)
  z = zero(eltype(a))
  a_x = @SMatrix [  z   -a[3]  a[2];
                   a[3]   z   -a[1];
                  -a[2]  a[1]   z   ]
end

#= # is this necessary?
function crossmat(a::SVector{3})
  z = zero(eltype(a))
  a_x = @SMatrix [  z   -a[3]  a[2];
                   a[3]   z   -a[1];
                  -a[2]  a[1]   z   ]
end
=#

#This function returns instead a standard matrix
function crossmat1(a::AbstractVector)
    z = zero(eltype(a))
    a_x = [  z   -a[3]  a[2];
            a[3]   z   -a[1];
           -a[2]  a[1]   z   ]
end

#= # is this necessary?
function crossmat1(a::SVector{3})
    z = zero(eltype(a))
    a_x = [  z   -a[3]  a[2];
            a[3]   z   -a[1];
           -a[2]  a[1]   z   ]
end
=#

# Attitude Estimation - Second order Poisson integration schemes
# ω : angular velocity of b wrt a expressed in b
# R = R^a_b; Rdot = R ω
function poissonUpdate(R::Rotation{3}, ω::AbstractVector, dt)
  ω_x = crossmat(ω)
  nrm = norm(ω)
  dtn = nrm * dt
  R * ( eye(SMatrix{3,3, Float64}) + ω_x * sin(dtn)/nrm +
        ω_x^2 * (1-cos(dtn))/nrm^2 )
end

function poissonUpdate(R::RotMatrix{3}, ω::AbstractVector, dt)
  ω_x = crossmat(ω)
  nrm = norm(ω)
  dtn = nrm * dt
  R * ( eye(SMatrix{3,3, Float64}) + ω_x * sin(dtn)/nrm +
        ω_x^2 * (1-cos(dtn))/nrm^2 )
end

# q : q^a_b; qdot = 0.5 * q * ω
# For a reference on quaternion integration, see for example
# https://link.springer.com/article/10.1007/s00707-013-0914-2
function poissonUpdate(q::Quat, ω::AbstractVector, dt)
  nrm = norm(ω)
  dtw = 0.5 * dt * nrm
  c = cos(dtw)
  s = sin(dtw)/nrm
  #qmult = Quat(c, s*ω[1], s*ω[2], s*ω[3])
  #qmult * q
  Quat( c*q.w - s*(ω[1]*q.x + ω[2]*q.y + ω[3]*q.z),
      c*q.x + s*(ω[1]*q.w + ω[2]*q.z - ω[3]*q.y),
      c*q.y - s*(ω[1]*q.z - ω[2]*q.w - ω[3]*q.x),
      c*q.z + s*(ω[1]*q.y - ω[2]*q.x + ω[3]*q.w) )
end

# Sequence of Poisson updates - Do not save trajectory
function poissonUpdate(q::Quat, ωs::Vector{SVector{3,Float64}}, dt)
  x = MVector(q.w, q.x, q.y, q.z)  # mutable, to store the current quaternion
  nrm = zero(Float64)
  dtw = zero(Float64)
  c = zero(Float64)
  s = zeros(Float64)
  for ω ∈ ωs
    nrm = norm(ω)
    dtw = 0.5 * dt * nrm
    c = cos(dtw)
    s = sin(dtw)/nrm
    (x[1], x[2], x[3], x[4]) =
      ( c*x[1] - s*(ω[1]*x[2] + ω[2]*x[3] + ω[3]*x[4]),
        c*x[2] + s*(ω[1]*x[1] + ω[2]*x[4] - ω[3]*x[3]),
        c*x[3] - s*(ω[1]*x[4] - ω[2]*x[1] - ω[3]*x[2]),
        c*x[4] + s*(ω[1]*x[3] - ω[2]*x[2] + ω[3]*x[1]) )
  end
  Quat(x[1],x[2],x[3],x[4])  # return the final quaternion
end

#=
# Less efficient computation
function poissonUpdate(q::Quat, ω::AbstractVector, dt)
  nrm = norm(ω)
  dtw = 0.5 * dt * nrm
  c = cos(dtw)
  s = sin(dtw)/nrm
  q_im = SVector(q.x, q.y, q.z)
  tmp = s * q.w * ω + (c * eye(SMatrix{3,3, Float64}) - s *  crossmat(ω)) * q_im
  Quat(
    c * q.w - s * (ω[1]*q.x + ω[2]*q.y + ω[3]*q.z), tmp[1], tmp[2], tmp[3]
  )
end

# ω : angular velocity of a wrt b expressed in a
# R = R^a_b; Rdot = - ω x R
function poissonUpdateBis(ω::AbstractVector, R::RotMatrix{3}, dt)
  ω_x = crossmat(ω)
  nrm = norm(ω)
  dtn = nrm * dt
  R1 = ( eye(SMatrix{3,3, Float64}) - ω_x * sin(dtn)/nrm +
        ω_x^2 * (1-cos(dtn))/nrm^2 ) * R
  R1
end
=#
