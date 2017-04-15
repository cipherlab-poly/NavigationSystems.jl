######################################
# Rotation matrices
######################################

"""
    O = closestOrthMatrix(M)
Compute the closest orthogonal matrix to M. This is guaranteed to provide the
closest rotation matrix if det(M) > 0.
"""
function closestOrthMatrix(M::Matrix)
  (polarfact(M)).U
end

# check this reference on closed form SVD for 2x2 matrices
# follosing still need to be sign correctd, returns a rot matrix always
# http://scicomp.stackexchange.com/questions/8899/robust-algorithm-for-2x2-svd
function closestOrthMatrix(M::SMatrix{2, 2, Float64})
  E = 0.5(M[1,1] + M[2,2]); F = 0.5(M[1,1] - M[2,2])
  G = 0.5(M[2,1] + M[1,2]); H = 0.5(M[2,1] - M[1,2])
  #Q = sqrt(E^2+H^2); R = sqrt(F^2+G^2)
  #sx = Q+R; sy = Q-R
  S = sign(E^2+H^2-F^2-G^2)
  a₁ = atan2(G, F); a₂ = atan2(H, E)
  Θ = 0.5(a₂-a₁); ϕ = 0.5(a₂+a₁)
  cϕ = cos(ϕ); sϕ = sin(ϕ)
  cΘ = cos(Θ); sΘ = sin(Θ)
  [cϕ -S*sϕ; sϕ S*cϕ] * [cΘ -sΘ; sΘ cΘ]
end

#=
# http://www.lucidarme.me/?p=4624
# See new references found for a probably better implementation with atan
# following function is probably wrong
function closestOrthMatrix2(M::SMatrix{2, 2, Float64})
  Θ = atan2(2M[1,1]*M[2,1]+2M[1,2]*M[2,2], M[1,1]^2+M[1,2]^2-M[2,1]^2-M[2,2]^2)
  cΘ = cos(Θ); sΘ = sin(Θ)
  ϕ = atan2(2M[1,1]*M[1,2]+2M[2,1]*M[2,2], M[1,1]^2-M[1,2]^2+M[2,1]^2-M[2,2]^2)
  cϕ = cos(ϕ); sϕ = sin(ϕ)
  u = (M[1,1]*cΘ + M[2,1]*sΘ)*cϕ + (M[1,2]*cΘ + M[2,2]*sΘ)*sϕ
  v = (M[1,1]*sΘ - M[2,1]*cΘ)*sϕ + (-M[1,2]*sΘ + M[2,2]*cΘ)*cϕ
  su = sign(u); sv = sign(v)
  [cΘ -sΘ; sΘ cΘ] * [su*cϕ -sv*sϕ; su*sϕ sv*cϕ]
end
=#
