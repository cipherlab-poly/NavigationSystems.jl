@testset "Navigation Systems Tests" begin
  @testset "Coordinates" begin
    c = NED(1.0,2,-3.0)

    @test c.n === 1.0
    @test c.e === 2.0
    @test c.d === -3.0
  end

  @testset "Rotation matrices" begin
    M = [-0.539715 2.39379; -0.872195 -0.371752]
    M1 = @SMatrix [-0.539715 2.39379; -0.872195 -0.371752]

    @test closestOrthMatrix(M) ≈ closestOrthMatrix(M1)
  end

  @testset "Rotational kinematics" begin
    ω = @SVector [2.1, 9.8, 4]
    ω1 = [2.1, 9.8, 4]
    ω_cross = @SMatrix [0.0 -4.0 9.8; 4.0 0.0 -2.1; -9.8 2.1 0.0]
    ω_cross1 = [0.0 -4.0 9.8; 4.0 0.0 -2.1; -9.8 2.1 0.0]

    @test crossmat(ω1) === ω_cross
    @test crossmat(ω) === ω_cross
    @test crossmat1(ω1) == ω_cross1
    @test crossmat1(ω) == ω_cross1

    r2 = RotZ(pi/7)
    @test poissonUpdate(r2, [0, 0, pi/8], 0.5) == RotZ(pi/7+pi/16)


    r1 = RotMatrix{3}([0.121957  0.261358  0.957506;
                       0.758496 -0.646757  0.0799281;
                       0.640164  0.716517 -0.277116])
    @test_approx_eq_eps poissonUpdate(r1, ω, 0.04) [-0.201498  0.347741  0.915682;
      0.550939 -0.732715  0.399492; 0.809855  0.584982 -0.0439436] 1e-6
  end
end
