let
    sc = Spacecraft()

    @test sc.power_generation == 0.0
    @test sc.power_draw_downlink == 0.0
    @test sc.data_generation_backorbit == 0.0
    @test sc.data_generation_imaging == 0.0
    @test sc.data_downlink == 0.0
end

let
    sc = Spacecraft(
        power_generation = 1.0,
        power_draw_downlink = 0.4,
        data_generation_backorbit = 0.3,
        data_generation_imaging = 0.2,
        data_downlink = 0.1,
    )

    problem = MDPProblem(spacecraft=sc)

    @test problem.spacecraft.power_generation == 1.0
    @test problem.spacecraft.power_draw_downlink == 0.4
    @test problem.spacecraft.data_generation_backorbit == 0.3
    @test problem.spacecraft.data_generation_imaging == 0.2
    @test problem.spacecraft.data_downlink == 0.1
end

let
    request = Request(id=1)

    @test request.id == 1
end