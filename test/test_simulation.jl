let
    # Simulation Horizon - 1 Day
    epc0 = Epoch(2019, 1, 1, 0, 0, 0, tsys=:UTC)
    epcf = Epoch(2019, 1, 1, 0, 15, 0, tsys=:UTC)

    oe   = [R_EARTH + 500e3, 0, 90.0, 0, 0, 0]
    eci0 = sOSCtoCART(oe, use_degrees=true)

    pmean = zeros(Float64, 6)
    pstd  = zeros(Float64, 6)

    num_orbits = 2

    @test_throws ArgumentError simulate_orbits(num_orbits, epc0, epcf, eci0, pmean, pstd, timestep=30)

    pstd  = ones(Float64, 6)
    true_orbit, perturbed_orbits, eci_errors = simulate_orbits(num_orbits, epc0, epcf, eci0, pmean, pstd, timestep=30)

    @test typeof(true_orbit)          == Orbit
    @test typeof(perturbed_orbits)    == Array{Orbit, 1}
    @test typeof(perturbed_orbits[1]) == Orbit
    @test length(perturbed_orbits)    == num_orbits
    @test all(eci_errors .!= 0)
end