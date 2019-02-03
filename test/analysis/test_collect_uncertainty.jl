let
    # Simulation Horizon - 1 Day
    epc0 = Epoch(2019, 1, 1, 0, 0, 0, tsys=:UTC)
    epcf = Epoch(2019, 1, 1, 1, 0, 0, tsys=:UTC)

    oe   = [R_EARTH + 500e3, 0, 90.0, 0, 0, 0]
    eci0 = sOSCtoCART(oe, use_degrees=true)

    pmean = zeros(Float64, 6)
    pstd  = zeros(Float64, 6)

    num_orbits = 3

    # pstd  = [5000/sqrt(3), 5000/sqrt(3), 5000/sqrt(3), 10/sqrt(3), 10/sqrt(3), 10/sqrt(3)]
    pstd = ones(Float64, 6)
    true_orbit, perturbed_orbits, eci_errors = simulate_orbits(num_orbits, epc0, epcf, eci0, pmean, pstd, timestep=1)

    # Load test images
    images = load_images("../data/landsat_test.json");

    true_opportunities, perturbed_opportunities, mean_diff, sdev_diff, miss_opp = compute_perturbed_opportunities(true_orbit, perturbed_orbits, images)

    @test typeof(true_opportunities)      == Array{Opportunity, 1}
    @test typeof(perturbed_opportunities) == Array{Array{Opportunity, 1}, 1}

    @debug mean_diff
    @debug sdev_diff
end

let
    opp1 = Opportunity(Epoch(2018, 1, 1, 12, 1, 0, 0), Epoch(2018, 1, 1, 12, 2, 0, 0))
    opp2 = Opportunity(Epoch(2018, 1, 1, 13, 1, 0, 0), Epoch(2018, 1, 1, 13, 2, 0, 0))

    missing_from_a, missing_from_b = find_missing_opportunity([opp1], [opp2])

    @test length(missing_from_a) == 1
    @test missing_from_a[1]      == opp2

    @test length(missing_from_b) == 1
    @test missing_from_b[1]      == opp1
end

let
    opp = Opportunity(Epoch(2018, 1, 1, 12, 0, 0, 0), Epoch(2018, 1, 1, 12, 0, 30, 0), dwell_time=10)

    collects = compute_collects_by_number([opp])

    @test length(collects) == 3
end