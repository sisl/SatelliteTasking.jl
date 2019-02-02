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

    true_collects, perturbed_collects, mean_diff, sdev_diff, miss_collect = compute_perturbed_collects(true_orbit, perturbed_orbits, images)

    @test typeof(true_collects)      == Array{Collect, 1}
    @test typeof(perturbed_collects) == Array{Array{Collect, 1}, 1}

    @debug mean_diff
    @debug sdev_diff
end

let
    col1 = Collect(Epoch(2018, 1, 1, 12, 1, 0, 0), Epoch(2018, 1, 1, 12, 2, 0, 0))
    col2 = Collect(Epoch(2018, 1, 1, 12, 3, 0, 0), Epoch(2018, 1, 1, 12, 4, 0, 0))

    missing_from_a, missing_from_b = find_missing_collect([col1], [col2])

    @test length(missing_from_a) == 1
    @test missing_from_a[1]      == col2

    @test length(missing_from_b) == 1
    @test missing_from_b[1]      == col1
end