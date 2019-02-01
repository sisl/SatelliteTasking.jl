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

    true_collects, perturbed_collects, mean_diff, sdev_diff  = compute_perturbed_collects(true_orbit, perturbed_orbits, images)

    @test typeof(true_collects)      == Array{Collect, 1}
    @test typeof(perturbed_collects) == Array{Array{Collect, 1}, 1}

    @debug mean_diff
    @debug sdev_diff
end