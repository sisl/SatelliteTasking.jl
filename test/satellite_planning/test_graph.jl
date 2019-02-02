let
    # Simulation Horizon - 1 Day
    epc0 = Epoch(2019, 1, 1, 0, 0, 0, tsys=:UTC)
    epcf = Epoch(2019, 1, 1, 3, 0, 0, tsys=:UTC)

    oe   = [R_EARTH + 500e3, 0, 90.0, 0, 0, 0]
    eci0 = sOSCtoCART(oe, use_degrees=true)

    orb = Orbit(epc0, epcf, eci0, dtmax=30)

    # Load test images
    images = load_images("../data/landsat_test.json");

    # Compute collect opportunities
    true_collects = find_all_collects(orb, images, sort=true)
    
    # Build graph of feasible transitions
    # graph = sp_construct_graph(true_collects, [])

    # @test typeof(graph) == Dict{Collect, Array{Collect, 1}}
end