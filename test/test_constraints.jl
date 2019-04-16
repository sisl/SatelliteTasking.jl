let
    # Simulation Horizon - 1 Day
    epc0 = Epoch(2019, 1, 1, 0, 0, 0, tsys=:UTC)
    epcf = Epoch(2019, 1, 1, 12, 0, 0, tsys=:UTC)

    oe   = [R_EARTH + 500e3, 0, 90.0, 0, 0, 0]
    eci0 = sOSCtoCART(oe, use_degrees=true)

    orb = Orbit(epc0, epcf, eci0, dtmax=30)

    # Load test images
    images = collect(load_images("../data/landsat_test.json")[[3, 5]])

    # Compute collect opportunities
    opportunities = find_all_opportunities(orb, images, sort=true)
    splitt_opps   = split_opportunities(opportunities, 1)

    opportunities = group_image_opportunities(splitt_opps)

    start_collect  = opportunities[images[1]][1][2]
    end_collect    = opportunities[images[end]][end][2]

    @test constraint_agility_single_axis(start_collect, end_collect) == true
    @test constraint_agility_single_axis(start_collect, start_collect) == false
end