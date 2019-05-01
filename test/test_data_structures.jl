let
    # Test Orbit constructor
    epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0) 
    oe0  = [R_EARTH + 500e3, 0.0, 90.0, 0, 0, 0]
    eci0 = sOSCtoCART(oe0, use_degrees=true)

    orb = Orbit(epc0, epc0+86400.0, eci0, dtmax=30)

    @test orb.t[end]   == 86400.0
    @test orb.epc[end] == Epoch(2019, 1, 2, 12)
end

let
    # Test Orbit constructor
    epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0, tsys=:UTC) 
    oe0  = [R_EARTH + 500e3, 0.0, 90.0, 0, 0, 0]
    eci0 = sOSCtoCART(oe0, use_degrees=true)

    orb = Orbit(epc0, epc0+60.0, eci0, dtmax=1)

    @test orb.t[end]   == 60.0

    epc = epc0 + 0.5

    # Test interpolation
    eci = interpolate(orb, epc)

    for i in 1:6
        @test isapprox(eci[i], ((orb.eci[i, 2] - orb.eci[i, 1])*0.5 + orb.eci[i, 1]), atol=10)
    end
end

let
    img = Image(10, 20, look_angle_min=6, look_angle_max=40, reward=2.0)

    @test img.lon            == 10
    @test img.lat            == 20
    @test img.look_angle_min == 6
    @test img.look_angle_max == 40
    @test img.reward         == 2.0
    @test typeof(img.id)     == UUID

    img2 = Image(0, 0)

    @test img.id != img2.id
end

let
    # Test Image loading
    img_data = abspath(string(@__DIR__), "./data/landsat_test.json")
    images = load_images(img_data)

    @test length(images) == 150
end

let
    # Test Opportunity constructor
    sow  = Epoch(2018, 1, 1, 12, 1, 0, 0)
    eow  = Epoch(2018, 1, 1, 12, 1, 1, 0)
    opp1 = Opportunity(sow, eow)
    opp2 = Opportunity(sow, eow)

    @test opp1.id != opp2.id
    @test opp1.orbit == opp2.orbit
end

let
    # Test Opportunity Equality
    sow  = Epoch(2018, 1, 1, 12, 1, 0, 0)
    eow  = Epoch(2018, 1, 1, 12, 1, 1, 0)
    oppl = Opportunity(sow, eow)
    oppr = deepcopy(oppl)

    @test isequal(oppl, oppr) == true
    @test oppl == oppr
end