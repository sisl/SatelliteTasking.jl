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
    img_data = abspath(string(@__DIR__), "../data/landsat_test.json")
    images = load_images(img_data)

    @test length(images) == 150
end

let
    # Test Collect constructor
    sow  = Epoch(2018, 1, 1, 12, 1, 0, 0)
    eow  = Epoch(2018, 1, 1, 12, 1, 1, 0)
    col1 = Collect(sow, eow)
    col2 = Collect(sow, eow)

    @test col1.id != col2.id
    @test col1.orbit_id == col2.orbit_id
end