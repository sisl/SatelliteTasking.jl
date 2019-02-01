let
    # Test image view
    ecef  = sGEODtoECEF([0.0, 0.0, R_EARTH + 500e3], use_degrees=true)
    image = Image(0, 0)

    look_angle, range = image_view_geometry(ecef, image)

    @test look_angle == 0.0
    @test isapprox(range, R_EARTH+500e3, atol=1e-9)
end

let
    ecef  = sGEODtoECEF([0.0, 0.0, R_EARTH + 500e3], use_degrees=true)
    image = Image(0, 0, look_angle_min=0.0)
    @test image_visible(ecef, image) == true

    ecef  = sGEODtoECEF([180.0, 0.0, R_EARTH + 500e3], use_degrees=true)
    image = Image(0, 0)
    @test image_visible(ecef, image) == false
end

let
    indices = [1, 2, 3, 47, 48, 50]

    gidx = SatelliteTasking.Collection.group_indices(indices)

    @test length(gidx[1]) == 3
    array_isapprox(gidx[1], [1, 2, 3])

    @test length(gidx[2]) == 2
    array_isapprox(gidx[2], [47, 48])

    @test length(gidx[3]) == 1
    array_isapprox(gidx[3], [50])
end

