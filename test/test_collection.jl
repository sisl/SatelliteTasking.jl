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

# let
#     @time images = load_images("../data/landsat_test_600.json", collect_duration=5.0);

#     # Configure simulation
#     epc0 = Epoch(2019, 1, 1, 0, 0, 0, tsys=:UTC) # Start of time span
#     epcf = Epoch(2019, 1, 1, 6, 0, 0, tsys=:UTC) # End of simulation time span

#     # Set Simulation Time Step
#     timestep = 1
#     dtmax    = 5

#     # Define Satellite Orbit
#     oe   = [R_EARTH + 500*1e3, 0, 90.0, 0, 0, 0]
#     eci0 = sOSCtoCART(oe, use_degrees=true)

#     # Numer of perturbed orbits to simulate

#     # Set Perturbation Values 
#     pos_error = 5000 # Position knowledge error [m]
#     vel_error = 5    # Velocity knowledge error [m/s]
#     orb_mean  = zeros(Float64, 6)
#     orb_sdev  = vcat((pos_error/sqrt(3)*ones(Float64, 3))..., (vel_error/sqrt(3)*ones(Float64, 3))...)

#     # Simulate true and perturbed orbits
#     true_orbit, perturbed_orbits, eci_errors = simulate_orbits(0, epc0, epcf, eci0, orb_mean, orb_sdev, timestep=timestep, dtmax=dtmax);
    
#     true_opportunities, perturbed_opportunities, mean_diff, sdev_diff, missing_opportunities = compute_perturbed_opportunities(true_orbit, [true_orbit], images, epc_step=3600);

#     for i in 1:length(true_opportunities)
#         println("Opp pair:")
#         println(true_opportunities[i])
#         println(perturbed_opportunities[1][i])
#     end

#     println(sdev_diff[1, :])
#     println(sdev_diff[2, :])
#     println(sdev_diff[3, :])
# end

