using Test
using Random
using LinearAlgebra
using OrdinaryDiffEq
using Logging
using UUIDs
using SatelliteDynamics

# Package Under Test
using SatelliteTasking
using SatelliteTasking.SatellitePlanning
using SatelliteTasking.Analysis

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
    opportunities = find_all_opportunities(orb, images, sort=true)
    collects      = compute_collects_by_number([opportunities[end]], 1)

    @test opportunities[end].sow <= collects[end].sow <= opportunities[end].eow
    @test opportunities[end].sow <= collects[end].eow <= opportunities[end].eow

    collects = compute_collects_by_number([opportunities[end]], 3)
    for collect in collects
        @test opportunities[end].sow <= collect.sow <= opportunities[end].eow
        @test opportunities[end].sow <= collect.eow <= opportunities[end].eow
    end
end

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
    opportunities = find_all_opportunities(orb, images, sort=true)
    collects      = compute_collects_by_number(opportunities, 1)

    # Build graph of feasible transitions
    graph = sp_construct_graph(collects, Function[], 0)

    @test typeof(graph) == Dict{Collect, Array{Collect, 1}}
end