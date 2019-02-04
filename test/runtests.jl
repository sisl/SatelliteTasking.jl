# Packages required for testing
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

# Set logging level
global_logger(SimpleLogger(stderr, Logging.Debug))

# Fix randomness during tests
Random.seed!(0)

# Check equality of two arrays
@inline function array_isapprox(x::AbstractArray{F},
                  y::AbstractArray{F};
                  rtol::F=sqrt(eps(F)),
                  atol::F=zero(F)) where {F<:AbstractFloat}

    # Easy check on matching size
    if length(x) != length(y)
        return false
    end

    for (a,b) in zip(x,y)
        @test isapprox(a,b, rtol=rtol, atol=atol)
    end
end

# Check if array equals a single value
@inline function array_isapprox(x::AbstractArray{F},
                  y::F;
                  rtol::F=sqrt(eps(F)),
                  atol::F=zero(F)) where {F<:AbstractFloat}

    for a in x
        @test isapprox(a, y, rtol=rtol, atol=atol)
    end
end

@inline function array_isapprox(x::AbstractArray{F}, y::AbstractArray{F}) where {F<:Integer}
    # Easy check on matching size
    if length(x) != length(y)
        return false
    end

    for (a,b) in zip(x,y)
        @test a == b
    end
end

# Check if array equals a single value
@inline function array_isapprox(x::AbstractArray{F}, y::F) where {F<:Integer}
    for a in x
        @test a == y
    end
end

@time @testset "SatelliteTasking Package Tests" begin
    testdir = joinpath(dirname(@__DIR__), "test")
    # Main Package Test Set
    @time @testset "SatelliteTasking.DataStructures" begin
        include(joinpath(testdir, "test_data_structures.jl"))
    end
    @time @testset "SatelliteTasking.Collection" begin
        include(joinpath(testdir, "test_collection.jl"))
    end
    @time @testset "SatelliteTasking.Constraints" begin
        include(joinpath(testdir, "test_constraints.jl"))
    end
    @time @testset "SatelliteTasking.Simulation" begin
        include(joinpath(testdir, "test_simulation.jl"))
    end

    # Satellite Planning Test Set
    @time @testset "SatelliteTasking.SatellitePlanning.Graph" begin
        include(joinpath(testdir, "satellite_planning", "test_graph.jl"))
    end
    @time @testset "SatelliteTasking.SatellitePlanning.MILP" begin
        include(joinpath(testdir, "satellite_planning", "test_milp.jl"))
    end
    @time @testset "SatelliteTasking.SatellitePlanning.MDP" begin
        include(joinpath(testdir, "satellite_planning", "test_mdp.jl"))
    end

    # Analysis Test Set
    @time @testset "SatelliteTasking.Analysis.CollectUncertainty" begin
        include(joinpath(testdir, "analysis", "test_collect_uncertainty.jl"))
    end
end