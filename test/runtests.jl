# Packages required for testing
using Test
using LinearAlgebra
using Logging
using Random
using SatelliteDynamics

# Package Under Test
using SatelliteTasking

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
    @time @testset "SatelliteTasking.Utils" begin
        include(joinpath(testdir, "test_utils.jl"))
    end
    @time @testset "SatelliteTasking.Visibility" begin
        include(joinpath(testdir, "test_visibility.jl"))
    end
end