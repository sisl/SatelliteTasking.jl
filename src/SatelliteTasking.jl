__precompile__(true)
module SatelliteTasking

# Julia Packages
using Printf
using JSON
using UUIDs
using Distributed
using LinearAlgebra
using Statistics
using Distributions
using Discretizers
using Parameters

# Package imports
using JuMP
using Gurobi

# Other dependencies
using SatelliteDynamics

# Includes
include("data_structures.jl")
include("models.jl")
include("utils.jl")
include("access.jl")
include("graph.jl")
include("milp.jl")
include("mdpcore.jl")
include("mdpfs.jl")
include("mdpmcts.jl")
include("analysis.jl")
include("visualization.jl")


# Export satellite planning submodule
# include(joinpath(".", "satellite_planning", "satellite_planning.jl"))

# Export analysis submodule
# include(joinpath(".", "analysis", "analysis.jl"))
# @reexport using SatelliteTasking.Analysis

end # End module