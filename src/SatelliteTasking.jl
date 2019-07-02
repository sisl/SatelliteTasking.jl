__precompile__(true)
module SatelliteTasking

using Reexport


# Base Module Includes
include("data_structures.jl")
include("access.jl")
include("models.jl")

# Export Values
@reexport using SatelliteTasking.DataStructures
@reexport using SatelliteTasking.Access
@reexport using SatelliteTasking.Models
# @reexport using SatelliteTasking.Simulation

# Direct includes
include("visualization.jl")


# Export satellite planning submodule
# include(joinpath(".", "satellite_planning", "satellite_planning.jl"))

# Export analysis submodule
# include(joinpath(".", "analysis", "analysis.jl"))
# @reexport using SatelliteTasking.Analysis

end # module