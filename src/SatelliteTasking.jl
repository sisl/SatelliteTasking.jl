__precompile__(true)
module SatelliteTasking

using Reexport

# Base Module Includes
include("data_structures.jl")
include("collection.jl")
include("simulation.jl")

# Export Values
@reexport using SatelliteTasking.DataStructures
@reexport using SatelliteTasking.Collection
@reexport using SatelliteTasking.Simulation

# Export satellite planning submodule
include(joinpath(".", "satellite_planning", "satellite_planning.jl"))

# Export analysis submodule
# println(joinpath("analysis", "analysis.jl"))
include(joinpath(".", "analysis", "analysis.jl"))

end # module
