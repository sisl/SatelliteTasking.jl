# Exports
export Spacecraft, MDPProblem
export Geolocation, Request, Station
export Opportunity, Collect, Contact

"""
Top level type type for geolocations 
"""
abstract type Geolocation end

"""
Tasking request definition
"""
struct Request <: Geolocation
    # General Parameters
    id::Int
    lon::Real
    lat::Real
    alt::Real
    ecef::Array{<:Real, 1}

    # Request specific
    reward::Real
    look_angle_max::Real
    duration::Real
end

function Request(;id::Int=0, lon::Real=0.0, lat::Real=0.0, alt::Real=0.0,
    reward::Real=0.0, look_angle_max::Real=0.0, duration::Real=0.0)

    # Compute ECEF coordinates of scene center
    ecef = sGEODtoECEF([lon, lat, alt], use_degrees=true)

    # Generate new image object
    Request(id, lon, lat, alt, ecef, reward, look_angle_max, duration)
end

"""
Ground station definition
"""
struct Station <: Geolocation
    # General Parameters
    id::Int
    lon::Real
    lat::Real
    alt::Real
    ecef::Array{<:Real, 1}

    # Request specific
    elevation_min::Real
    name::String
end

function Station(;id::Int=0, lon::Real=0.0, lat::Real=0.0, alt::Real=0.0,
    elevation_min::Real=0.0, name::String="undefined")

    # Compute ECEF coordinates of scene center
    ecef = sGEODtoECEF([lon, lat, alt], use_degrees=true)

    # Generate new image object
    Station(id, lon, lat, alt, ecef, elevation_min, name)
end

"""
Action opportunity
"""
abstract type Opportunity end

"""
Image collection opportunity
"""
@with_kw struct Collect <: Opportunity
    # Opporttunity collects
    id::Int = 0
    geolocation_id::Int = 0
    duration::Real = 0.0
    window_open::Real = 0.0
    window_close::Real = 0.0

    # Collect specific
    reward::Real = 0.0
end

"""
Ground station contact definition
"""
@with_kw struct Contact <: Opportunity
    # Opporttunity collects
    id::Int = 0
    geolocation_id::Int = 0
    duration::Real = 0.0
    window_open::Real = 0.0
    window_close::Real = 0.0
end


"""
Spacecraft model containing parameter settings and values
"""
@with_kw struct Spacecraft
    power_generation::Real = 0.0
    power_draw_downlink::Real = 0.0
    power_draw_imaging::Real = 0.0
    data_generation_backorbit::Real = 0.0
    data_generation_imaging::Real = 0.0
    data_downlink::Real = 0.0
end

"""
Define a Markov Decision Proccess Problem 
"""
@with_kw struct MDPProblem
    # Problem Parameters
    t_start::Epoch = Epoch(2019, 1, 1, 0, 0, 0)
    t_end::Epoch   = Epoch(2019, 1, 2, 0, 0, 0)

    # Constraint List
    constraint_list::Array{Function,1} = Function[]

    # Spacecraft Parameters
    spacecraft::Spacecraft   = Spacecraft()
    tle::Union{TLE, Nothing} = nothing
    locations::Array{<:Geolocation, 1} = Geolocation[]
    opportunities::Array{<:Opportunity, 1} = Opportunity[]

    # Solver Parameters - General 
    gamma::Real  = 0.0
    depth::Int   = 10
    breadth::Int = 0

    # Solver Parameters - MCTS
    alpha::Real = 1.0
    max_iterations::Int = 10
    c::Real = 1.0
end