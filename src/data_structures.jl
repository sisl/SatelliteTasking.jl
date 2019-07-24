# Exports
export Orbit
export interpolate
export Spacecraft
export Location
export Request
export load_requests
export set_collect_duration
export GroundStation
export load_stations
export Opportunity
export Done
export Noop
export Sunpoint
export Collect
export Contact
export PlanningProblem
export add_locations
export clear_locations
export clear_opportunities
export clear_all

#########
# Orbit #
#########

"""
Orbit provides a data structure which wraps SatelliteDynamics propagators

Arguments:
- `epc0::Union{Epoch, Real}` Initial Epoch of propagation
- `epcf::Union{Epoch, Real}` Final Epoch of propagation
- `eci0::Array{<:Real, 1}` Initial inertial state at initial Epoch [m; m/s]
- `id::Integer` Identification number for this orbit
- `timestep::Real` Timestep used for propgation output
- `dtmax::Real` Maximum allowable integrator step
- `mass::Real` Satellite mass [kg]
- `area_drag::Real` Wetted wind-facing area [m^2]
- `coef_drag::Real` Coefficient of drag
- `area_srp::Real` Area illuminated by solar radiation pressure [m^2]
- `coef_srp::Real` Coefficient of reflectivity
- `n_grav::Integer` Gravity model degree
- `m_grav::Integer` Gravity model order
- `drag::Bool` Enable drag perturbation in propgation
- `srp::Bool` Enable solar radiation pressure perturbation in propgation
- `moon::Bool` Enable thirdbody moon perturbation in propgation
- `sun::Bool` Enable thirdbody sun perturbation in propgation
- `relativity::Bool` Enable relativity perturbation in propgation

Attributes:
- `id::Int32` Identification number for the specific satellite
- `t::Array{Float64, 1}` Elapsed time of propagator output
- `epc::Array{Epoch, 1}` Epoch of propagator output
- `eci::Array{Float64, 2}` Earth Centered Inertial state of satellite output. State entries alightex with column output.
- `ecef::Array{Float64, 2}` Earth Centered Earth Fixed state of statellite output. State entries alightex with column output.
- `oe::Array{Float64, 2}` Osculating orbital element state of statellite output. State entries alightex with column output.
"""
mutable struct Orbit
    id::Int32
    t::Array{Float64, 1}
    epc::Array{Epoch, 1}
    eci::Array{Float64, 2}
    ecef::Array{Float64, 2}
    oe::Array{Float64, 2}

    function Orbit(epc0::Union{Epoch, Real},  epcf::Union{Epoch, Real}, eci0::Array{<:Real, 1}; 
                id=0::Integer, timestep=1.0::Real, dtmax=30.0::Real,
                mass=100.0::Real, area_drag=1.0::Real, coef_drag=2.3::Real, 
                area_srp=1.0::Real, coef_srp=1.8::Real, 
                n_grav=20::Integer, m_grav=20::Integer, 
                drag=true::Bool, srp=true::Bool, moon=true::Bool, 
                sun=true::Bool, relativity=true::Bool)
        
        # Initialize State Vector
        orb  = EarthInertialState(epc0, eci0, dt=timestep,
            mass=mass, 
            area_drag=area_drag, coef_drag=coef_drag, 
            area_srp=area_srp, coef_srp=coef_srp, 
            n_grav=n_grav, m_grav=m_grav, 
            drag=drag, srp=srp,
            moon=moon, sun=sun,
            relativity=relativity
        )

        # Propagate the orbit
        t, epc, eci = sim!(orb, epcf)

        ecef = zeros(Float64, 6, length(epc))
        oe   = zeros(Float64, 6, length(epc))

        for i in 1:length(t)
            ecef[:, i] = sECItoECEF(epc[i], eci[:, i])
            oe[:, i]   = sCARTtoOSC(eci[:, i], use_degrees=true)
        end

        new(id, t, epc, eci, ecef, oe)
    end
end

function Base.show(io::IO, orb::Orbit)

    s = @sprintf "Orbit(%d)" orb.id

    print(io, s)
end

# Helpfer function 
function searchsortednearest(a,x)
    idx = searchsortedfirst(a,x)
    if (idx==1); return idx; end
    if (idx>length(a)); return length(a); end
    if (a[idx]==x); return idx; end
    if (abs(a[idx]-x) < abs(a[idx-1]-x))
       return idx
    else
       return idx-1
    end
 end

"""
Interpolate orbit information to a given `Epoch`.

This interpolation is not particularly accurate, care should be used to determine
that the use-case can handle the interprolation errors introduced.

Arguments:
- `orbit::Orbit` Orbit data to interpolate
- `epc::Union{Epoch, Real}` Epoch to interpolate state output to.

Returns:
- `eci::Arrray{Float64, 1}` Earth intertial state information interpolated to time of Epoch
"""
function interpolate(orbit::Orbit, epc::Union{Epoch, Real})
    # Check validity of input
    if epc < orbit.epc[1]
        throw(ArgumentError("Invaid interpolation time. Cannot interpolate before start of Orbit."))
    end

    if epc > orbit.epc[end]
        throw(ArgumentError("Invaid interpolation time. Cannot interpolate after end of Orbit."))
    end

    # Exit early if epoch is end of interpolation (corner case)
    if epc == orbit.epc[end]
        return orbit.eci[:, end]
    end

    # Get Uppler and lower time bounds
    idx_l = searchsortedfirst(orbit.epc, epc)
    idx_u = idx_l + 1

    # Epoch
    epc_l = orbit.epc[idx_l]
    epc_u = orbit.epc[idx_u]

    # Earth Intertial State
    eci_l = orbit.eci[:, idx_l]
    eci_u = orbit.eci[:, idx_u]

    # Interpolate
    eci = zeros(typeof(eci_l[1]), 6)
    for i in 1:length(eci)
        eci[i] = (eci_u[i] - eci_l[i])/(epc_u - epc_l)*(epc - epc_l) + eci_l[i]
    end

    return eci
end

##############
# Spacecraft #
##############

"""
Object to store spacecraft information

Model Parameters:
- `powergen_sunpoint::Float64` Power generated when sunpointed
- `powergen_image::Float64` Power generated when imaging (should be negative)
- `powergen_downlink::Float64` Power generated when downlinking (should be negative)
- `datagen_backorbit::Float64` Data generated in backorbit
- `datagen_image::Float64` Data generated when imaging
- `datagen_downlink::Float64` Data generated when downlinking (should be negative)
- `slew_rate::Float64` Spacecraft slew rate (degrees per second)
"""
@with_kw mutable struct Spacecraft
    # Core Parameters
    id::Integer
    tle::TLE

    # Resource Model
    powergen_sunpoint::Float64 = 0.0
    powergen_image::Float64 = 0.0
    powergen_downlink::Float64 = 0.0
    datagen_backorbit::Float64 = 0.0
    datagen_image::Float64 = 0.0
    datagen_downlink::Float64 = 0.0
    slew_rate::Float64 = 1.0
end

############
# Location #
############

abstract type Location end

function Base.isequal(ll::Location, rl::Location)
    return (
        (ll.id == rl.id) &&
        (ll.lat == rl.lat) &&
        (ll.lon == rl.lon) &&
        (ll.ecef == rl.ecef) && true
    )
end
Base.:(==)(ll::Location, rl::Location) = Base.isequal(ll, rl)

function Base.show(io::IO, loc::Location)

    s = @sprintf "%s(ID: %s, Lon: %.3f, Lat: %0.3f)" typeof(loc) string(loc.id) loc.lon loc.lat

    print(io, s)
end

###########
# Request #
###########

"""
Request represents a single area/point on Earth to be observed. It represents a location
for which it is possible to collect the entire area in a singel collect

Attributes:
- `lon::Float64` Longitude of request center [deg]
- `lat::Float64` Lattitude of iamge center [deg]
- `ecef::Array{Float64, 1}` Earth-Centered-Earth-Fixed Coordinates of request center [m]
- `look_angle_min::Float64` Minimum look angle of valid collect
- `look_angle_max::Float64` Maximum look angle of valid collect
- `collect_duration` Required time to well of request for a feasible collection
- `require_zero_doppler::Bool` Require that the collect occurs centered at zero Doppler offset
- `reward::Float64` Reward for for request collection
- `id::Union{UUID, Integer}` Unique request identifier
"""
mutable struct Request <: Location
    lon::Float64
    lat::Float64
    ecef::Array{Float64, 1}
    look_angle_min::Float64
    look_angle_max::Float64
    require_zero_doppler::Bool
    collect_duration::Float64
    reward::Float64
    id::Union{UUID, Integer}

    function Request(lon::Real, lat::Real ; id::Union{UUID, Integer}=uuid4(),
                    look_angle_min::Real=5.0, look_angle_max::Real=55.0, 
                    require_zero_doppler::Bool=true, collect_duration::Real=1.0,
                    reward::Real=1.0, use_degrees::Bool=true)

        # Compute ECEF coordinates of scene center
        ecef = sGEODtoECEF([lon, lat], use_degrees=use_degrees::Real)

        # Generate new request object
        new(lon, lat, ecef, look_angle_min, look_angle_max,
            require_zero_doppler, collect_duration, reward, id)
    end
end

function Base.show(io::IO, req::Request)

    s = @sprintf "Request(ID: %s, Lon: %.3f, Lat: %0.3f, Reward: %.2f, Duration: %.2f)" string(req.id) req.lon req.lat req.reward req.collect_duration

    print(io, s)
end

"""
Loads request data from a JSON file.

THe JSON file is expected to contain an array of JSON objects where each object 
has the following format:

```json
{
    "lon": 22,
    "lat": 12.23,
    "reward": 42.0,
    "look_angle_min": 5.0,
    "look_angle_max": 55.0
}
```

Arguments:
- `file::String` Filepath to JSON file encoding request 
"""
function load_requests(file::String; id_offset::Integer=0, require_zero_doppler::Bool=true)
    # Parse data file
    data = JSON.parsefile(file)

    requests = Array{Request, 1}(undef, length(data))

    for (i, req) in enumerate(data)
        lon = req["lon"]
        lat = req["lat"]
        look_angle_min = req["look_angle_min"]
        look_angle_max = req["look_angle_max"]
        id  = i + id_offset

        requests[i] = Request(lon, lat, id=id,
            look_angle_min=look_angle_min,
            look_angle_max=look_angle_max
        )
    end

    println("Loaded $(length(requests)) requests.")

    return requests
end

"""
Set the request collection duration for satellites in the constellation.
"""
function set_collect_duration(requests::Array{<:Request, 1}, collect_duration::Real)
    if collect_duration <= 0.0
        throw(ArgumentError("Invalid collection duration."))
    end

    for request in requests
        request.collect_duration = collect_duration
    end
end

#################
# GroundStation #
#################

"""
GroundStation represents a location which the satellite can use to communicate
with the ground.

Attributes:
- `lon::Float64` Longitude of request center [deg]
- `lat::Float64` Lattitude of iamge center [deg]
- `ecef::Array{Float64, 1}` Earth-Centered-Earth-Fixed Coordinates of request center [m]
- `elevation_min::Float64` Minimum look angle of valid collect
- `id::Union{UUID, Integer}` Unique station identifier
"""
mutable struct GroundStation <: Location
    lon::Float64
    lat::Float64
    ecef::Array{Float64, 1}
    elevation_min::Float64
    duration_min::Float64
    id::Union{UUID, Integer}

    function GroundStation(lon::Real, lat::Real ; id::Union{UUID, Integer}=uuid4(),
                    elevation_min::Real=5.0, duration_min::Real=0.0,
                    use_degrees::Bool=true)

        # Compute ECEF coordinates of scene center
        ecef = sGEODtoECEF([lon, lat], use_degrees=use_degrees)

        # Generate new request object
        new(lon, lat, ecef, elevation_min, duration_min, id)
    end
end

# function Base.isequal(ll::GroundStation, rl::GroundStation)
#     return (
#         (ll.id == rl.id) &&
#         (ll.lat == rl.lat) &&
#         (ll.lon == rl.lon) &&
#         (ll.ecef == rl.ecef) && true
#     )
# end
# Base.:(==)(ll::GroundStation, rl::GroundStation) = Base.isequal(ll, rl)

function Base.show(io::IO, sta::GroundStation)

    s = @sprintf "GroundStation(ID: %s, Lon: %.3f, Lat: %0.3f)" string(sta.id) sta.lon sta.lat

    print(io, s)
end


"""
Loads request data from a JSON file.

THe JSON file is expected to contain an array of JSON objects where each object 
has the following format:

```json
{
    "lon": 22,
    "lat": 12.23,
    "elevation_min": 5.0
}
```

Arguments:
- `file::String` Filepath to JSON file encoding request 
"""
function load_stations(file::String, id_offset::Integer=0)
    data = JSON.parsefile(file)

    # Initialize array of Requests
    stations = Array{GroundStation, 1}(undef, length(data))

    for (i, sta) in enumerate(data)
        lon = sta["lon"]
        lat = sta["lat"]
        elevation_min = sta["elevation_min"]
        id  = i + id_offset


        stations[i] = GroundStation(lon, lat, id=id, elevation_min=elevation_min)
    end

    println("Loaded $(length(stations)) stations.")

    return stations
end

"""
Set minimum access duration for all stations
"""
function set_station_duration_min(stations::Array{<:GroundStation, 1}, duration_min::Real)
    if duration_min <= 0.0
        throw(ArgumentError("Invalid minimum contact duration duration."))
    end

    for station in stations
        station.duration_min = duration_min
    end
end

###############
# Opportunity #
###############

"""
Opportunity represents an "opportunity" for action for a spacecraft

Common Attributes:
- `id::Integer` Unique Identifier for action
- `orbit::TLE` Orbit object associated with this collect
- `location::Location` Location object on the ground for interaction
- `t_start::Union{Epoch, Real}` Start of acquisition window 
- `t_mid::Union{Epoch, Real}` Mid-time of acquisition window
- `t_end::Union{Epoch, Real}` End of possible acquisition window
- `duration::Real` Length of action
"""
abstract type Opportunity end

"""
"""
@with_kw mutable struct Done <: Opportunity
    t_start::Union{Epoch, Real} = Inf
end

"""
No-operation action. Advances state without changing other values
"""
@with_kw mutable struct Noop <: Opportunity
    t_start::Union{Epoch, Real} = 0.0
end

function Base.show(io::IO, noop::Noop)

    s = @sprintf "Noop(time: %s)" noop.t_start

    print(io, s)
end

@with_kw mutable struct Sunpoint <: Opportunity
    t_start::Union{Epoch, Real} = 0.0
end

function Base.show(io::IO, sp::Sunpoint)

    s = @sprintf "Sunpoint(time: %s)" sp.t_start

    print(io, s)
end

"""
Opportunity type for a Request collection.
"""
@with_kw mutable struct Collect <: Opportunity
    id::Integer
    spacecraft::Spacecraft
    location::Union{Request, Nothing}
    t_start::Union{Epoch, Real}
    t_mid::Union{Epoch, Real}
    t_end::Union{Epoch, Real}
    duration::Float64
    reward::Float64
    nr::Int

    function Collect(t_start::Union{Epoch, Real}, t_end::Union{Epoch, Real};
        id::Integer=0, spacecraft=nothing::Union{Spacecraft, Nothing}, 
        location=nothing::Union{Location, Nothing},
        nr::Int=0)

        duration = t_end - t_start
        t_mid    = t_start + duration/2.0
        new(id, spacecraft, location, t_start, t_mid, t_end, duration, location.reward, nr)
    end
end

function JSON.lower(col::Collect)
    return Dict(
        "id" => col.id,
        "location" => col.location.id,
        "lon" => col.location.lon,
        "lat" => col.location.lat,
        "t_start" => string(col.t_start),
        "t_end"  => string(col.t_end),
        "duration" => col.duration
    )    
end

function Base.isequal(cl::Collect, cr::Collect)
    return (
        (cl.id == cr.id) &&
        ((cl.spacecraft == nothing && cr.spacecraft == nothing) ||
         ((cl.spacecraft != nothing && cr.spacecraft != nothing &&
          (cl.spacecraft.id == cr.spacecraft.id)))) &&
        ((cl.location == nothing && cr.location == nothing) ||
         ((cl.location != nothing && cr.location != nothing &&
          (cl.location.id == cr.location.id)))) &&
        (cl.t_start == cr.t_start) &&
        (cl.t_end == cr.t_end) && true
    )
end
Base.:(==)(cl::Collect, cr::Collect) = Base.isequal(cl, cr)

function Base.show(io::IO, col::Collect)

    if typeof(col.t_start) == Epoch
        s = @sprintf "Collect(ID: %s, Spacecraft %s, Request: %s, t_start: %s, t_end: %s, Reward: %.2f, Remaining: %d)" string(col.id) string(col.spacecraft.id) string(col.location.id) col.t_start col.t_end col.reward col.nr
    else
        s = @sprintf "Collect(ID: %s, Spacecraft %s, Request: %s, t_start: %.3f, t_end: %.3f, Reward: %.2f, Remaining: %d)" string(col.id) string(col.spacecraft.id) string(col.location.id) col.t_start col.t_end col.reward col.nr
    end

    print(io, s)
end

"""
Opportunity type for a ground contact.
"""
mutable struct Contact <: Opportunity
    id::Integer
    spacecraft::Spacecraft
    location::Union{GroundStation, Nothing}
    t_start::Union{Epoch, Real}
    t_mid::Union{Epoch, Real}
    t_end::Union{Epoch, Real}
    duration::Float64
    nr::Int

    function Contact(t_start::Union{Epoch, Real}, t_end::Union{Epoch, Real};
        id::Integer=0, spacecraft=nothing::Union{Spacecraft, Nothing}, 
        location=nothing::Union{Location, Nothing}, nr::Int=0)

        duration = t_end - t_start
        t_mid    = t_start + duration/2.0
        new(id, spacecraft, location, t_start, t_mid, t_end, duration, nr)
    end
end

function JSON.lower(con::Contact)
    return Dict(
        "id" => con.id,
        "location" => con.location.id,
        "lon" => con.location.lon,
        "lat" => con.location.lat,
        "t_start" => string(con.t_start),
        "t_end"  => string(con.t_end),
        "duration" => con.duration
    )    
end

function Base.isequal(cl::Contact, cr::Contact)
    return (
        (cl.id == cr.id) &&
        ((cl.spacecraft == nothing && cr.spacecraft == nothing) ||
         ((cl.spacecraft != nothing && cr.spacecraft != nothing &&
          (cl.spacecraft.id == cr.spacecraft.id)))) &&
        ((cl.location == nothing && cr.location == nothing) ||
         ((cl.location != nothing && cr.location != nothing &&
          (cl.location.id == cr.location.id)))) &&
        (cl.t_start == cr.t_start) &&
        (cl.t_end == cr.t_end) && true
    )
end
Base.:(==)(cl::Contact, cr::Contact) = Base.isequal(cl, cr)

function Base.show(io::IO, con::Contact)

    minutes = floor(con.duration/60)
    seconds = con.duration - 60*minutes

    if typeof(con.t_start) == Epoch
        s = @sprintf "Contact(ID: %s,  Spacecraft %s, Station: %s, t_start: %s, t_end: %s, Duration: %02dm %02ds, Remaining: %d)" string(con.id) string(con.spacecraft.id) string(con.location.id) con.t_start con.t_end minutes seconds con.nr
    else
        s = @sprintf "Contact(ID: %s,  Spacecraft %s, Station: %s, t_start: %.3f, t_end: %.3f, Duration: %02dm %02ds, Remaining: %d)" string(con.id) string(con.spacecraft.id) string(con.location.id) con.t_start con.t_end minutes seconds con.nr
    end

    print(io, s)
end

####################
# Planning Problem #
####################

@with_kw mutable struct PlanningProblem
    # General Planning Settings
    t_start::Epoch
    t_end::Epoch

    # Spacecraft Settings
    spacecraft::Array{Spacecraft, 1} = Spacecraft[]
    constraints::Array{Function,1} = Function[]

    # Planning Locations
    locations::Array{<:Location, 1} = Location[]
    stations::Array{GroundStation, 1} = GroundStation[]
    requests::Array{Request, 1} = Request[]

    # Planning Opportunities
    opportunities::Array{<:Opportunity, 1} = Opportunity[]
    contacts::Array{Contact, 1} = Contact[]
    collects::Array{Collect, 1} = Collect[]

    # Lookup Tables
    lt_locations     = Dict{Union{Integer, UUID}, Location}()
    lt_opportunities = Dict{Union{Integer, UUID}, Opportunity}()
    lt_contacts      = Dict{Union{Integer, UUID}, Contact}()
    lt_collects      = Dict{Union{Integer, UUID}, Collect}()
    lt_loc_opps = Dict{Union{Integer, UUID}, Array{Union{Integer, UUID}, 1}}()

    # Solver Parameters - General 
    solve_gamma::Real   = 1.0
    solve_depth::Int    = 3
    solve_breadth::Int  = 3
    solve_horizon::Real = 90.0

    # Solver Parameters - MCTS
    mcts_rollout_iterations::Int = 10
    mcts_c::Real = 1.0
end

"""
Clear problem locations.
"""
function clear_locations(problem::PlanningProblem)
    # Reset Locations
    problem.locations = Location[]
    problem.stations = GroundStation[]
    problem.requests = Request[]
end

"""
Clear problem opportunities.
"""
function clear_opportunities(problem::PlanningProblem)
    # Reset Opportunities
    problem.opportunities = Opportunity[]
    problem.contacts = Contact[]
    problem.collects = Collect[]

    problem.lt_opportunities = Dict{Union{Integer, UUID}, Opportunity}()
    problem.lt_contacts      = Dict{Union{Integer, UUID}, Contact}()
    problem.lt_collects      = Dict{Union{Integer, UUID}, Collect}()
    problem.lt_loc_opps      = Dict{Union{Integer, UUID}, Array{Union{Integer, UUID}, 1}}()
end

"""
Clear all problem opportunities
"""
function clear_all(problem::PlanningProblem)
    clear_locations(problem)
    clear_opportunities(problem)
end

"""
Add locations to planning problem.
"""
function add_locations(problem::PlanningProblem, locations::Array{<:Location})
    if length(locations) == 0
        @warn "No locations in array. Nothing to add..."
        return
    end

    # TODO: Check for duplicate IDs in requests (assumes already unique)

    # Add location to problem
    push!(problem.locations, locations...)

    # Create separate arrays of contacts and requests
    push!(problem.stations, filter(x -> typeof(x) == GroundStation, locations)...)
    push!(problem.requests, filter(x -> typeof(x) == Request, locations)...)

    # Update location lookup table
    for loc in problem.locations
        problem.lt_locations[loc.id] = loc
    end

    return
end