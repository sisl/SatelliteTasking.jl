__precompile__(true)
module DataStructures # Satellite Tasking Data Structures

using SatelliteDynamics.Time
using SatelliteDynamics.Astrodynamics
using SatelliteDynamics.ReferenceSystems
using SatelliteDynamics.Coordinates
using SatelliteDynamics.Simulation
using SatelliteDynamics.SGPModels

using LinearAlgebra
using Statistics
using JSON
using UUIDs
using Printf

# Exports
export Orbit
export interpolate
export Location
export Request
export load_requests
export set_collect_duration
export GroundStation
export load_stations
export Opportunity
export Collect
export Contact

#########
# Orbit #
#########

"""
Orbit provides a data structure which wraps SatelliteDynamics propagators

Arguments:
- `epc0::Epoch` Initial Epoch of propagation
- `epcf::Epoch` Final Epoch of propagation
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

    function Orbit(epc0::Epoch,  epcf::Epoch, eci0::Array{<:Real, 1}; 
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
- `epc::Epoch` Epoch to interpolate state output to.

Returns:
- `eci::Arrray{Float64, 1}` Earth intertial state information interpolated to time of Epoch
"""
function interpolate(orbit::Orbit, epc::Epoch)
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
function load_requests(file::String; require_zero_doppler::Bool=true)
    # Parse data file
    data = JSON.parsefile(file)

    requests = Array{Request, 1}(undef, length(data))

    for (i, req) in enumerate(data)
        lon = req["lon"]
        lat = req["lat"]
        look_angle_min = req["look_angle_min"]
        look_angle_max = req["look_angle_max"]
        id  = i

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
function load_stations(file::String)
    data = JSON.parsefile(file)

    # Initialize array of Requests
    stations = Array{GroundStation, 1}(undef, length(data))

    for (i, sta) in enumerate(data)
        lon = sta["lon"]
        lat = sta["lat"]
        elevation_min = sta["elevation_min"]
        id  = i


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
- `t_start::Epoch` Start of acquisition window 
- `t_mid::Epoch` Mid-time of acquisition window
- `t_end::Epoch` End of possible acquisition window
- `duration::Real` Length of action
"""
abstract type Opportunity end
"""
Opportunity type for a Request collection.
"""
mutable struct Collect <: Opportunity
    id::Integer
    orbit::Union{TLE, Nothing}
    location::Union{Request, Nothing}
    t_start::Epoch
    t_mid::Epoch
    t_end::Epoch
    duration::Float64
    reward::Float64

    function Collect(t_start::Epoch, t_end::Epoch;
        id::Integer=0, orbit=nothing::Union{Orbit, TLE, Nothing}, 
        location=nothing::Union{Location, Nothing})

        duration = t_end - t_start
        t_mid    = t_start + duration/2.0
        new(id, orbit, location, t_start, t_mid, t_end, duration, location.reward)
    end
end

Base.:(==)(ol::Collect, or::Collect) = Base.isequal(ol, or)

function Base.show(io::IO, col::Collect)

    s = @sprintf "Collect(ID: %s, Request: %s, t_start: %s, t_end: %s, Reward: %.2f)" string(col.id) string(col.location.id) col.t_start col.t_end col.reward

    print(io, s)
end

"""
Opportunity type for a ground contact.
"""
mutable struct Contact <: Opportunity
    id::Integer
    orbit::Union{TLE, Nothing}
    location::Union{GroundStation, Nothing}
    t_start::Epoch
    t_mid::Epoch
    t_end::Epoch
    duration::Float64

    function Contact(t_start::Epoch, t_end::Epoch;
        id::Integer=0, orbit=nothing::Union{Orbit, TLE, Nothing}, 
        location=nothing::Union{Location, Nothing})

        duration = t_end - t_start
        t_mid    = t_start + duration/2.0
        new(id, orbit, location, t_start, t_mid, t_end, duration)
    end
end

function Base.show(io::IO, con::Contact)

    minutes = floor(con.duration/60)
    seconds = con.duration - 60*minutes

    s = @sprintf "Contact(ID: %s, Station: %s, t_start: %s, t_end: %s, Duration: %02dm %02ds)" string(con.id) string(con.location.id) con.t_start con.t_end minutes seconds

    print(io, s)
end

end # Close DataStructures