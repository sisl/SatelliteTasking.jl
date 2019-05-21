__precompile__(true)
module DataStructures # Satellite Tasking Data Structures

using SatelliteDynamics.Time
using SatelliteDynamics.Astrodynamics
using SatelliteDynamics.ReferenceSystems
using SatelliteDynamics.Coordinates
using SatelliteDynamics.Simulation

using LinearAlgebra
using Statistics
using JSON
using UUIDs
using Printf

#########
# Orbit #
#########

export Orbit
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

export interpolate
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

##########
# Location #
##########

abstract type Location end

#########
# Image #
#########

export Image
"""
Image represents a single area/point on Earth to be observed. It represents a location
for which it is possible to collect the entire area in a singel collect

Attributes:
- `lon::Float64` Longitude of image center [deg]
- `lat::Float64` Lattitude of iamge center [deg]
- `ecef::Array{Float64, 1}` Earth-Centered-Earth-Fixed Coordinates of image center [m]
- `look_angle_min::Float64` Minimum look angle of valid collect
- `look_angle_max::Float64` Maximum look angle of valid collect
- `collect_duration` Required time to well of image for a feasible collection
- `require_zero_doppler::Bool` Require that the collect occurs centered at zero Doppler offset
- `reward::Float64` Reward for for image collection
- `id::UUID` Unique image identifier
"""
mutable struct Image <: Location
    lon::Float64
    lat::Float64
    ecef::Array{Float64, 1}
    look_angle_min::Float64
    look_angle_max::Float64
    require_zero_doppler::Bool
    collect_duration::Float64
    reward::Float64
    id::UUID

    function Image(lon::Real, lat::Real ; id::UUID=uuid4(),
                    look_angle_min::Real=5.0, look_angle_max::Real=55.0, 
                    require_zero_doppler::Bool=false, collect_duration::Real=1.0,
                    reward::Real=1.0, use_degrees::Bool=true)

        # Compute ECEF coordinates of scene center
        ecef = sGEODtoECEF([lon, lat], use_degrees=use_degrees::Real)

        # Generate new image object
        new(lon, lat, ecef, look_angle_min, look_angle_max,
            require_zero_doppler, collect_duration, reward, id)
    end
end

function Base.isequal(ll::Image, rl::Image)
    return (
        (ll.id == rl.id) &&
        (ll.lat == rl.lat) &&
        (ll.lon == rl.lon) &&
        (ll.ecef == rl.ecef) && true
    )
end
Base.:(==)(ll::Image, rl::Image) = Base.isequal(ll, rl)

function Base.show(io::IO, img::Image)

    s = @sprintf "Image(Ptr: %s Lon: %.3f Lat: %0.3f Reward: %f)" string(UInt64(pointer_from_objref(img)),base=16) img.lon img.lat img.reward

    print(io, s)
end

export load_images
"""
Loads image data from a JSON file.

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
- `file::String` Filepath to JSON file encoding image 
"""
function load_images(file::String; collect_duration::Real=1.0)
    data = JSON.parsefile(file)

    # Initialize array of Images
    n_img  = length(data["images"])
    images = Array{Image, 1}(undef, n_img)

    for (i, img) in enumerate(data["images"])
        lon = img["lon"]
        lat = img["lat"]
        look_angle_min = img["look_angle_min"]
        look_angle_max = img["look_angle_max"]
        reward = img["reward"]
        id  = uuid4()

        # Read from file if present
        img_keys = keys(data)

        if "id" in img_keys
            id = UUID(img["id"])
        end

        if "collect_duration" in img_keys
            collect_duration = img["collect_duration"]
        end

        images[i] = Image(lon, lat, id=id, 
                        look_angle_min=look_angle_min,
                        look_angle_max=look_angle_max, 
                        collect_duration=collect_duration,
                        require_zero_doppler=false, reward=reward)
    end

    return images
end

#################
# GroundStation #
#################

export GroundStation
"""
GroundStation represents a location which the satellite can use to communicate
with the ground.

Attributes:
- `lon::Float64` Longitude of image center [deg]
- `lat::Float64` Lattitude of iamge center [deg]
- `ecef::Array{Float64, 1}` Earth-Centered-Earth-Fixed Coordinates of image center [m]
- `elevation_min::Float64` Minimum look angle of valid collect
- `id::UUID` Unique station identifier
"""
mutable struct GroundStation <: Location
    lon::Float64
    lat::Float64
    ecef::Array{Float64, 1}
    elevation_min::Float64
    id::UUID

    function GroundStation(lon::Real, lat::Real ; id::UUID=uuid4(),
                    elevation_min::Real=5.0, use_degrees::Bool=true)

        # Compute ECEF coordinates of scene center
        ecef = sGEODtoECEF([lon, lat], use_degrees=use_degrees)

        # Generate new image object
        new(lon, lat, ecef, elevation_min, id)
    end
end

function Base.isequal(ll::GroundStation, rl::GroundStation)
    return (
        (ll.id == rl.id) &&
        (ll.lat == rl.lat) &&
        (ll.lon == rl.lon) &&
        (ll.ecef == rl.ecef) && true
    )
end
Base.:(==)(ll::GroundStation, rl::GroundStation) = Base.isequal(ll, rl)

function Base.show(io::IO, img::GroundStation)

    s = @sprintf "GroundStation(Ptr: %s Lon: %.3f Lat: %0.3f)" string(UInt64(pointer_from_objref(img)),base=16) img.lon img.lat 

    print(io, s)
end


export load_stations
"""
Loads image data from a JSON file.

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
- `file::String` Filepath to JSON file encoding image 
"""
function load_stations(file::String)
    data = JSON.parsefile(file)

    # Initialize array of Images
    n_stations = length(data)
    stations   = Array{GroundStation, 1}(undef, n_stations)

    for (i, sta) in enumerate(data)
        lon = sta["lon"]
        lat = sta["lat"]
        elevation_min = sta["elevation_min"]
        id  = uuid4()

        # Read from file if present
        sta_keys = keys(data)

        if "id" in sta_keys
            id = UUID(sta["id"])
        end

        stations[i] = GroundStation(lon, lat, id=id, elevation_min=elevation_min)
    end

    return stations
end

###############
# Opportunity #
###############

export Opportunity
"""
Opportunity represents a single 

Attributes:
- `id::UUID` Unique collect identifier
- `orbit::Orbit` Orbit object associated with this collect
- `location::Location` Location object on the ground for interaction
- `collect_duration` Required dewell time to make collection
- `sow::Epoch` Start of acquisition window 
- `mid::Epoch` Mid-time of acquisition window
- `eow::Epoch` End of possible acquisition window
- `duration::` Length of acquisition window
"""
mutable struct Opportunity
    id::UUID
    orbit::Union{Orbit, Nothing}
    location::Union{Location, Nothing}
    sow::Epoch
    mid::Epoch
    eow::Epoch
    duration::Float64
    collect_duration::Float64

    function Opportunity(sow::Epoch, eow::Epoch;
                         id=uuid4()::UUID, orbit=nothing::Union{Orbit, Nothing}, 
                         location=nothing::Union{Location, Nothing},
                         collect_duration=0::Real)

        duration = eow - sow
        mid      = sow + duration/2.0
        new(id, orbit, location, sow, mid, eow, duration, collect_duration)
    end
end

function Base.isequal(ol::Opportunity, or::Opportunity)
    return (
        (ol.id == or.id) &&
        (ol.location == or.location) &&
        (ol.sow == or.sow) &&
        (ol.mid == or.mid) &&
        (ol.eow == or.eow) &&
        (ol.duration == or.duration) &&
        (ol.collect_duration == or.collect_duration) && true
    )
end
Base.:(==)(ol::Opportunity, or::Opportunity) = Base.isequal(ol, or)

function Base.show(io::IO, opp::Opportunity)

    orbit = "nothing"
    if opp.orbit != nothing 
        orbit = string(opp.orbit.id)
    end

    image = "nothing"
    if opp.location != nothing 
        image = string(UInt64(pointer_from_objref(opp.location)), base=16)
    end
    
    s = @sprintf "Opportunity(Ptr: %s, a Orbit: %s, Location: %s, Start: %s, End: %s, Duration: %.2f)" string(UInt64(pointer_from_objref(opp)), base=16) orbit image string(opp.sow) string(opp.eow) opp.duration

    print(io, s)
end

end # Close DataStructures