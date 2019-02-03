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
struct Orbit
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
        
            # Initialze orbit and simulate trajectory
        t, epc, eci = propagate_orbit(epc0, eci0, epcf, timestep=timestep, dtmax=max(dtmax, timestep),
                                    mass=mass, area_drag=area_drag, coef_drag=coef_drag, 
                                    area_srp=area_srp, coef_srp=coef_srp, 
                                    n_grav=n_grav, m_grav=m_grav, 
                                    drag=drag, srp=srp, moon=moon, sun=sun, relativity=relativity)

        ecef = zeros(Float64, 6, length(epc))
        oe   = zeros(Float64, 6, length(epc))

        for i in 1:length(t)
            ecef[:, i] = sECItoECEF(epc[i], eci[:, i])
            oe[:, i]   = sCARTtoOSC(eci[:, i], use_degrees=true)
        end

        new(id, t, epc, eci, ecef, oe)
    end
end

#########
# Image #
#########

export Image
"""
Image represents a single area/point on Earth to be observed. It represents a target
for which it is possible to collect the entire area in a singel collect

Attributes:
- `lon::Float64` Longitude of image center [deg]
- `lat::Float64` Lattitude of iamge center [deg]
- `ecef::Array{Float64, 1}` Earth-Centered-Earth-Fixed Coordinates of image center [m]
- `look_angle_min::Float64` Minimum look angle of valid collect
- `look_angle_max::Float64` Maximum look angle of valid collect
- `dwell_time` Required time to well of image for a feasible collection
- `require_zero_doppler::Bool` Require that the collect occurs centered at zero Doppler offset
- `reward::Float64` Reward for for image collection
- `id::UUID` Unique image identifier
"""
struct Image
    lon::Float64
    lat::Float64
    ecef::Array{Float64, 1}
    look_angle_min::Float64
    look_angle_max::Float64
    require_zero_doppler::Bool
    dwell_time::Float64
    reward::Float64
    id::UUID

    function Image(lon::Real, lat::Real ; id::UUID=uuid4(),
                    look_angle_min=5.0::Real, look_angle_max=55.0::Real, 
                    require_zero_doppler=false::Bool, dwell_time=1.0::Real,
                    reward=1.0::Real, use_degrees=true::Bool)

        # Compute ECEF coordinates of scene center
        ecef = sGEODtoECEF([lon, lat], use_degrees=use_degrees::Real)

        # Generate new image object
        new(lon, lat, ecef, look_angle_min, look_angle_max,
            require_zero_doppler, dwell_time, reward, id)
    end
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
function load_images(file::String; dwell_time=1.0::Real)
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

        if "dwell_time" in img_keys
            dwell_time = img["dwell_time"]
        end

        images[i] = Image(lon, lat, id=id, 
                        look_angle_min=look_angle_min,
                        look_angle_max=look_angle_max, 
                        dwell_time=dwell_time,
                        require_zero_doppler=false, reward=reward)
    end

    return images
end

###########
# Opportunity #
###########

export Opportunity
"""
Opportunity represents a single 

Attributes:
- `id::UUID` Unique collect identifier
- `orbit::Orbit` Orbit object associated with this collect
- `image::Image` Image object associated with this collect
- `dwell_time` Required dewell time to make collection
- `sow::Epoch` Start of acquisition window 
- `mid::Epoch` Mid-time of acquisition window
- `eow::Epoch` End of possible acquisition window
- `duration::` Length of acquisition window
"""
mutable struct Opportunity
    id::UUID
    orbit::Union{Orbit, Nothing}
    image::Union{Image, Nothing}
    sow::Epoch
    mid::Epoch
    eow::Epoch
    duration::Float64
    dwell_time::Float64

    function Opportunity(sow::Epoch, eow::Epoch;
                         id=uuid4()::UUID, orbit=nothing::Union{Orbit, Nothing}, 
                         image=nothing::Union{Image, Nothing},
                         dwell_time=0::Real)

        duration = eow - sow
        mid      = sow + duration/2.0
        new(id, orbit, image, sow, mid, eow, duration, dwell_time)
    end
end

###########
# Collect #
###########

export Collect
"""
Image represents a single 

Attributes:
- `id::UUID` Unique collect identifier
- `orbit::Orbit` Orbit object associated with this collect
- `image::Image` Image object associated with this collect
- `opportunity::Opportunity` Opportunity associated with this collect
- `start::Epoch` Start of acquisition window 
- `end::Epoch` End of possible acquisition window
- `duration::Float64` Length of acquisition window
"""
mutable struct Collect
    id::UUID
    orbit::Union{Orbit, Nothing}
    image::Union{Image, Nothing}
    opportunity::Union{Opportunity, Nothing}
    sow::Epoch
    mid::Epoch
    eow::Epoch
    duration::Float64

    function Collect(sow::Epoch, eow::Epoch;
                         id=uuid4()::UUID, orbit=nothing::Union{Orbit, Nothing}, 
                         image=nothing::Union{Image, Nothing}, opportunity=nothing::Union{Opportunity, Nothing},
                         dwell_time=0::Real)

        duration = eow - sow
        mid      = sow + duration/2.0
        new(id, orbit, image, opportunity, sow, mid, eow, duration)
    end
end

end # Close DataStructures