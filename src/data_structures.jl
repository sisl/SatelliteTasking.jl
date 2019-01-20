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
- `area_drag::Real`: Wetted wind-facing area [m^2]
- `coef_drag::Real`: Coefficient of drag
- `area_srp::Real`: Area illuminated by solar radiation pressure [m^2]
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
    reward::Float64
    id::UUID

    function Image(lon::Real, lat::Real ; id::UUID=uuid4(),
                    look_angle_min=5.0::Real, look_angle_max=55.0::Real, 
                    require_zero_doppler=false::Bool, reward=1.0::Real, use_degrees=true::Bool)

        # Compute ECEF coordinates of scene center
        ecef = sGEODtoECEF([lon, lat], use_degrees=use_degrees::Real)

        # Generate new image object
        new(lon, lat, ecef, look_angle_min, look_angle_max,
            require_zero_doppler, reward, id)
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
- `orbit_id::UUID` Identification ID of the orbit which produced this collect
- `image_id::UUID` Identification UUID of the image associated with this collect
- `sow::Epoch` Start of acquisition window 
- `eow::Epoch` End of possible acquisition window
- `duration::` Length of acquisition window
"""
mutable struct Collect
    id::UUID
    orbit_id::Integer
    image_id::UUID
    sow::Epoch
    eow::Epoch
    duration::Float64

    function Collect(sow::Epoch, eow::Epoch;
                         id=uuid4()::UUID, orbit_id=0::Integer, image_id=uuid4()::UUID)

        duration = eow - sow
        new(id, orbit_id, image_id, sow, eow, duration)
    end
end

end # Close DataStructures