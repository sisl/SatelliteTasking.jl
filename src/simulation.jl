# Exports
# export simulate_orbits

# export simulate_orbits
# """
# Simulate multiple satellite orbits with perturbed initial conditions.

# First simulates the true orbit without any perturbation error then simulates the
# given number of orbits each with randomly sampled perturbations on the initial
# state.

# Positional Arguments:
# - `num_orbits::Integer` Number of perturbed orbits to simulate
# - `epc0::Epoch` Initial Epoch of propagation
# - `epcf::Epoch` Final Epoch of propagation
# - `eci0::Array{<:Real, 1}` Initial inertial state at initial Epoch [m; m/s]
# - `pmean::Array{<:Real, 1}` Vector of mean bias of initial state error from `eci0`
# - `pstd::Array{<:Real, 1}` Vector of standard deviation of initial state error from `eci0`

# Keyword Arguments:
# - `timestep::Real` Timestep used for propgation output
# - `dtmax::Real` Maximum allowable integrator step
# - `mass::Real` Satellite mass [kg]
# - `area_drag::Real` Wetted wind-facing area [m^2]
# - `coef_drag::Real` Coefficient of drag
# - `area_srp::Real` Area illuminated by solar radiation pressure [m^2]
# - `coef_srp::Real` Coefficient of reflectivity
# - `n_grav::Integer` Gravity model degree
# - `m_grav::Integer` Gravity model order
# - `drag::Bool` Enable drag perturbation in propgation
# - `srp::Bool` Enable solar radiation pressure perturbation in propgation
# - `moon::Bool` Enable thirdbody moon perturbation in propgation
# - `sun::Bool` Enable thirdbody sun perturbation in propgation
# - `relativity::Bool` Enable relativity perturbation in propgation

# Returns:
# - `true_orbit::Orbit` True unperturbed orbit
# - `perturbed_orbits::Array{Orbit, 1}` Orbit objects containing simulated orbits with perturbed initial conditions.
# """
# function simulate_orbits(num_orbits::Integer, 
#                         epc0::Epoch, epcf::Epoch, eci0::Array{<:Real, 1}, 
#                         pmean::Array{<:Real, 1}, pstd::Array{<:Real, 1};
#                         timestep=1.0::Real, dtmax=30.0::Real,
#                         mass=100.0::Real, area_drag=1.0::Real, coef_drag=2.3::Real, 
#                         area_srp=1.0::Real, coef_srp=1.8::Real, 
#                         n_grav=20::Integer, m_grav=20::Integer, 
#                         drag=true::Bool, srp=true::Bool, moon=true::Bool, 
#                         sun=true::Bool, relativity=true::Bool)

#     # Test to make sure inputs are valid
#     if all(pstd .== 0.0)
#         throw(ArgumentError("There must be at least non-zero element in standard deviation vector"))
#     end

#     # Simulate True Orbit
#     @debug "Simulating truth orbit"
#     true_orbit = Orbit(epc0, epcf, eci0, timestep=timestep, dtmax=max(dtmax, timestep),
#                         mass=mass, area_drag=area_drag, coef_drag=coef_drag, 
#                         area_srp=area_srp, coef_srp=coef_srp, 
#                         n_grav=n_grav, m_grav=m_grav, 
#                         drag=drag, srp=srp, moon=moon, sun=sun, relativity=relativity)

#     # Create covatiance matrix from standard deviations
#     covar    = diagm(0 => pstd.^2)
#     dist_eci = MvNormal(pmean, covar) # Multivariate normal distribution generating errors

#     # Simulate perturbed orbits
#     perturbed_orbits = Orbit[]
#     eci_errors       = rand(dist_eci, num_orbits)

#     for i in 1:num_orbits
#         println("Simulating perturbed orbit $i")
#         eci_perturbed = eci0 + eci_errors[:, i]
#         orb = Orbit(epc0, epcf, eci_perturbed, timestep=timestep, dtmax=max(dtmax, timestep),
#                     mass=mass, area_drag=area_drag, coef_drag=coef_drag, 
#                     area_srp=area_srp, coef_srp=coef_srp, 
#                     n_grav=n_grav, m_grav=m_grav, 
#                     drag=drag, srp=srp, moon=moon, sun=sun, relativity=relativity)
        
#         push!(perturbed_orbits, orb)
#     end

#     return true_orbit, perturbed_orbits, eci_errors
# end