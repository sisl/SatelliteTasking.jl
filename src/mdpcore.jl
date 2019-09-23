###########################
# Action Space Definition #
###########################

# Terminal States
# POMDPs.isterminal(mdp::TaskingMDP, s::TaskingState) = (length(mdp.action_lookup[a.sow]) == 0)
POMDPs.isterminal(problem::SatPlanningProblem, s::SatMDPState) = s.done
POMDPs.discount(problem::SatPlanningProblem) = 1.0

# Generate next state from state and action
function POMDPs.generate_s(problem::SatPlanningProblem, state::SatMDPState, action::Opportunity, rng::AbstractRNG)
    # # If done return same state, just done
    # if typeof(action) == Done
    #     return SatMDPState(time=state.time,
    #             last_action=Done(t_start=state.time),
    #             last_cdo_action=state.last_cdo_action,
    #             requests=copy(state.requests),
    #             power=state.power,
    #             data=state.data)
    # end

    # Current state time
    time0 = state.time
    time  = state.time

    # Store last viable action
    last_cdo_action = state.last_cdo_action

    # Resource
    power = state.power
    data  = state.data

    # Copy observed requests
    requests = copy(state.requests)

    # Finish state
    done = false

    # Advance time to next action
    time = action.t_start

    # Update data generate
    power_generated = 0.0
    data_generated  = problem.spacecraft[1].datagen_backorbit*(time - time0)

    if typeof(action) == Noop
        # Do nothing if Noop
    elseif typeof(action) == Sunpoint
        # Charge for duration of sunpoint
        power_generated += problem.spacecraft[1].powergen_sunpoint*(time - time0)

    elseif typeof(action) == Collect
        # Update last action
        last_cdo_action = action

        collect_generation = action.duration * problem.spacecraft[1].datagen_image

        if data + collect_generation < 1.0
            # Only perform collect if we have capacilty
            data_generated  += collect_generation
            power_generated += action.duration * problem.spacecraft[1].powergen_image

            push!(requests, action.location)
        end

    elseif typeof(action) == Contact
        # Update last action
        last_cdo_action = action

        data_generated  += action.duration * problem.spacecraft[1].datagen_downlink
        power_generated += action.duration * problem.spacecraft[1].powergen_downlink
    else
        throw(ErrorException("Unknown Action Type: $(string(action))."))
    end

    # Update power and data numbers
    power += power_generated
    data  += data_generated

    # Ensure resource limits are enforced
    if power > 1.0 power = 1.0 end

    if power < 0.0 
        power = 0.0 
        done = true
    end

    if data > 1.0 data = 1.0 end
    if data < 0.0 data = 0.0 end

    return SatMDPState(time=time,
            last_action=action,
            last_cdo_action=last_cdo_action,
            requests=requests,
            power=power,
            data=data)
end

# State reward function
function POMDPs.reward(problem::SatPlanningProblem, state::SatMDPState, action::Opportunity)
    reward = 0.0

    # Update data generate
    power_generated = 0.0
    data_generated  = problem.spacecraft[1].datagen_backorbit*(action.t_start - state.time)

    if typeof(action) == Collect
        # Resources consumed by action
        data_generated  += action.duration * problem.spacecraft[1].datagen_image
        power_generated += action.duration * problem.spacecraft[1].powergen_image
        
        # Only reward if we have enough spare data space
        if (state.data + data_generated < 1.0)
            if !(action.location in state.requests)
                reward += action.reward
            end
            # elseif problem.solve_allow_repeats == true
            #     reward += action.reward*0.25
            # end
        else
            # Penalize collection when full up
            # reward -= action.reward
        end

    elseif typeof(action) == Contact
        # reward += 0.01*action.duration

        # Update data generation
        data_generated  += action.duration * problem.spacecraft[1].datagen_downlink
        power_generated += action.duration * problem.spacecraft[1].powergen_downlink
    elseif typeof(action) == Noop
        reward += 0.0
    elseif typeof(action) == Sunpoint
        # Charge for duration of sunpoint
        power_generated += problem.spacecraft[1].powergen_sunpoint*(action.t_start - state.time)

        # Have tiny bit of reward for charging
        # reward += (action.t_start - state.time)/state.last_cdo_action.duration
        reward += 0.0001*(action.t_start - state.time)
    end

    # Update resources to see if violations occur
    power = state.power + power_generated
    data  = state.data  + data_generated

    # Penalize running out of power
    if power <= 0
        reward -= 10000
    end

    # Don't penalize overcharge
    # if power > 1
    #     reward -= 10.0
    # end

    # Don't penalize full capacity
    # if data > 1
    #     reward -= 10.0
    # end

    # Don't penalize having everything down
    # if data < 0
    #     reward -= 10.0
    # end

    return reward
end

# Action Space
POMDPs.actions(problem::SatPlanningProblem) = problem.actions
function POMDPs.actions(problem::SatPlanningProblem, state::SatMDPState)
    return problem.lt_feasible_actions[(state.last_cdo_action.id, state.current_action.id)]
end

# Initial state funciton
function POMDPs.initialstate(problem::SatPlanningProblem, rng::AbstractRNG)

    # Get Initial Opportunity
    init_opp = problem.opportunities[1]
    
    # Create Initial State
    state = SatMDPState(time=init_opp.t_start, last_action=init_opp, last_cdo_action=init_opp)

    return state
end