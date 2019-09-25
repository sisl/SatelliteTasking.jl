###########################
# Action Space Definition #
###########################

# Terminal States
function POMDPs.isterminal(problem::SatPlanningProblem, s::SatMDPState)
    done = false

    # Have remaining actions 
    if length(problem.lt_feasible_actions[(s.last_cdo_action.id, s.last_action.id)]) == 0
        done = true
    end

    # if s.power <= 0
    #     done = true
    # end
    
    return done
end
POMDPs.discount(problem::SatPlanningProblem) = problem.solve_discount

# Generate next state from state and action
function POMDPs.gen(problem::SatPlanningProblem, state::SatMDPState, action::Opportunity, rng::AbstractRNG)

    # Current state time
    time0 = state.time
    time  = state.time

    # Store last viable action
    last_cdo_action = state.last_cdo_action

    # Resource
    power = state.power
    data  = state.data

    # Copy observed request_ids
    request_ids = copy(state.request_ids)

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

            push!(request_ids, action.location.id)
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
    end

    if data > 1.0 data = 1.0 end
    if data < 0.0 data = 0.0 end

    # Compute next state and reward for next state
    sp = SatMDPState(time=time,
        last_action=action,
        last_cdo_action=last_cdo_action,
        request_ids=request_ids,
        power=power,
        data=data
    )

    # spr = reward(problem, state, action)

    return (sp=sp,)
end

# State reward function
function POMDPs.reward(problem::SatPlanningProblem, state::SatMDPState, action::Opportunity)
    r = 0.0

    # println("Reward Action: $action")

    # Update data generate
    power_generated = 0.0
    data_generated  = problem.spacecraft[1].datagen_backorbit*(action.t_start - state.time)

    if typeof(action) == Collect
        # Resources consumed by action
        data_generated  += action.duration * problem.spacecraft[1].datagen_image
        power_generated += action.duration * problem.spacecraft[1].powergen_image
        
        # Only reward if we have enough spare data space
        if (state.data + data_generated < 1.0)
            # println("Action Location ID: $(action.location.id) - Requests: $(state.request_ids)")
            if !(action.location.id in state.request_ids)
                # r += action.reward

                # Semi-markov reward update
                t_diff = abs(action.t_start - state.time)
                r += action.reward*POMDPs.discount(problem)^(t_diff)
            end
            # elseif problem.solve_allow_repeats == true
            #     r += action.reward*0.25
            # end
        else
            # Penalize collection when full up
            # r -= 1000
        end

    elseif typeof(action) == Contact
        # r += 0.01*action.duration

        # Update data generation
        data_generated  += action.duration * problem.spacecraft[1].datagen_downlink
        power_generated += action.duration * problem.spacecraft[1].powergen_downlink
    elseif typeof(action) == Noop
        r += 0.0
    elseif typeof(action) == Sunpoint
        # Charge for duration of sunpoint
        power_generated += problem.spacecraft[1].powergen_sunpoint*(action.t_start - state.time)

        # Have tiny bit of reward for charging
        # r += (action.t_start - state.time)/state.last_cdo_action.duration
        r += 0.0001*(action.t_start - state.time)
    end

    # Update resources to see if violations occur
    power = state.power + power_generated
    data  = state.data  + data_generated

    # Penalize running out of power
    if power <= 0.3
        r -= 10000
    end

    # if 0.2 < power <= 0.4
    #     r -= 10
    # end

    # Penalize being full on data
    if data >= 0.75
        r -= 1000
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

    # println("Reward Reward: $r")
    return r
end

# Action Space
POMDPs.actions(problem::SatPlanningProblem) = problem.actions
function POMDPs.actions(problem::SatPlanningProblem, state::SatMDPState)
    return problem.lt_feasible_actions[(state.last_cdo_action.id, state.last_action.id)]
end

# Initial state funciton
function POMDPs.initialstate(problem::SatPlanningProblem, rng::AbstractRNG)

    # Get Initial Opportunity
    init_opp = problem.opportunities[1]
    
    # Create Initial State
    state = SatMDPState(time=init_opp.t_start, last_action=init_opp, last_cdo_action=init_opp)

    return state
end