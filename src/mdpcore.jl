# Exports
export MDPState
export bs_find_firstgte
export find_actions
export mdp_reward
export mdp_step
export mdp_state_actions

@with_kw struct MDPState
    time::Union{Epoch, Float64}
    last_action::Opportunity
    requests::Array{Request, 1} = Request[]
    power::Float64 = 1.0
    data::Float64 = 0.0
    done::Bool = false
end

function Base.hash(state::MDPState, h::UInt)
    return hash(state.time, hash(length(state.requests), hash(state.power, hash(state.data, hash(state.done, hash(:Epoch, h))))))
end

function Base.show(io::IO, mdp::MDPState)

    s = @sprintf "MDPState(Time: %s, requests: %d, power: %.2f, data: %.2f, done: %s)" string(mdp.time) length(mdp.requests) mdp.power mdp.data string(mdp.done)

    print(io, s)
end

# Comparison operators for MDP state, relies on same allocation of opportunities
function Base.:(==)(state_left::MDPState, state_right::MDPState)
    return ((state_left.time == state_right.time) &&
            (state_left.last_action == state_right.last_action) &&
            (state_left.requests == state_right.requests) &&
            (state_left.power == state_right.power) &&
            (state_left.data == state_right.data) &&
            (state_left.done == state_right.done))
end

function Base.:(!=)(state_left::MDPState, state_right::MDPState)
    return ((state_left.time != state_right.time) ||
            (state_left.last_action != state_right.last_action) ||
            (state_left.requests != state_right.requests) ||
            (state_left.power != state_right.power) ||
            (state_left.data != state_right.data) ||
            (state_left.done != state_right.done))
end

"""
Binary search method to find the first opportunity with an index greater than
the target simulation time (Epoch or elapsed time)
"""
function bs_find_firstgte(problem::PlanningProblem, target::Union{Epoch, Real})
    # Set search indices
    l = 1
    r = length(problem.opportunities)
    m = 0
    
    # Perform binary search to get first opportunity with start time
    # greater or equal to the target value
    while l < r
        m = floor(Int, (l+r)/2)
        # println("l: $l, r: $r, avg: $((l+r)/2), m: $m, t: $(problem.opportunities[m].t_start)")
        
        if problem.opportunities[m].t_start < target
            l = m + 1
        else
            r = m
        end
    end
    
    return m
end

"""
Find all possible (not necessarily feasible) actions given a state and a
planning horizon.
"""
function find_actions(problem::PlanningProblem, state::MDPState, horizon::Real)
    # Find index of next possible action
    idx_start = bs_find_firstgte(problem::PlanningProblem, state.time)

    # Ensure index is for next action after current state
    while problem.opportunities[idx_start].t_start <= state.time
        idx_start += 1
    end

    # Find index of first action outside of horizon
    idx_end = bs_find_firstgte(problem::PlanningProblem, state.time+horizon)
    
    # Get index of last action inside horizon
    idx_end = idx_end - 1
    
    # Return all opportunities inside planning horizon
    return problem.opportunities[idx_start:idx_end]
end

function mdp_reward(problem::PlanningProblem, state::MDPState, action::Opportunity)
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
                if problem.mdp_reward_scarcity == true
                    reward += action.reward*(1+1.0/(action.nr+1.0))
                else
                    reward += action.reward
                end
            elseif problem.solve_allow_repeats == true
                reward += action.reward*0.25
            end
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

        reward += 0.01
    end

    # Update resources to see if violations occur
    power = state.power + power_generated
    data  = state.data  + data_generated

    # Penalize resource over usage
    if power <= 0
        # reward -= 1000
        reward -= 1
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

function mdp_step(problem::PlanningProblem, state::MDPState, action::Opportunity)
    # If done return same state, just done
    if typeof(action) == Done
        return MDPState(time=state.time,
                last_action=Done(t_start=state.time),
                requests=copy(state.requests),
                power=state.power,
                data=state.data,
                done=true)
    end

    # Current state time
    time0 = state.time
    time  = state.time

    # Store last viable action
    last_action = state.last_action

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
        last_action = action

        collect_generation = action.duration * problem.spacecraft[1].datagen_image

        if data + collect_generation < 1.0
            # Only perform collect if we have capacilty
            data_generated  += collect_generation
            power_generated += action.duration * problem.spacecraft[1].powergen_image

            push!(requests, action.location)
        end

    elseif typeof(action) == Contact
        # Update last action
        last_action = action

        data_generated  += action.duration * problem.spacecraft[1].datagen_downlink
        power_generated += action.duration * problem.spacecraft[1].powergen_downlink
    else
        throw(ErrorException("Action $(string(action)) is not of a known type."))
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

    # If no other possible future actions problem state is done 
    if length(find_actions(problem, state, problem.solve_horizon)) == 0
        done = true
    end

    return MDPState(time=time,
            last_action=last_action,
            requests=requests,
            power=power,
            data=data,
            done=done)
    
end

function mdp_state_actions(problem::PlanningProblem, state::MDPState)
    # Limit search to future opportunities
    fopps = find_actions(problem, state, problem.solve_horizon)

    # List of all possible actions
    actions = Opportunity[]

    # If done no more actions available
    if typeof(state.last_action) == Done
        return actions
    end

    # Populate list of possible actions
    for opp in fopps
        if problem.solve_breadth > 0 && length(actions) == problem.solve_breadth
            break
        end

        # Note state.last_action should _always_ be a Collect or Contact
        if typeof(state.last_action) == Collect || typeof(state.last_action) == Contact
            # Valuate transition to see if valid
            valid = true

            # Check to see if each transition constraint is valid
            for constraint in problem.constraints
                if valid == false
                    break
                end

                valid = valid && constraint(state.last_action, opp)
            end

            # If valid transition add to edges
            if valid == true
                push!(actions, opp)
            end
        else
            throw(ErrorException("Unexpected last_action type: $(typeof(state.last_action))"))
        end
    end

    # If no actions left to transition to we're done
    if length(actions) == 0 || typeof(state.last_action) == Done
        return Opportunity[Done(t_start=state.time)]
    end

    # Add default actions at the end not to mess width action breadth-limit
    next_action_time = fopps[1].t_start
    # next_action_time = problem.lt_opportunities[state.last_action.id+1].t_start

    push!(actions, Sunpoint(t_start=next_action_time))
    # push!(actions, Noop(t_start=next_action_time))

    # Return set of possible actions
    return actions
end

function mdp_reachable_states(problem::PlanningProblem, state::MDPState, action::Opportunity)
    return MDPState[mdp_step(problem, state, action)]
end