# Exports
export MDPState
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
    return hash(state.time, hash(length(state.images), hash(state.power, hash(state.data, hash(state.done, hash(:Epoch, h))))))
end

function Base.show(io::IO, mdp::MDPState)

    s = @sprintf "MDPState(Time: %s, requests: %d, power: %.2f, data: %.2f, done: %s)" string(mdp.time) length(mdp.requests) mdp.power mdp.data string(mdp.done)

    print(io, s)
end

"""
Find actions for planning problem given state and time horizon.

Uses binary search of sorted opportunity list to find element.
"""
function mdp_find_actions(problem::PlanningProblem, state::MDPState, horizon::Real)
    # Get Index of Current time from State
    t_start = state.last_action.t_start

    l = 1
    r = length(problem.opportunities)
    m = floor((l + r) / 2)
    while l <= r
        if problem.opportunities[m].t_start < horizon
            l = m + 1
        elseif problem.opportunities[m].t_start > horizon
            r = m - 1
        else
            return m
        end
    end
    return nothing
end

function mdp_reward(problem::PlanningProblem, state::MDPState, action::Opportunity)
    reward = 0.0

    if typeof(action) == Collect && !(action.location in state.requests)

        # Resources consumed by action
        data_generated  = action.duration*problem.spacecraft[1].datagen_image
        power_generated = action.duration*problem.spacecraft[1].powergen_image

        # Only reward if we haven't overfilled
        if (state.data + data_generated < 1.0)
            reward += action.reward*(1.0 + problem.reward_alpha)^(-abs(action.t_start - state.time))
        end
    elseif typeof(action) == Contact
        reward += 0.01*action.duration
    elseif typeof(action) == Noop
        reward += 0.0
    end

    if state.power <= 0
        reward -= 1000000
    end

    return reward
end

function mdp_step(problem::PlanningProblem, state::MDPState, action::Opportunity)
    # If done return same state, just done
    if typeof(action) == Done
        return MDPState(time=state.time,
                last_action=Done(),
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

    if typeof(action) == Noop
        # Do nothing if Noop
    elseif typeof(action) == Sunpoint
        # Charge
    elseif typeof(action) == Collect
        # Update last action
        last_action = action

        # Only perform collect if we have capacilty
        data_generated  = action.duration*problem.spacecraft[1].datagen_image
        power_generated = action.duration*problem.spacecraft[1].powergen_image

        if data + data_generated < 1.0
            power += power_generated
            data  += data_generated

            push!(requests, action.location)
        end

    elseif typeof(action) == Contact
        # Update last action
        last_action = action

        data_generated  = action.duration*problem.spacecraft[1].datagen_downlink
        power_generated = action.duration*problem.spacecraft[1].powergen_downlink
    else
        throw(ErrorException("Action $(string(action)) is not of a known type."))
    end

    # Ensure resource limits are enforced
    if power > 1.0 power = 1.0 end
    if power < 0.0 power = 0.0 end

    if data > 1.0 data = 1.0 end
    if data < 0.0 data = 0.0 end

    # If no other possible future actions problem state is done 
    if length(filter(x -> x.t_start > state.time, problem.opportunities)) == 0
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
    fopps = filter(x -> x.t_start > state.time, problem.opportunities)

    if length(fopps) == 0
        return Opportunity[Done()]
    end

    # List of all possible actions
    actions = Opportunity[]

    # Populate list of possible actions
    for opp in fopps
        if problem.solve_breadth > 0 && length(actions) == problem.solve_breadth
            break
        end

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
            push!(actions, opp)
        end
    end

    # Add default actions at the end not to mess width action breadth-limit
    push!(actions, Noop(t_start=fopps[1].t_start))
    push!(actions, Sunpoint(t_start=fopps[1].t_start))

    # Return set of possible actions
    return actions
end

function mdp_reachable_states(problem::PlanningProblem, state::MDPState, action::Opportunity)
    return MDPState[mdp_step(problem, state, action)]
end