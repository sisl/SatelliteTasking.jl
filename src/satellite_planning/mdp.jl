__precompile__(true)
module MDP

# Julia Packages
using Random
using Printf

using SatelliteDynamics.Time: Epoch, mjd
using SatelliteDynamics.OrbitDynamics: eclipse_conical
using SatelliteTasking.DataStructures: Image, Opportunity, Orbit


########################
# General MDP Planning #
########################

# Default spacecraft model
sc_model_default = Dict{String,Float64}(
    "power_generation" => 0.0,
    "power_draw_downlink" => 0.0,
    "power_draw_imaging"  => 0.0,
    "data_generation_backorbit" => 0.0,
    "data_generation_imaging" => 0.0,
    "data_downlink" => 0.0
)

export MDPState
"""
MDPState 

Attributes:
- `time::Epoch` Current planning time
- `locations::BitArray{1}` Boolean bit array of whether a specific image has 
    been vistited. Index indicates image.
"""
mutable struct MDPState
    time::Epoch
    last_action::Union{Opportunity, Symbol}
    images::Array{Image, 1}
    power::Float64
    data::Float64
    done::Bool
end

function Base.show(io::IO, state::MDPState)

    s = @sprintf "MDPState(time: %s, power: %.2f, data: %.2f, image_count: %d)" state.time state.power state.data length(state.images)

    print(io, s)
end

# Comparison operators for MDP state, relies on same allocation of opportunities
# function Base.:(==)(state_left::MDPState, state_right::MDPState)
#     return ((state_left.time == state_right.time) && (state_left.locations == state_right.locations))
# end

# function Base.:(!=)(state_left::MDPState, state_right::MDPState)
#     return ((state_left.time != state_right.time) || (state_left.locations != state_right.locations))
# end

# State constructor from image length
function mdp_compute_actions(state::MDPState, collect_opportunities::Array{Opportunity, 1}; 
            constraint_list::Array{Function, 1}=Function[], breadth=0::Real)

    # Limit search to future opportunities
    future_opportunities = filter(x -> x.sow > state.time, collect_opportunities)

    # Populate list of possible actions
    actions = Union{Opportunity, Symbol}[]

    # Add image collection opportunities
    for end_opportunity in future_opportunities
        # Exit early if we've found enough actions for our breadth search
        if breadth > 0 && length(actions) == breadth
            break
        end

        # Check if valid transition
        valid_transition = true
        if typeof(state.last_action) == Opportunity
            # Check constraints 
            for constraint in constraint_list
                # If not valid transititon break early
                if valid_transition == false
                    continue
                end

                # Use logical AND ot evaluate continued feasibility
                valid_transition = valid_transition && constraint(state.last_action, end_opportunity)
            end

            # If valid add action to list of valid state actions
            if valid_transition
                push!(actions, end_opportunity)
            end
        else
            # Last action was a Symbol
            if state.last_action == :SUNPOINT
                push!(actions, end_opportunity)
            end
        end
    end

    # Only add forced actions if other actions are left
    if length(future_opportunities) > 0
        # Add sunpoint action as valid action
        push!(actions, :SUNPOINT)
    end

    return actions
end

"""
Compute the reward for being in a given state.
"""
function mdp_reward(state::MDPState, action::Union{Symbol, Opportunity}; alpha::Real=0.0)

    # Bias function

    r = 0.0
    
    if typeof(action) == Opportunity
        if !(action.location in state.images)
            r = action.location.reward*(1.0 + alpha)^(-abs(action.sow - state.time))
        end
    end

    # Check power state
    if state.power <= 0.0
        r += -1000000
    end

    return r
end

function extract_images(opportunities::Array{Opportunity, 1})
    images = Image[]

    for opp in opportunities
        if !(opp.location in images)
            push!(images, opp.location)
        end
    end

    return images
end

################################
# Deterministic Forward Search #
################################

function reachable_states(state::MDPState, action::Union{Symbol, Opportunity}, 
            collect_opportunities::Array{Opportunity, 1}, orbit::Orbit,
            sc_model::Dict{String, Float64}=sc_model_default,
            next_action_time::Dict{Epoch, Epoch}=Dict{Epoch, Epoch}())
    # Compute all possible reachable states given acition
    # (In this case it's deterministic) 

    time0 = state.time
    time  = state.time

    # Spacecraft resources
    power = state.power
    data  = state.data

    # Shallow copy observed images
    images = copy(state.images)

    if action == :SUNPOINT
        # Sunpoint transition
        # idx  = findfirst(o -> o.sow > time, collect_opportunities)
        # time = collect_opportunities[idx].sow
        time = next_action_time[time0]

        # Generate power for duration of step
        power += sc_model["power_generation"]*(time - time0)

        # Recharge
        # println("Investigating sunpoint state.")
        # println("Power before: $power")
        # power += sc_model["power_generation"]*(collect_opportunities[idx].sow-state.time)
        # println("Power after: $power")
    else
        # Advance time
        time = action.sow

        # Decrement power
        power += action.location.collect_duration*sc_model["power_draw_imaging"]

        # Increment data 
        data += action.location.collect_duration*sc_model["data_generation_imaging"]

        push!(images, action.location)
    end
    sp = MDPState(time, action, images, power, data, false)

    return MDPState[sp]
end

function forward_search(state::MDPState, d::Int, 
            collect_opportunities::Array{Opportunity, 1}, orbit::Orbit,
            sc_model::Dict{String, Float64}=sc_model_default; 
            constraint_list::Array{Function, 1}=Function[],
            next_action_time::Dict{Epoch, Epoch}=Dict{Epoch, Epoch}(),
            depth::Real=3, breadth::Real=3, gamma::Real=1.0, alpha::Real=0.0)

    if d == 0
        return :SUNPOINT, 0
    end

    # Initialize current optimal action and value
    astar, vstar = :DONE, -Inf 

    for a in mdp_compute_actions(state, collect_opportunities, constraint_list=constraint_list, breadth=breadth)
        v = mdp_reward(state, a, alpha=alpha)
        
        for sp in reachable_states(state, a, collect_opportunities, orbit, sc_model, next_action_time)
            # Continue forward search of space
            ap, vp = forward_search(sp, d-1, collect_opportunities, orbit, sc_model,
                        constraint_list=constraint_list,
                        next_action_time=next_action_time,
                        depth=depth, breadth=breadth, gamma=gamma)
            
            # Update state reward.
            v = v + gamma*vp
        end

        if v > vstar
            astar, vstar = a, v
        end
    end

    return astar, vstar
end

function mdp_forward_step(state::MDPState, action::Union{Opportunity, Symbol}, 
            collect_opportunities::Array{Opportunity, 1}, orbit::Orbit,
            next_action_time::Dict{Epoch, Epoch}=Dict{Epoch, Epoch}(),
            sc_model::Dict{String, Float64}=sc_model_default)
    
    # Current state time
    time0 = state.time
    time  = state.time

    # Resources
    power = state.power
    data  = state.data

    # Shallow copy observed images
    images = copy(state.images)

    # Transition based on action type
    if action == :SUNPOINT
        # Sunpoint transition
        # idx  = findfirst(o -> o.sow > time, collect_opportunities)
        # time = collect_opportunities[idx].sow
        time = next_action_time[time0]

        # Generate power for duration of step
        power += sc_model["power_generation"]*(time - time0)

    elseif typeof(action) == Opportunity
        # Opportunity Transition

        # Advance time
        time = action.sow

        # Decrement power
        power += action.location.collect_duration*sc_model["power_draw_imaging"]

        # Increment data 
        data += action.location.collect_duration*sc_model["data_generation_imaging"]

        # Add newly observed image to list
        push!(images, action.location)

    else
        throw(ErrorException("Unknown action type $(string(action))"))
    end

    # Initialize next state
    state_next = MDPState(time, action, images, power, data, state.done)

    return state_next
end

# MDP forward search algorithm
export mdp_solve_forward_search
"""
Solve MDP using basic forward search algorithm.
"""
function mdp_solve_forward_search(opportunities::Array{Opportunity, 1}, 
            constraint_list::Array{Function, 1}=Function[];
            depth::Real=3, breadth::Real=3, gamma::Real=1.0, alpha::Real=0.0, 
            sc_model::Dict{String, Float64}=sc_model_default)

    # Sort Opportunities
    sort!(opportunities, by = x -> x.sow)
    
    # Extract image list
    images = extract_images(opportunities)

    # Next Action list 
    next_action_time = Dict{Epoch, Epoch}()
    for idx in 1:(length(opportunities)-1)
        next_action_time[opportunities[idx].sow] = opportunities[idx+1].sow
    end

    # Orbit
    orbit = opportunities[1].orbit

    # Set Initial state
    init_opp = opportunities[collect(keys(opportunities))[findmin(collect([o.sow for o in opportunities]))[2]]]
    state    = MDPState(init_opp.sow, init_opp, Image[init_opp.location], 1.0, 0.0, false)

    # Store plan (action history) and state history
    plan   = Union{Opportunity, Symbol, Nothing}[]
    states = Union{MDPState, Nothing}[]

    # Take initial step
    action, value = forward_search(state, depth, opportunities, orbit, sc_model, 
                        constraint_list=constraint_list, depth=depth, 
                        breadth=breadth, gamma=gamma, alpha=alpha, next_action_time=next_action_time)

    # Add initial state to search
    push!(states, state)
    push!(plan, action)

    # Transition state forward
    state = mdp_forward_step(state, action, opportunities, orbit, next_action_time, sc_model)

    # println("State: $(string(state))")

    while true
        # Compute next action
        action, value = forward_search(state, depth, opportunities, orbit, sc_model, 
                            constraint_list=constraint_list, depth=depth, 
                            breadth=breadth, gamma=gamma, alpha=alpha, next_action_time=next_action_time)

        # println("Next action: $(string(action))")

        if action == :DONE
            state.done = true
            break
        end

        # Add state and action to plan 
        push!(states, state)
        push!(plan, action)

        # Transition state forward
        state = mdp_forward_step(state, action, opportunities, orbit, next_action_time, sc_model)

        # println("State: $(string(state))")
    end

    # Compute plan reward
    image_list = copy(state.images)
    reward     = 0
    for image in image_list
        reward += image.reward
    end
        
    return states, plan, reward, image_list
end

###########################
# Monte Carlo Tree Search #
###########################

# function mcts_generate_rollout(opportunities::Array{Opportunity, 1}, constraint_list::Array{Function, 1}; breadth::Integer=10)

#     rollout_policy = Dict{Opportunity, Array{Tuple{Opportunity, <:Real}}}()

#     for start_opportunity in opportunities

#         # Initialize Array
#         rollout_policy[start_opportunity] = Array{Tuple{Opportunity, <:Real}, 1}[]

#         # Populate list of possible actions
#         actions = mdp_compute_actions(start_opportunity, opportunities, constraint_list, breadth=breadth)

#         num_actions = length(actions)
#         for a in actions
#             push!(rollout_policy[start_opportunity], (a, 1.0/num_actions))
#         end
#     end

#     return rollout_policy
# end

# function mcts_get_optimal_action(state::MDPState, N::Dict{Opportunity, Dict{Opportunity, Int32}}, Q::Dict{Opportunity, Dict{Opportunity, Float64}}; c=1.0)
#     # Find optimal action
#     amax, vmax = nothing, -Inf

#     Ns = sum([N[state.time][as] for as in keys(Q[state.time])])
#     for a in keys(Q[state.time])

#         if Ns == 0
#             # If Ns is zero use this otherwise log(Ns) = -Inf
#             v = Q[state.time][a]
#         else
#             v = Q[state.time][a] + c*sqrt(log(Ns)/N[state.time][a])
#         end

#         if v > vmax
#             amax, vmax = a, v
#         end
#     end

#     return amax
# end

# function mcts_tree_search(state::MDPState, d::Integer, rollout_policy::Dict{Opportunity, Array{Tuple{Opportunity, <:Real}}};
#             opportunities::Array{Opportunity, 1},
#             probabilities::Union{Dict{Opportunity, <:Real}, Nothing}=nothing,
#             constraint_list::Array{Function, 1}, 
#             position_lookup::Dict{Image, <:Integer},
#             image_lookup::Dict{<:Integer, Image}, 
#             T::Array{Opportunity, 1},
#             N::Dict{Opportunity, Dict{Opportunity, Int32}},
#             Q::Dict{Opportunity, Dict{Opportunity, Float64}},
#             c::Real=1.0, gamma::Real=0.95, breadth::Integer=10)
   
#     if d == 0
#         return 0.0
#     end

#     # If state not seen, add it and rollout policy 
#     if !(state.time in T)
#         # println("Encountered new state $(string(UInt64(pointer_from_objref(state.time)), base=16)).")

#         # Initiallies N & Q arrays
#         N[state.time] = Dict{Opportunity, Int32}()
#         Q[state.time] = Dict{Opportunity, Float64}()

#         for a in mdp_compute_actions(state.time, opportunities, constraint_list, breadth=breadth)
#             # Initialize state action reward value if not
#             N[state.time][a] = 0
#             Q[state.time][a] = 0
#         end

#         # Add state to observed states
#         push!(T, state.time)

#         # Return reward for state using rollout_policy
#         return mcts_rollout(state, d, rollout_policy, 
#                 gamma=gamma, 
#                 probabilities=probabilities, 
#                 position_lookup=position_lookup, 
#                 image_lookup=image_lookup)
#     else    
#         # println("Encountered state $(string(UInt64(pointer_from_objref(state.time)), base=16)) - $(sum([N[state.time][as] for as in keys(Q[state.time])])) times previously.")
#     end

#     # Find optimal action
#     amax = mcts_get_optimal_action(state, N, Q, c=c)    

#     # No future actions return nothing
#     if amax == nothing
#         return 0.0
#     end

#     # Simulate next state
#     next_state = mdp_forward_step(state, amax, position_lookup, probabilities=probabilities)
       
#     # Reward for action
#     # r = amax.location.reward*next_state.locations[position_lookup[amax.location]]
#     r = mdp_reward(next_state, image_lookup, position_lookup)

#     q = r + gamma*mcts_tree_search(next_state, d-1, rollout_policy, 
#                     opportunities=opportunities, probabilities=probabilities,
#                     constraint_list=constraint_list,
#                     position_lookup=position_lookup,
#                     image_lookup=image_lookup, 
#                     T=T, Q=Q, N=N, c=c, gamma=gamma, breadth=breadth)
    
#     N[state.time][amax] = N[state.time][amax] + 1
#     Q[state.time][amax] = Q[state.time][amax] + (q - Q[state.time][amax])/N[state.time][amax]

#     return q
# end

# function sample_rollout(state::MDPState, rollout_policy::Dict{Opportunity, Array{Tuple{Opportunity, <:Real}}})
#     # Generate random number
#     p = rand()

#     # Set cumulative probability
#     cum_prob = 0.0

#     for ap in rollout_policy[state.time]
#         cum_prob += ap[2]

#         if p <= cum_prob
#             return ap[1]
#         end
#     end

#     return nothing
# end

# function mcts_rollout(state, d, rollout_policy::Dict{Opportunity, Array{Tuple{Opportunity, <:Real}}}; 
#             probabilities::Union{Dict{Opportunity, <:Real}, Nothing}=nothing,
#             position_lookup::Dict{Image, <:Integer},
#             image_lookup::Dict{<:Integer, Image}, 
#             gamma::Real=0.95)
            
#     if d == 0
#         return 0
#     end

#     # Sample rollout policy
#     action = sample_rollout(state, rollout_policy)

#     # If no possible future states return 0
#     if action == nothing
#         return 0
#     end

#     # Simulate next state
#     next_state = mdp_forward_step(state, action, position_lookup, probabilities=probabilities)
    
#     # Reward for action
#     # r = action.location.reward*next_state.locations[position_lookup[action.location]]
#     r = mdp_reward(next_state, image_lookup, position_lookup)

#     return r + gamma*mcts_rollout(next_state, d-1, rollout_policy, 
#                         gamma=gamma,
#                         probabilities=probabilities, 
#                         position_lookup=position_lookup,
#                         image_lookup=image_lookup)
# end


# export mdp_solve_mcts
# """
# Solve MDP witht Monte Carlo Tree Search
# """
# function mdp_solve_mcts(opportunities::Array{Opportunity, 1}, 
#     constraint_list::Array{Function, 1}, probabilities::Union{Dict{Opportunity, <:Real}, Nothing}=nothing;
#     depth::Real=10, breadth::Integer=10, gamma::Real=0.95, c::Real=0.75, max_iterations::Integer=10)

#     # Extract image list
#     images = extract_images(opportunities)

#     # Compute image lookup
#     image_lookup, position_lookup = create_lookups(images)

#     # Initialize MCTS variables
#     plan = Union{Opportunity, Nothing}[]
#     T = Opportunity[]                                    # Visited States 
#     N = Dict{Opportunity, Dict{Opportunity, Int32}}()    # State-Action Exploration Count
#     Q = Dict{Opportunity, Dict{Opportunity, Float64}}()  # State-Action Reward

#     # Generate rollout policy
#     rollout_policy = mcts_generate_rollout(opportunities, constraint_list, breadth=breadth)

#     # Set Initial state
#     init_opp = opportunities[collect(keys(opportunities))[findmin(collect([o.sow for o in opportunities]))[2]]]
#     state    = MDPState(init_opp, images)

#     # println("Current state: $(state.time.sow)")

#     # Run MCTS for fix number of iterations
#     for i in 1:max_iterations
#         # println("Tree search Iteration $i")
#         q = mcts_tree_search(state, depth, rollout_policy, 
#                     opportunities=opportunities,
#                     probabilities=probabilities,
#                     constraint_list=constraint_list,
#                     position_lookup=position_lookup,
#                     image_lookup=image_lookup,
#                     T=T,
#                     Q=Q,
#                     N=N,
#                     c=c,
#                     gamma=gamma,
#                     breadth=breadth)
#     end

#     # Get Optimal action for state
#     action = mcts_get_optimal_action(state, N, Q, c=0)

#     # Add action to plan
#     push!(plan, action)

#     # Transition state forward probabilistically
#     state = mdp_forward_step(state, action, position_lookup, probabilities=probabilities)
#     # println("Current state: $(state.time.sow)")

#     while true
#         # Run MCTS for fix number of iterations
#         for i in 1:max_iterations
#             # println("Tree search Iteration $i")
#             q = mcts_tree_search(state, depth, rollout_policy, 
#                         opportunities=opportunities,
#                         probabilities=probabilities,
#                         constraint_list=constraint_list,
#                         position_lookup=position_lookup,
#                         image_lookup=image_lookup,
#                         T=T,
#                         Q=Q,
#                         N=N,
#                         c=c,
#                         gamma=gamma,
#                         breadth=breadth)
#         end

#         # Get Optimal action for state
#         action = mcts_get_optimal_action(state, N, Q, c=0)

#         if action == nothing
#             @debug "No next action at time: $(state.time) Stopping search. $(mdp_compute_actions(state.time, opportunities, constraint_list, breadth=breadth))"
#             # Exit early
#             break
#         end

#         # Add action to plan
#         push!(plan, action)

#         # Transition state forward probabilistically
#         state = mdp_forward_step(state, action, position_lookup, probabilities=probabilities)
#         @debug "Current state: $(state.time.sow)"
#     end


#     # Compute plan reward
#     reward = 0
#     image_list = Image[]
    
#     for (i, image_collected) in enumerate(state.locations)
#         if image_collected
#             image   = image_lookup[i]
#             reward += image.reward
#             push!(image_list, image)
#         end
#     end
        
#     return plan, reward, image_list
# end

end # End MDP module