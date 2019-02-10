__precompile__(true)
module MDP

# Julia Packages
using LinearAlgebra

using SatelliteDynamics.Time: Epoch
using SatelliteTasking.DataStructures: Image, Opportunity, Collect


export MDPState
"""
State 

Attributes:
- `time::Epoch` Current planning time
- `images::Array{Bool,1}` Boolean array of whether a specific image has been collected. Index indicates image.
"""
mutable struct MDPState
    time::Collect
    images::Array{Bool, 1}
end

function MDPState(epc::Collect, image_list::Array{Image,1})    
    return MDPState(epc, falses(length(image_list)))
end

export mdp_construct_actions
"""
Construct the possible action space for all states. The actions are all dependent
on time alone, (not which images have already been collected).

Arguments:
- `collects::Array{Collect, 1}` List of possible actions to take 
- `constraint_list::Array{Collect, 1}` List of constraintt functions that limit feasibility in transition.
- `horizon::Real` Horizon to consider building actions. If 0 and infinite horizon will be used. Other times reduce problem complexity.

Returns:
- `actions::Dict{Epoch, Array{Collect, 1}}` List of possible actions to take 
"""
function mdp_construct_actions(collects::Array{Collect, 1}, constraint_list::Array{Function, 1}; horizon=0::Real)
    
    # Sort Collects to ensure they are in time-asecnding order
    sort!(collects, by = x -> x.sow)

    # Initialize validate transition graph
    actions = Dict{Collect, Array{Collect, 1}}()

    for start_collect in collects
        # List of valid edges/transitions for start_collect
        actions[start_collect] = Array{Collect, 1}[]
        
        for end_collect in collects
            # Decide to insert edge
            if (start_collect == end_collect || 
                start_collect.image == end_collect.image ||
                end_collect.sow < start_collect.eow)
                # Skip insertion if same collect, image, or starts before the current
                # collection ends (no instantaneous manevers)
                continue
            elseif horizon > 0 && end_collect.sow > (start_collect.eow + horizon)
                # Since we know collects are sorted we can break building the
                # transition grarph if the distance from the next to the next start
                # is greater than the look-ahead horizon
                break
            else
                # Set transition valid by default
                valid_transition = true

                for constraint in constraint_list
                    # If not valid transition break early
                    if valid_transition == false
                        continue
                    end

                    # Use logical and to evaulate path feasibility on graph
                    valid_transition = valid_transition && constraint(start_collect, end_collect)
                end

                if valid_transition
                    push!(actions[start_collect], end_collect)
                end
            end
        end
    end
    
    return actions
end

export mdp_image_lookup
"""
Given a list of images create the dictionary lookup of indices 

Arguments:
- `image_list::Array{Image, 1}` List of images (should match MDPState list)

Returns:
- `image_lookup::Dict{Int64, Image}` Dictionary mapping array indices to image objects.
"""
function mdp_image_lookup(image_list::Array{Image,1})
    image_lookup = Dict{Int64, Image}()
    
    for (i, img) in enumerate(image_list)
        image_lookup[i] = img
    end
    
    return image_lookup
end

export compute_collect_probability
"""
Compute the probability that a collect will be possible given the apparent
distribution of collection probabilities.
"""
function compute_collect_probability(collects::Array{Collect}, opportunities::Array{Array{Opportunity,1},1})
    collect_probability = Dict{Collect, Float64}()
    for col in collects
        num_present, num_missing = 0, 0
        for opportunity_list in opportunities
            if length(filter(x -> x.image == col.image && x.sow <= col.sow && x.eow >= col.eow, opportunity_list)) == 1
                num_present += 1
            else
                num_missing += 1
            end
        end
        
        # collect_probability[col] = num_present/(num_missing+num_present)
        collect_probability[col] = 1
    end
    
    return collect_probability
end

export mdp_reward
"""
Compute reward of current MDP state

Arguments:
- `s::MDPState` State of Markov decision process
"""
function mdp_reward(s::MDPState, image_lookup::Dict{<:Integer, Image})
    r = 0.0
    for (i, img) in enumerate(s.images)
        if img
            r += image_lookup[img].reward
        end
    end
    
    return r
end

export mdp_future_states
"""
Compute possible future takes the system may assume given the action.
"""
function mdp_future_states(state::MDPState, action::Collect, int_lookup::Dict{Image, <:Integer})
    failure      = MDPState(state.time, copy(state.images))
    failure.time = action
    
    success = MDPState(failure.time, copy(failure.images))
    success.images[int_lookup[action.image]] = true
    
    return MDPState[success, failure]
end

export mdp_transition
"""
Return transition probabilities given the set of currrent states, action, and possible future states.
"""
function mdp_transition(sp::MDPState, s::MDPState, a::Collect, collect_probabilities::Dict{Collect, Float64})
    p = collect_probabilities[a]
    
    return p
end

export mdp_select_action
"""
MDP forward search algorithm to select the next decision.
"""
function mdp_select_action(s::MDPState, d::Integer, actions::Dict{Collect, Array{Collect, 1}}, 
                            image_lookup::Dict{<:Integer, Image}, int_lookup::Dict{Image, <:Integer},
                            collect_probabilities::Dict{Collect, Float64};
                            gamma::Real=0.7)

    if d == 0
        return nothing, 0
    end
    
    as, vs = nothing, -Inf
    
    for a in actions[s.time]
        v = mdp_reward(s, image_lookup)
        
        for sp in mdp_future_states(s, a, int_lookup)
            ap, vp = mdp_select_action(sp, d-1, actions, image_lookup, int_lookup, collect_probabilities, gamma=gamma)
            v      = v + gamma*mdp_transition(sp, s, a, collect_probabilities)*vp
        end
        
        if v > vs
            as, vs = a, v
        end
    end
    
    return as, vs
end

export mdp_transition_state
"""
Transitions the MDP state to the next state given the action.
"""
function mdp_transition_state(s::MDPState, a::Collect, int_lookup::Dict{Image, <:Integer})
    sn = MDPState(a, deepcopy(s.images))
    
    # Set state to collected
    sn.images[int_lookup[a.image]] = true
    
    return sn
end

export mdp_forward_search
"""
Generate tasking policy through forward search algorithm.
"""
function mdp_forward_search(collects::Array{Collect, 1}, constraint_list, 
                            images::Array{Image,1}, collect_probabilities; 
                            horizon=0::Real, search_depth=5::Integer, gamma=0.7::Real)
    # Compute action set
    actions = mdp_construct_actions(collects, constraint_list, horizon=7200.0)
    
    # Compute Initial start
    col_init = collect(keys(actions))[findmin(collect([a.sow for a  in keys(actions)]))[2]]
    s_init   = MDPState(col_init, images)
    
    # Create lookups
    image_lookup = mdp_image_lookup(images)
    int_lookup = Dict{Image, Int64}([v => k for (k,v) in image_lookup])
    
    # Compute initial action
    a, v = mdp_select_action(s_init, search_depth, actions, image_lookup, int_lookup, collect_probabilities, gamma=gamma)
#     println("Optimal action: $a")
    s    = mdp_transition_state(s_init, a, int_lookup)
    
    plan = Union{Collect, Nothing}[a]
    
    # Perform forward search 
    while true
        a, v = mdp_select_action(s, search_depth, actions, image_lookup, int_lookup, collect_probabilities, gamma=gamma)  
        
        if a == nothing
            # Exit early if no possible action left
            break
        end
#         println("State: $s\nReward:$v\nOptimal Action: $a\n")
        
        # Update State
        push!(plan, a)

        # Check validity
        valid_transition = true
        for constraint in constraint_list            
            # Use logical and to evaulate path feasibility on graph
            valid_transition = valid_transition && constraint(s.time, a)
        end

        if valid_transition == false
            println("Invalid transition....wtf?!?!?")
        end
        
        # Transition state
        s = mdp_transition_state(s, a, int_lookup)
    end
    
    # Compute plan reward
    reward = 0
    image_list = Image[]
    
    num_pos = 0
    num_neg = 0
    
    for (i, img_state) in enumerate(s.images)
        if img_state
            img     = image_lookup[i]
            reward += img.reward
            push!(image_list, img)
            num_pos += 1
        else
            num_neg += 1
        end
    end
    
    println("Number of all images: $(length(s.images)), Collected: $num_pos, Missed: $num_neg")
    
    return plan, reward, image_list
end

# """
# Compute transitions from 

# Arguments:
# - `collects::Array{Collect, 1}` Array of collects to compute feasible transitoins for
# - `constraint_list` List of constraints 
# - `horizon::Real` Look-ahead horizon to compute possible transitions

# Returns:
# - `states::Array{Collect,1}` List of states of the decision problem
# - `transitions::Dict{Collect, Array{Collect, 1}}` List of possible transitions for each state
# """
# function mdp_compute_transitions(collects::Array{Collect, 1}, constraint_list; horizon=0::Real)

#     # Sort Collects to ensure they are in time-asecnding order
#     sort!(collects, by = x -> x.sow)

#     # Compute possible transitions
#     transitions = Dict{Collect, Array{Collect, 1}}()

#     for start_collect in collects
#         # List of valid edges/transitions for start_collect
#         transitions[start_collect] = Array{Collect, 1}[]
        
#         for end_collect in collects
#             # Decide to insert edge
#             if (start_collect == end_collect || 
#                 start_collect.image == end_collect.image ||
#                 end_collect.sow < start_collect.eow)
#                 # Skip insertion if same collect, image, or starts before the current
#                 # collection ends (no instantaneous manevers)
#                 continue
#             elseif horizon > 0 && end_collect.sow > (start_collect.eow + horizon)
#                 # Since we know collects are sorted we can break building the
#                 # transition grarph if the distance from the next to the next start
#                 # is greater than the look-ahead horizon
#                 break
#             else
#                 # Set transition valid by default
#                 valid_transition = true

#                 for constraint in constraint_list
#                     # If not valid transition break early
#                     if valid_transition == false
#                         continue
#                     end

#                     # Use logical and to evaulate path feasibility on transitions
#                     valid_transition = valid_transition && constraint(start_collect, end_collect)
#                 end

#                 if valid_transition
#                     push!(transitions[start_collect], end_collect)
#                 end
#             end
#         end
#     end

#     states = sort!(collect(keys(transitions)), by = x -> x.sow)

#     return states, transitions
# end

# function _mdp_update_state_values(Uk::Dict{Collect, <:Real}, Ukp::Dict{Collect, <:Real}, states::Array{Collect, 1}, transitions::Dict{Collect, Array{Collect, 1}})
#     # Perform initial update
#     for s in states

#         r_s = nothing # State reward
#         if length(transitions[s]) == 0
#             r_s = 0.0
#         end 

#         for a in transitions[s]
#             # Transition probability to action state is always 1.0
#             r_a = a.image.reward + Uk[a]

#             # Update optimal state value 
#             if r_s == nothing || r_a > r_s
#                 r_s = r_a
#             end

#         Ukp[s] = r_s

#         end
#     end
    
# end

# function _mdp_extract_policy(U::Dict{Collect, <:Real}, states::Array{Collect, 1}, transitions::Dict{Collect, Array{Collect, 1}})
#     policy = Dict{Collect, Union{Collect, Nothing}}()

#     for s in states
#         r_s = nothing
#         if length(transitions[s]) == 0
#             r_s = 0.0
#         end

#         a_max = nothing
#         for a in transitions[s]
#             # Transition probability to action state is always 1.0
#             r_a = a.image.reward + U[a]

#             # Update optimal state value 
#             if r_s == nothing || r_a > r_s
#                 r_s   = r_a
#                 a_max = a
#             end
#         end

#         policy[s] = a_max
#     end

#     return policy
# end

# function mdp_value_iteration(states::Array{Collect, 1}, transitions::Dict{Collect, Array{Collect, 1}}, eps=1e-6::Real)

#     # Initialize Initial value function 
#     Uk  = Dict{Collect, Float64}(s => 0.0 for s in states)
#     Ukp = Dict{Collect, Float64}(s => 0.0 for s in states)
#     k   = 1

#     # Perofrm first round of value ierations
#     _mdp_update_state_values(Uk, Ukp, states, transitions)

#     # Iterate while not-converged
#     resid = norm([Ukp[s] - Uk[s] for s in states])
#     while resid > eps
#         # Update values in UkP
#         for (k,v) in Ukp
#             Uk[k] = v
#         end

#         _mdp_update_state_values(Uk, Ukp, states, transitions)
#         resid = norm([Ukp[s] - Uk[s] for s in states])
#         k += 1

#         @debug "Finished iteration $k. Bellman residual: $resid"
#     end

#     @debug "Finished value iteration. Converged after $k iterations with Bellman residual: $resid"

#     # Extract optimal policy from state values
#     policy = _mdp_extract_policy(Ukp, states, transitions) 

#     # Extract optimal path from policy
#     max_s    = findmax(Ukp)[2]
#     opt_path = Union{Collect,Nothing}[max_s]

#     while policy[opt_path[end]] != nothing
#         push!(opt_path, policy[opt_path[end]])
#     end

#     # # Remove last element since it's going to be nothing
#     # pop!(path)

#     return Ukp, opt_path
# end

# export sp_mdp_policy
# """
# Solve for optimal collect plan using Markov Decision Process Approach

# Arguments:
# - `collects::Array{Collect,1}` Array of collects to plan the optimal tasking schedule for
# - `constraints::Array{Any, 1}` Array of constraint function which may limit feasible transitions
# - `horizon::Real` Look-ahead horizon for constructing graphs. Transition further than this apart are not considered. Not used if 0
# - `allow_repeats::Bool` Allow images to be collected multiple times over the course of a plan

# Returns:
# - `collect_policy::Array{Collect}` List of collects to take in the order which they should be taken
# """
# function sp_mdp_policy(collects::Array{Collect, 1}, constraint_list; horizon=0::Real, allow_repeats=false::Bool, eps=1e-6::Real)
#     # Compute States and Transitions
#     states, transitions = mdp_compute_transitions(collects, constraint_list, horizon=horizon)

#     # Solve graph for taskign policy
#     Ukp, path = mdp_value_iteration(states, transitions, eps)

#     reward     = findmax(Ukp)[2]
#     image_list = []
#     for col in path
#         if !(col.image in image_list)
#             push!(image_list, col.image)
#         end
#     end

#     return path, reward, image_list
# end

end # End MDP module