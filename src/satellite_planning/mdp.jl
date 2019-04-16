__precompile__(true)
module MDP

# Julia Packages
using LinearAlgebra

using SatelliteDynamics.Time: Epoch
using SatelliteTasking.DataStructures: Image, Opportunity


export MDPState
"""
State 

Attributes:
- `time::Epoch` Current planning time
- `locations::Array{Bool,1}` Boolean array of whether a specific image has been opportunityed. Index indicates image.
"""
mutable struct MDPState
    time::Opportunity
    locations::Array{Bool, 1}
end

function MDPState(epc::Opportunity, image_list::Array{Image,1})    
    return MDPState(epc, falses(length(image_list)))
end

export mdp_construct_actions
"""
Construct the possible action space for all states. The actions are all dependent
on time alone, (not which locations have already been opportunityed).

Arguments:
- `opportunities::Array{Opportunity, 1}` List of possible actions to take 
- `constraint_list::Array{Opportunity, 1}` List of constraintt functions that limit feasibility in transition.
- `horizon::Real` Horizon to consider building actions. If 0 and infinite horizon will be used. Other times reduce problem complexity.

Returns:
- `actions::Dict{Epoch, Array{Opportunity, 1}}` List of possible actions to take 
"""
function mdp_construct_actions(opportunities::Array{Opportunity, 1}, constraint_list::Array{Function, 1}; horizon=0::Real)
    
    # Sort Opportunitys to ensure they are in time-asecnding order
    sort!(opportunities, by = x -> x.sow)

    # Initialize validate transition graph
    actions = Dict{Opportunity, Array{Opportunity, 1}}()

    for start_opportunity in opportunities
        # List of valid edges/transitions for start_opportunity
        actions[start_opportunity] = Array{Opportunity, 1}[]
        
        for end_opportunity in opportunities
            # Decide to insert edge
            if (start_opportunity == end_opportunity || 
                start_opportunity.location == end_opportunity.location ||
                end_opportunity.sow < start_opportunity.eow)
                # Skip insertion if same opportunity, image, or starts before the current
                # opportunity ends (no instantaneous manevers)
                continue
            elseif horizon > 0 && end_opportunity.sow > (start_opportunity.eow + horizon)
                # Since we know opportunities are sorted we can break building the
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
                    valid_transition = valid_transition && constraint(start_opportunity, end_opportunity)
                end

                if valid_transition
                    push!(actions[start_opportunity], end_opportunity)
                end
            end
        end
    end
    
    return actions
end

export mdp_image_lookup
"""
Given a list of locations create the dictionary lookup of indices 

Arguments:
- `image_list::Array{Image, 1}` List of locations (should match MDPState list)

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

export compute_opportunity_probability
"""
Compute the probability that a opportunity will be possible given the apparent
distribution of opportunity probabilities.
"""
function compute_opportunity_probability(opportunities::Array{Opportunity}, observed_opporttunities::Array{Array{Opportunity,1},1})
    opportunity_probability = Dict{Opportunity, Float64}()
    for opp in opportunities
        num_present, num_missing = 0, 0
        for opportunity_list in observed_opporttunities
            if length(filter(x -> x.location == opp.location && x.sow <= opp.sow && x.eow >= opp.eow, opportunity_list)) == 1
                num_present += 1
            else
                num_missing += 1
            end
        end
        
        # opportunity_probability[opp] = num_present/(num_missing+num_present)
        opportunity_probability[opp] = 1
    end
    
    return opportunity_probability
end

export mdp_reward
"""
Compute reward of current MDP state

Arguments:
- `s::MDPState` State of Markov decision process
"""
function mdp_reward(s::MDPState, image_lookup::Dict{<:Integer, Image})
    r = 0.0
    for (i, img) in enumerate(s.locations)
        if img
            r += image_lookup[i].reward
        else
        # if !img
            r -= image_lookup[i].reward
        end
    end
    
    return r
end

export mdp_future_states
"""
Compute possible future takes the system may assume given the action.
"""
function mdp_future_states(state::MDPState, action::Opportunity, int_lookup::Dict{Image, <:Integer})
    failure      = MDPState(state.time, copy(state.locations))
    failure.time = action
    
    success = MDPState(failure.time, copy(failure.locations))
    success.locations[int_lookup[action.location]] = true
    
    return MDPState[success, failure]
end

export mdp_transition
"""
Return transition probabilities given the set of currrent states, action, and possible future states.
"""
function mdp_transition(sp::MDPState, s::MDPState, a::Opportunity, opportunity_probabilities::Dict{Opportunity, Float64})
    p = opportunity_probabilities[a]
    
    return p
end

export mdp_select_action
"""
MDP forward search algorithm to select the next decision.
"""
function mdp_select_action(s::MDPState, d::Integer, actions::Dict{Opportunity, Array{Opportunity, 1}}, 
                            image_lookup::Dict{<:Integer, Image}, int_lookup::Dict{Image, <:Integer},
                            opportunity_probabilities::Dict{Opportunity, Float64};
                            gamma::Real=0.7)

    if d == 0
        return nothing, 0
    end
    
    as, vs = nothing, -Inf
    
    for a in actions[s.time]
        v = mdp_reward(s, image_lookup)
        
        for sp in mdp_future_states(s, a, int_lookup)
            ap, vp = mdp_select_action(sp, d-1, actions, image_lookup, int_lookup, opportunity_probabilities, gamma=gamma)
            v      = v + gamma*mdp_transition(sp, s, a, opportunity_probabilities)*vp
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
function mdp_transition_state(s::MDPState, a::Opportunity, int_lookup::Dict{Image, <:Integer})
    sn = MDPState(a, deepcopy(s.locations))
    
    # Set state to opportunityed
    sn.locations[int_lookup[a.location]] = true
    
    return sn
end

export mdp_forward_search
"""
Generate tasking policy through forward search algorithm.
"""
function mdp_forward_search(opportunities::Array{Opportunity, 1}, constraint_list, 
                            locations::Array{Image,1}, opportunity_probabilities; 
                            horizon=0::Real, search_depth=5::Integer, gamma=0.7::Real)
    # Compute action set
    actions = mdp_construct_actions(opportunities, constraint_list, horizon=7200.0)
    
    # Compute Initial start
    opp_init = collect(keys(actions))[findmin(collect([a.sow for a  in keys(actions)]))[2]]
    s_init   = MDPState(opp_init, locations)
    
    # Create lookups
    image_lookup = mdp_image_lookup(locations)
    int_lookup = Dict{Image, Int64}([v => k for (k,v) in image_lookup])
    
    # Compute initial action
    a, v = mdp_select_action(s_init, search_depth, actions, image_lookup, int_lookup, opportunity_probabilities, gamma=gamma)
#     println("Optimal action: $a")
    s    = mdp_transition_state(s_init, a, int_lookup)
    
    plan = Union{Opportunity, Nothing}[a]
    
    # Perform forward search 
    while true
        a, v = mdp_select_action(s, search_depth, actions, image_lookup, int_lookup, opportunity_probabilities, gamma=gamma)  
        
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
    
    for (i, img_state) in enumerate(s.locations)
        if img_state
            img     = image_lookup[i]
            reward += img.reward
            push!(image_list, img)
            num_pos += 1
        else
            num_neg += 1
        end
    end
    
    println("Number of all locations: $(length(s.locations)), Opportunityed: $num_pos, Missed: $num_neg")
    
    return plan, reward, image_list
end

# """
# Compute transitions from 

# Arguments:
# - `opportunities::Array{Opportunity, 1}` Array of opportunities to compute feasible transitoins for
# - `constraint_list` List of constraints 
# - `horizon::Real` Look-ahead horizon to compute possible transitions

# Returns:
# - `states::Array{Opportunity,1}` List of states of the decision problem
# - `transitions::Dict{Opportunity, Array{Opportunity, 1}}` List of possible transitions for each state
# """
# function mdp_compute_transitions(opportunities::Array{Opportunity, 1}, constraint_list; horizon=0::Real)

#     # Sort Opportunitys to ensure they are in time-asecnding order
#     sort!(opportunities, by = x -> x.sow)

#     # Compute possible transitions
#     transitions = Dict{Opportunity, Array{Opportunity, 1}}()

#     for start_opportunity in opportunities
#         # List of valid edges/transitions for start_opportunity
#         transitions[start_opportunity] = Array{Opportunity, 1}[]
        
#         for end_opportunity in opportunities
#             # Decide to insert edge
#             if (start_opportunity == end_opportunity || 
#                 start_opportunity.location == end_opportunity.location ||
#                 end_opportunity.sow < start_opportunity.eow)
#                 # Skip insertion if same opportunity, image, or starts before the current
#                 # opportunity ends (no instantaneous manevers)
#                 continue
#             elseif horizon > 0 && end_opportunity.sow > (start_opportunity.eow + horizon)
#                 # Since we know opportunities are sorted we can break building the
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
#                     valid_transition = valid_transition && constraint(start_opportunity, end_opportunity)
#                 end

#                 if valid_transition
#                     push!(transitions[start_opportunity], end_opportunity)
#                 end
#             end
#         end
#     end

#     states = sort!(opportunity(keys(transitions)), by = x -> x.sow)

#     return states, transitions
# end

# function _mdp_update_state_values(Uk::Dict{Opportunity, <:Real}, Ukp::Dict{Opportunity, <:Real}, states::Array{Opportunity, 1}, transitions::Dict{Opportunity, Array{Opportunity, 1}})
#     # Perform initial update
#     for s in states

#         r_s = nothing # State reward
#         if length(transitions[s]) == 0
#             r_s = 0.0
#         end 

#         for a in transitions[s]
#             # Transition probability to action state is always 1.0
#             r_a = a.location.reward + Uk[a]

#             # Update optimal state value 
#             if r_s == nothing || r_a > r_s
#                 r_s = r_a
#             end

#         Ukp[s] = r_s

#         end
#     end
    
# end

# function _mdp_extract_policy(U::Dict{Opportunity, <:Real}, states::Array{Opportunity, 1}, transitions::Dict{Opportunity, Array{Opportunity, 1}})
#     policy = Dict{Opportunity, Union{Opportunity, Nothing}}()

#     for s in states
#         r_s = nothing
#         if length(transitions[s]) == 0
#             r_s = 0.0
#         end

#         a_max = nothing
#         for a in transitions[s]
#             # Transition probability to action state is always 1.0
#             r_a = a.location.reward + U[a]

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

# function mdp_value_iteration(states::Array{Opportunity, 1}, transitions::Dict{Opportunity, Array{Opportunity, 1}}, eps=1e-6::Real)

#     # Initialize Initial value function 
#     Uk  = Dict{Opportunity, Float64}(s => 0.0 for s in states)
#     Ukp = Dict{Opportunity, Float64}(s => 0.0 for s in states)
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
#     opt_path = Union{Opportunity,Nothing}[max_s]

#     while policy[opt_path[end]] != nothing
#         push!(opt_path, policy[opt_path[end]])
#     end

#     # # Remove last element since it's going to be nothing
#     # pop!(path)

#     return Ukp, opt_path
# end

# export sp_mdp_policy
# """
# Solve for optimal opportunity plan using Markov Decision Process Approach

# Arguments:
# - `opportunities::Array{Opportunity,1}` Array of opportunities to plan the optimal tasking schedule for
# - `constraints::Array{Any, 1}` Array of constraint function which may limit feasible transitions
# - `horizon::Real` Look-ahead horizon for constructing graphs. Transition further than this apart are not considered. Not used if 0
# - `allow_repeats::Bool` Allow locations to be opportunityed multiple times over the course of a plan

# Returns:
# - `opportunity_policy::Array{Opportunity}` List of opportunities to take in the order which they should be taken
# """
# function sp_mdp_policy(opportunities::Array{Opportunity, 1}, constraint_list; horizon=0::Real, allow_repeats=false::Bool, eps=1e-6::Real)
#     # Compute States and Transitions
#     states, transitions = mdp_compute_transitions(opportunities, constraint_list, horizon=horizon)

#     # Solve graph for taskign policy
#     Ukp, path = mdp_value_iteration(states, transitions, eps)

#     reward     = findmax(Ukp)[2]
#     image_list = []
#     for opp in path
#         if !(opp.location in image_list)
#             push!(image_list, opp.location)
#         end
#     end

#     return path, reward, image_list
# end

end # End MDP module