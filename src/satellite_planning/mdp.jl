__precompile__(true)
module MDP

# Julia Packages
using LinearAlgebra

using SatelliteTasking.DataStructures: Image, Opportunity, Collect

"""
Compute transitions from 

Arguments:
- `collects::Array{Collect, 1}` Array of collects to compute feasible transitoins for
- `constraint_list` List of constraints 
- `horizon::Real` Look-ahead horizon to compute possible transitions

Returns:
- `states::Array{Collect,1}` List of states of the decision problem
- `transitions::Dict{Collect, Array{Collect, 1}}` List of possible transitions for each state
"""
function mdp_compute_transitions(collects::Array{Collect, 1}, constraint_list; horizon=0::Real)

    # Sort Collects to ensure they are in time-asecnding order
    sort!(collects, by = x -> x.sow)

    # Compute possible transitions
    transitions = Dict{Collect, Array{Collect, 1}}()

    for start_collect in collects
        # List of valid edges/transitions for start_collect
        transitions[start_collect] = Array{Collect, 1}[]
        
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

                    # Use logical and to evaulate path feasibility on transitions
                    valid_transition = valid_transition && constraint(start_collect, end_collect)
                end

                if valid_transition
                    push!(transitions[start_collect], end_collect)
                end
            end
        end
    end

    states = sort!(collect(keys(transitions)), by = x -> x.sow)

    return states, transitions
end

function _mdp_update_state_values(Uk::Dict{Collect, <:Real}, Ukp::Dict{Collect, <:Real}, states::Array{Collect, 1}, transitions::Dict{Collect, Array{Collect, 1}})
    # Perform initial update
    for s in states

        r_s = nothing # State reward
        if length(transitions[s]) == 0
            r_s = 0.0
        end 

        for a in transitions[s]
            # Transition probability to action state is always 1.0
            r_a = a.image.reward + Uk[a]

            # Update optimal state value 
            if r_s == nothing || r_a > r_s
                r_s = r_a
            end

        Ukp[s] = r_s

        end
    end
    
end

function _mdp_extract_policy(U::Dict{Collect, <:Real}, states::Array{Collect, 1}, transitions::Dict{Collect, Array{Collect, 1}})
    policy = Dict{Collect, Union{Collect, Nothing}}()

    for s in states
        r_s = nothing
        if length(transitions[s]) == 0
            r_s = 0.0
        end

        a_max = nothing
        for a in transitions[s]
            # Transition probability to action state is always 1.0
            r_a = a.image.reward + U[a]

            # Update optimal state value 
            if r_s == nothing || r_a > r_s
                r_s   = r_a
                a_max = a
            end
        end

        policy[s] = a_max
    end

    return policy
end

function mdp_value_iteration(states::Array{Collect, 1}, transitions::Dict{Collect, Array{Collect, 1}}, eps=1e-6::Real)

    # Initialize Initial value function 
    Uk  = Dict{Collect, Float64}(s => 0.0 for s in states)
    Ukp = Dict{Collect, Float64}(s => 0.0 for s in states)
    k   = 1

    # Perofrm first round of value ierations
    _mdp_update_state_values(Uk, Ukp, states, transitions)

    # Iterate while not-converged
    resid = norm([Ukp[s] - Uk[s] for s in states])
    while resid > eps
        # Update values in UkP
        for (k,v) in Ukp
            Uk[k] = v
        end

        _mdp_update_state_values(Uk, Ukp, states, transitions)
        resid = norm([Ukp[s] - Uk[s] for s in states])
        k += 1

        @debug "Finished iteration $k. Bellman residual: $resid"
    end

    @debug "Finished value iteration. Converged after $k iterations with Bellman residual: $resid"

    # Extract optimal policy from state values
    policy = _mdp_extract_policy(Ukp, states, transitions) 

    # Extract optimal path from policy
    max_s    = findmax(Ukp)[2]
    opt_path = Union{Collect,Nothing}[max_s]

    while policy[opt_path[end]] != nothing
        push!(opt_path, policy[opt_path[end]])
    end

    # # Remove last element since it's going to be nothing
    # pop!(path)

    return Ukp, opt_path
end

export sp_mdp_policy
"""
Solve for optimal collect plan using Markov Decision Process Approach

Arguments:
- `collects::Array{Collect,1}` Array of collects to plan the optimal tasking schedule for
- `constraints::Array{Any, 1}` Array of constraint function which may limit feasible transitions
- `horizon::Real` Look-ahead horizon for constructing graphs. Transition further than this apart are not considered. Not used if 0
- `allow_repeats::Bool` Allow images to be collected multiple times over the course of a plan

Returns:
- `collect_policy::Array{Collect}` List of collects to take in the order which they should be taken
"""
function sp_mdp_policy(collects::Array{Collect, 1}, constraint_list; horizon=0::Real, allow_repeats=false::Bool, eps=1e-6::Real)
    # Compute States and Transitions
    states, transitions = mdp_compute_transitions(collects, constraint_list, horizon=horizon)

    # Solve graph for taskign policy
    Ukp, path = mdp_value_iteration(states, transitions, eps)

    reward     = findmax(Ukp)[2]
    image_list = []
    for col in path
        if !(col.image in image_list)
            push!(image_list, col.image)
        end
    end

    return path, reward, image_list
end

end # End MDP module