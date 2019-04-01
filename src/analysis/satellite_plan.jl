__precompile__(true)
module SatellitePlan

# Julia Imports
using Statistics
using SatelliteDynamics.Time: Epoch

# Package Imports
using SatelliteTasking.DataStructures: Orbit, Image, Opportunity, Collect

# Evaluate plan

export sp_check_feasibility
"""
Independently check a collect plan (array of collects in order) against a set of
constraint functions to determine whether all transitions are feasible.

Arguments:
- `collects::Array{<:Any, 1}` Collect plan to validate
- `colelcts::Array{Function, 1}` Array of constraint functions to validate plan against

Returns:
- `valid::Bool` Returns `true` if collect plan is valid, returns `false` otherwise.
"""
function sp_check_feasibility(collects::Array{<:Any,1}, constraint_list::Array{Function,1})
    # Independently check that all constraints are still met
    for i in 1:length(collects)-1
        col_start = collects[i]
        col_end   = collects[i+1]

        # Set transition valid by default
        valid_transition = true

        for constraint in constraint_list
            # Use logical and to evaulate path feasibility on graph
            valid_transition = valid_transition && constraint(col_start, col_end)
        end

        if !valid_transition
            println("Collect Failed:\n Start: $col_start\n End: $col_end\n")
            return false
        end
    end
    
    return true
end

export sp_evaluate_plan
"""
Evaluate performance of plan against against a set of opportunities.

Arguments:
- `plan::Array{<:Any, 1}` Tasking plan to evaluate feasibility and reward for
- `opportunities::Array{Opportunity, 1}` List of opportunities to evaluate collections against.

Returns:
- `realized_reward` Total reward achieved by the plan
- `feasible_collects` List of collects which are feasible and taken
- `infeasible_collects` List of collects which are no longer possible due to uncertainties
"""
function sp_evaluate_plan(plan::Array{<:Any,1}, opportunities::Array{Opportunity,1})
    realized_reward     = 0.0
    feasible_collects   = Collect[]
    infeasible_collects = Collect[]
    
    unique_images = Image[]
    
    for collect in plan        
        if length(filter(opp -> opp.image == collect.image 
                            && opp.sow <= collect.sow <= opp.eow
                            && opp.sow <= collect.eow <= opp.eow, opportunities)) >= 1
            if !(collect.image in unique_images) 
                realized_reward += collect.image.reward
                push!(feasible_collects, collect)
                push!(unique_images, collect.image)
            end
        else
            push!(infeasible_collects, collect)
        end
    end
    
    return realized_reward, feasible_collects, infeasible_collects
end


export sp_compute_collect
"""
Evaluate performance of plan against against a set of opportunities.

Arguments:
- `plan::Array{<:Any, 1}` Tasking plan to evaluate feasibility and reward for
- `opportunities::Array{Array{Opportunity,1}, 1}` Set of opportunities to evaluate realized reward against

Returns:
- `realized_reward` Reward realized by following the collect plan for each set of opportunities
- `mean_reward` Mean reward realized by following the tasking plan 
- `sdev_reward` Standard deviation of the reward realized when following the tasking plan
"""
function sp_compute_collect(plan::Array{<:Any,1}, perturbed_opportunities::Array{Array{Opportunity,1},1})
    realized_rewards = Float64[]
    
    for opp_set in perturbed_opportunities
        realized_reward, fcollects, ifcollects = sp_evaluate_plan(plan, opp_set)
        
        push!(realized_rewards, realized_reward)
    end
    
    mean_reward  = mean(realized_rewards)
    sdev_reward = std(realized_rewards)
    
    return realized_rewards, mean_reward, sdev_reward 
end

end # End CollectUncertainty module
