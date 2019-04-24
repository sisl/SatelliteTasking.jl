__precompile__(true)
module MILP

# Package imports
using JuMP
using Gurobi

# SatelliteTasking imports
using SatelliteTasking.DataStructures: Image, Opportunity
using SatelliteTasking.Collection: group_image_opportunities

export sp_milp_policy
"""
Solve for optimal opportunity plan using mixed-integer linear programming.

Arguments:
- `opportunities::Array{Opportunity,1}` Array of opportunities to plan the optimal tasking schedule for
- `constraints::Array{Any, 1}` Array of constraint function which may limit feasible transitions
- `horizon::Real` Look-ahead horizon for constructing constraints. Transitions further than this apart are not considered. Not used if 0
- `allow_repeats::Bool` Allow images to be opportunityed multiple times over the course of a plan

Returns:
- `opportunity_policy::Array{Opportunity}` List of opportunities to take in the order which they should be taken
"""
function sp_milp_policy(opportunities::Array{Opportunity, 1}, constraint_list; horizon::Real=0, allow_repeats::Bool=false, timeout::Real=900)
    
    # Sort Opportunitys to ensure they are in time-asecnding order
    opportunities = sort!(opportunities, by = x -> x.sow)

    # Initialize MILP problem
    # milp = Model(solver=GurobiSolver())
    # milp = Model(with_optimizer(Gurobi.Optimizer, Presolve=0, Heuristics=0.0, TimeLimit=900))
    milp = Model(with_optimizer(Gurobi.Optimizer, TimeLimit=900))

    # Initialize Variables
    @variable(milp, x[1:length(opportunities)], Bin)

    # Add Objective
    @objective(milp, Max, sum(opp.location.reward*x[i] for (i,opp) in enumerate(opportunities)))

    # Define non-repetition constraints if necessary
    if allow_repeats == false
        # Group opportunities by image
        image_opportunities = group_image_opportunities(opportunities) # Group opportunities by image
        
        # Add constraints to limit one opportunity per image
        for img in keys(image_opportunities)
            # Only add constraints for when there is more than one possible opportunity
            if length(image_opportunities[img]) > 1
                @constraint(milp, sum(x[i] for i in collect(e[1] for e in image_opportunities[img])) <= 1)
            end
        end
    end

    # Add satellite model-derived constraints
    for i in 1:length(opportunities)
        for j in  i:length(opportunities)
            # Since all constraints are reciprocal they only need to be checked in one direction
            opp_start = opportunities[i]
            opp_end   = opportunities[j]

            # Skipp adding constraint if different satellites
            if opp_start.orbit.id != opp_end.orbit.id
                continue
            end
            
            # Skip adding constraints if a planning horizon is being used
            if horizon > 0 && opp_end.sow > (opp_start.eow + horizon)
                # Condition to exit early is only considering transitions within a certain horizon may be invalid
                continue
            end

            if opp_start == opp_end || opp_start.location == opp_end.location j < i
                # Skip if the same opportunity, or same image because this is already covered
                continue
            else
                # Transition is default valid
                valid_transition = true
                
                # Only evaluate transitions if time constraint doesn't matter
                for cons in constraint_list
                    # Use logical and to evaulate path feasibility
                    # Because this is a binary comparison (independent of time), and constraints are nominally
                    # evaluated with start dependent on the end, we only add a constraint if neither transition is valid
                    # Otherwise the problem would be over constrained just due to the final opportunity not
                    # being able to take images before the current
                    valid_transition = valid_transition && (cons(opp_start, opp_end) || cons(opp_end, opp_start))
                    
                    # Add constraint as soon as invalid to short-circuit additional evaluations
                    if !valid_transition
                        # Boolean logic constraint permitting taking start opportunity, but
                        @constraint(milp, x[i] + x[j] <= 1)
                        # println("Invalid due to agility")
                        continue
                    end
                end
            end
        end
    end

    # Solve Problem
    optimize!(milp)

    # Extract Optimal Path, Reward, and Images
    opportunity_idx    = []
    opportunities_taken = Opportunity[]
    for i in 1:length(opportunities)
        if value(x[i]) != 0.0 # If opportunity taken
            push!(opportunity_idx, i)
            push!(opportunities_taken, opportunities[i])
        end
    end

    # Get images taken
    image_list = [opportunities[opp_idx].location for opp_idx in opportunity_idx]

    # Get Path
    path = sort!(opportunities_taken, by = o -> o.sow)

    # Get reward
    reward = objective_value(milp)

    return path, reward, image_list
end

end # End MILP module