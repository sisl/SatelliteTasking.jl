__precompile__(true)
module MILP

# Package imports
using JuMP
using Gurobi

# SatelliteTasking imports
using SatelliteTasking.DataStructures: Image, Opportunity, Collect
using SatelliteTasking.Collects: group_image_collects

export sp_milp_policy
"""
Solve for optimal collect plan using mixed-integer linear programming.

Arguments:
- `collects::Array{Collect,1}` Array of collects to plan the optimal tasking schedule for
- `constraints::Array{Any, 1}` Array of constraint function which may limit feasible transitions
- `horizon::Real` Look-ahead horizon for constructing constraints. Transitions further than this apart are not considered. Not used if 0
- `allow_repeats::Bool` Allow images to be collected multiple times over the course of a plan

Returns:
- `collect_policy::Array{Collect}` List of collects to take in the order which they should be taken
"""
function sp_milp_policy(collects::Array{Collect, 1}, constraint_list; horizon=0::Real, allow_repeats=false::Bool)
    
    # Sort Collects to ensure they are in time-asecnding order
    collects = sort!(collects, by = x -> x.sow)

    # Initialize MILP problem
    milp = Model(solver=GurobiSolver(Presolve=0, Heuristics=0.0))

    # Initialize Variables
    @variable(milp, x[1:length(collects)], Bin)

    # Add Objective
    @objective(milp, Max, sum(col.image.reward*x[i] for (i,col) in enumerate(collects)))

    # Define non-repetition constraints if necessary
    if allow_repeats == false
        # Group collects by image
        image_collects = group_image_collects(collects) # Group collects by image
        
        # Add constraints to limit one collect per image
        for img in keys(image_collects)
            # Only add constraints for when there is more than one possible collect
            if length(image_collects[img]) > 1
                @constraint(milp, sum(x[i] for i in collect(e[1] for e in image_collects[img])) <= 1)
            end
        end
    end

    # Add satellite model-derived constraints
    for i in 1:length(collects)
        for j in  i:length(collects)
            # Since all constraints are reciprocal they only need to be checked in one direction
            col_start = collects[i]
            col_end   = collects[j]
            
            # Skip adding constraints if a planning horizon is being used
            if horizon > 0 && col_end.sow > (col_start.eow + horizon)
                # Condition to exit early is only considering transitions within a certain horizon may be invalid
                continue
            end

            if col_start == col_end || col_start.image == col_end.image || col_start.opportunity == col_end.opportunity || j < i
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
                    valid_transition = valid_transition && (cons(col_start, col_end) || cons(col_end, col_start))
                    
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
    status = solve(milp)

    # Extract Optimal Path, Reward, and Images
    collects_taken = []
    for i in 1:length(collects)
        if getvalue(x[i]) != 0.0 # If opportunity taken
            push!(collects_taken, i)
        end
    end

    # Get images taken
    image_list = [col.image for col in collects_taken]

    # Get Path
    path = sort!(collects_taken, by = x -> x.sow)

    # Get reward
    reward = getobjectivevalue(milp)

    return path, reward, image_list
end

end # End MILP module