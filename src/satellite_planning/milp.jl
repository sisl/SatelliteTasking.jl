__precompile__(true)
module MILP

# Package imports
using JuMP
using Gurobi

# SatelliteTasking imports
using SatelliteTasking.DataStructures: Image, Opportunity, Collect

export sp_construct_milp
"""
Solve for optimal collect plan using mixed-integer linear programming.

Arguments:
- `collects::Array{Collect,1}` Array of collects to plan the optimal tasking schedule for
- `constraints::Array{Any, 1}` Array of constraint function which may limit feasible transitions
- `horizon::Real` Look-ahead horizon for constructing constraints. Transitions further than this apart are not considered. Not used if 0
- `allow_repeats::Bool` Allow images to be collected multiple times over the course of a plan

Returns:
- `milp` JuMP Optimization problem
"""
function sp_construct_milp(collects::Array{Collect, 1}, constraint_list::Array{Function, 1}; horizon=0::Real, allow_repeats=false::Bool)
    # Initialize MILP problem
    milp = Model(solver=GurobiSolver())

    # Sort Collects to ensure they are in time-asecnding order
    sort!(mcollects, by = x -> x.sow)

    # constraint_list = Function[constraint_agility_single_axis]
    constraint_list = Function[constraint_agility_single_axis]

    horizon = 0
    allow_repeats = false

    # Initialize Variables
    @variable(milp, x[1:length(mcollects)], Bin)

    # Add Objective
    @objective(milp, Max, sum(col.image.reward*x[i] for (i,col) in enumerate(mcollects)))

    # Define non-repetition constraints if necessary
    if allow_repeats == false
        # Group collects by image
        image_collects = group_image_collects(mcollects) # Group collects by image
        
        # Add constraints to limit one collect per image
        for img in keys(image_collects)
            # Only add constraints for when there is more than one possible collect
            if length(image_collects[img]) > 1
                @constraint(milp, sum(x[i] for i in collect(e[1] for e in image_collects[img])) <= 1)
            end
        end
    end

    # Add satellite model-derived constraints
    for (i, col_start) in enumerate(mcollects)
        for (j, col_end) in enumerate(mcollects)
            if horizon > 0 && col_end.sow > (col_start.eow + horizon)
                # Condition to exit early is only considering transitions within a certain horizon may be invalid
                break
            end

            if col_start == col_end || col_start.image == col_end.image
                # Skip if the same opportunity, or same imaage
                continue
            else
                # Transition is default valid
                valid_transition = true
                
                # Can't go backwards in time
                if col_end.sow < col_start.eow || col_start.sow < col_end.eow
                    valid_transition = false 
                else
                    # Only evaluate transitions if time constraint doesn't matter
                    for constraint in constraint_list
                        # Use logical and to evaulate path feasibility
                        # Because this is a binary comparison (independent of time), and constraints are nominally
                        # evaluated with start dependent on the end, we only add a constraint if neither transition is valid
                        # Otherwise the problem would be over constrained just due to the final opportunity not
                        # being able to take images before the current
                        valid_transition = valid_transition && (constraint(col_start, col_end) || constraint(col_end, opp_start))
                    end
                end

                if !valid_transition
                    # Boolean logic constraint permitting taking start opportunity, but
                    @constraint(milp, x[i] + x[j] <= 1)
                end
            end
        end
    end
end

#     # Define binary decision variables on opportunities, and lookups
#     model_vars = {}
#     idx_lookup = {}
#     opp_lookup = {}
#     for i, opp in enumerate(opportunities):
#         model_vars[i]   = model.addVar(f'x_{i:d}', vtype='BINARY')
#         idx_lookup[i]   = opp
#         opp_lookup[opp] = i

#     # Add non-repeat constraints
#     if not repeat_images:
#         # Construct dict of images and all associated opportunities
#         image_dict = get_image_opportunity_lookup(opportunities)
    
#         # Add single-image constraint
#         for i, image in enumerate(image_dict.keys()):
#             model.addCons(_scip.quicksum([model_vars[opp_lookup[opp]] for opp in image_dict[image]]) <= 1, name=f'i_{i:d}')

#     # Add satellite constraints
#     for i, opp_start in enumerate(opportunities):
#         for j, opp_end in enumerate(opportunities):
#             if horizon and opp_end.aos > (opp_start.los + horizon):
#                 # Condition to exit early is only considering transitions within a certain horizon may be invalid
#                 break

#             if opp_start == opp_end or opp_start.image == opp_end.image:
#                 # Skip if the same opportunity, or same imaage
#                 continue
#             else:
#                 # Transition is default valid
#                 valid_transition = True

#                 for constraint in constraint_list:
#                     # Use logical and to evaulate path feasibility
#                     # Because this is a binary comparison (independent of time), and constraints are nominally
#                     # evaluated with start dependt on the end, we only add a constraint if neither transition is valid
#                     # Otherwise the problem would be over constrained just due to the final opportunity not
#                     # being able to take images before the current
#                     valid_transition = valid_transition and (constraint(opp_start, opp_end) or constraint(opp_end, opp_start))

#                 if not valid_transition:
#                     # Boolean logic constraint permitting taking start opportunity, but
#                     # if the start opportunity is taken the end cannot be
#                     model.addCons(_scip.quicksum([model_vars[opp_lookup[opp_start]], model_vars[opp_lookup[opp_end]]]) <= 1, name=f'b_{i:d}_{j:d}')

    

#     # Construct objective
#     objective = _scip.quicksum([opp.image.reward * model_vars[i] for i, opp in enumerate(opportunities)])
#     model.setObjective(objective, 'maximize')

#     return model, idx_lookup, opp_lookup

export sp_solve_milp
"""
Solve a JuMP MILP tasking problem and return the optimal policy, total reward, and images collected.
"""
function sp_solve_milp(milp)
    solve(milp)



    # Extract policy information
end

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
    # Construct graph
    milp  = sp_construct_milp(collects, constraint_list, horizon, allow_repeats)

    # Solve graph for taskign policy
    path, reward, image_list = sp_solve_milp(milp)

    return path, reward, image_list
end

end # End MILP module