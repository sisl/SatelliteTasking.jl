export satellite_plan_milp

function satellite_plan_milp(problem::PlanningProblem; sat_id::Integer=1,
    contact_reward::Real=0, timeout::Real=900,
    enforce_first::Bool=false)

    # Initialize MILP problem
    milp = Model(with_optimizer(Gurobi.Optimizer, Presolve=0, TimeLimit=timeout, OutputFlag=0))

    # Initialize Variables
    @variable(milp, x[keys(problem.lt_collects)], Bin)

    # Add Objective
    @objective(milp, Max, sum(opp.reward*x[id] for (id, opp) in problem.lt_collects))

    # Enforce constraint that first opportunity in time must be taken to
    # ensure consistency with other planning methods
    if enforce_first == true
        @constraint(milp, x[1] == 1)
    end

    # Define non-repeat constraints if allowed
    if problem.solve_allow_repeats == false
        # Group opportunities by image
        for req in problem.requests
            # Get all collects for that request
            collects = collect(filter(x -> x.location.id == req.id, problem.opportunities))

            if length(collects) > 1
                @constraint(milp, sum(x[col.id] for col in collects) <= 1)
            end
        end
    end

    for (sopp_id, start_opp) in problem.lt_collects
        for (eopp_id, end_opp) in problem.lt_collects
            
            if (sopp_id == eopp_id) || (start_opp.t_start >= end_opp.t_end)
               # Conditions to skip insertion
            else
                # Check to see if each transition constraint is valid
                for constraint in problem.constraints
                    if constraint(start_opp, end_opp) == false
                        @constraint(milp, x[sopp_id] + x[eopp_id] <= 1)
                        break # Exit inner loop early
                    end
                end
            end
        end
    end

    # Solve problem
    optimize!(milp)

    # Extract plan and reward from MILP
    plan = Opportunity[]
    reward = 0.0

    for idx in keys(problem.lt_collects)
        if value(x[idx]) != 0.0
            push!(plan, problem.lt_collects[idx])
        end
    end

    # Sort plan in ascending order to be viable
    sort!(plan, by = x -> x.t_start)

    reward = objective_value(milp)

    return plan, reward
end