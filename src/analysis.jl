# Exports
export analyze_plan
export plot_resources

function analyze_plan(problem::PlanningProblem, plan::Array{<:Opportunity, 1})
    
    # Evaluate any constraint violations
    feasible = true
    for i in 1:length(plan)-1
        col_start = plan[i]
        col_end   = plan[i+1]

        # Only check validity of opportunity transitions
        if (typeof(col_start) == Collect ||
            typeof(col_start) == Contact) &&
            (typeof(col_end) == Collect ||
            typeof(col_end) == Contact)
            # Set transition valid by default
            valid = true

            for constraint in problem.constraints
                # Use logical and to evaulate path feasibility on graph
                if !constraint(col_start, col_end)
                    println("Collect Failed:\n Start: $col_start\n End: $col_end\n")
                    feasible = false
                    break
                end
            end
        end
    end

    unique_requests = Request[]
    duplicate_requests = Request[]

    n_contacts = 0
    n_requests = 0
    n_dup_requests = 0

    reward = 0.0

    for opp in plan
        if typeof(opp) == Collect
            n_requests += 1
            if !(opp.location in unique_requests)
                reward += opp.reward
                push!(unique_requests, opp.location)
            else
                push!(duplicate_requests, opp.location)
                n_dup_requests += 1
            end
        elseif typeof(opp) == Contact
            n_contacts += 1
        end
    end
    n_unique_reqests = length(unique_requests)
    
    return feasible, reward, n_unique_reqests, n_dup_requests, n_requests, n_contacts
end

function plot_resources(states::Array{MDPState, 1}; width::Real=25, height::Real=12.5, hsep::Real=2.0, vsep::Real=2.0, xmax::Real=24)
    # Plot history of power

    # Initialize array
    epochs = Epoch[]
    times  = Float64[]
    power  = Float64[]
    data   = Float64[]

    # Aggregate data
    for state in states
        push!(epochs, state.time)
        push!(times, (state.time - epochs[1])/3600.0)
        push!(power, state.power)
        push!(data, state.data)
    end

    hsepstr = @sprintf "%.2f" hsep
    vsepstr = @sprintf "%.2f" vsep
    wsepstr = @sprintf "%.2f" width
    hsepstr = @sprintf "%.2f" height

    g = GroupPlot(1, 2, groupStyle = "horizontal sep = $(hsepstr)cm, vertical sep = $(vsepstr)cm")
    push!(g, Axis(Plots.Linear(times, power), 
        xlabel="Elapsed Time [Hrs]", 
        ylabel="Power Fraction", 
        title="Satellite Power History",
        xmin=0.0,
        xmax=xmax,
        ymin=0.0,
        ymax=1.0,
        style = "width=$(wsepstr)cm,height=$(hsepstr)cm"
    ))

    push!(g, Axis(Plots.Linear(times, data), 
        xlabel="Elapsed Time [Hrs]", 
        ylabel="Data Storage", 
        title="Satellite Data History",
        xmin=0.0,
        xmax=xmax,
        ymin=0.0,
        ymax=1.0,
        style = "width=$(wsepstr)cm,height=$(hsepstr)cm"
    ))
end