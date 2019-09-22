# Exports

export sample_requests_uniform
export sample_requests_file
export SolveSettings
export sat_sim
export mass_sim
export sat_resource_sim
export mass_resource_sim

function sample_requests_uniform(num_samples::Int=0;lon_min::Real=-180, lon_max::Real=180,
                                lat_min::Real=-90, lat_max::Real=90,
                                look_angle_min::Real=10, look_angle_max::Real=55,
                                collect_duration_min::Real=30, collect_duration_max::Real=30,
                                reward_min::Real=1.0, reward_max::Real=1.0)

    if num_samples <= 0
        throw(ErrorException("Number of samples must be positive."))
    end

    if reward_min <= 0
        throw(ErrorException("Minimum reward must be positive."))
    end

    if reward_max < reward_min
        throw(ErrorException("Maximum reward must be greater or equal to minimum reward."))
    end

    if collect_duration_min <= 0
        throw(ErrorException("Minimum collect duration must be positive."))
    end

    if collect_duration_max < collect_duration_min
        throw(ErrorException("Maximum collect duration must be greater or equal to minimum collect duration."))
    end

    # Array of requests to sample
    requests = Request[]

    # Sample locations
    lons = rand(Uniform(lon_min, lon_max), num_samples)
    lats = rand(Uniform(lat_min, lat_max), num_samples)

    if collect_duration_min == collect_duration_max
        cdur = [collect_duration_min for idx in 1:num_samples]
    else
        cdur = rand(Uniform(collect_duration_min, collect_duration_max), num_samples)
    end

    if reward_min == reward_max
        rewards = [reward_min for idx in 1:num_samples]
    else
        rewards = rand(Uniform(reward_min, reward_max), num_samples)
    end

    for i in 1:num_samples
        push!(requests, Request(lons[i], lats[i], look_angle_min=look_angle_min,
                look_angle_max=look_angle_max, collect_duration=cdur[i], 
                reward=rewards[i], use_degrees=true, id=i))
    end

    return requests
end


function sample_requests_file(num_samples::Int=0, filepath::String="";
            look_angle_min::Real=10, look_angle_max::Real=55,
            collect_duration_min::Real=30, collect_duration_max::Real=30,
            reward_min::Real=1.0, reward_max::Real=1.0)

    if reward_min <= 0
        throw(ErrorException("Minimum reward must be positive."))
    end

    if reward_max < reward_min
        throw(ErrorException("Maximum reward must be greater or equal to minimum reward."))
    end

    if collect_duration_min <= 0
        throw(ErrorException("Minimum collect duration must be positive."))
    end

    if collect_duration_max < collect_duration_min
        throw(ErrorException("Maximum collect duration must be greater or equal to minimum collect duration."))
    end

    json_requests = JSON.parsefile(filepath)
    num_images = length(json_requests)

    if num_samples > num_images
        throw(ErrorException("Number of requested samples greater than data points in file"))
    end
    
    @debug "Found $num_images images in file"
    indexes = collect(Int, 1:num_images)

    # Presample collect durations and rewards
    if collect_duration_min == collect_duration_max
        cdur = [collect_duration_min for idx in 1:num_samples]
    else
        cdur = rand(Uniform(collect_duration_min, collect_duration_max), num_samples)
    end

    if reward_min == reward_max
        rewards = [reward_min for idx in 1:num_samples]
    else
        rewards = rand(Uniform(reward_min, reward_max), num_samples)
    end

    # Array of requests to sample
    requests = Request[]

    observed_idx = []
    while length(requests) < num_samples
        idx = rand(indexes, 1)[1]
        if !(idx in observed_idx)
            # Create request
            i = length(requests) + 1
            lon = json_requests[idx]["lon"]
            lat = json_requests[idx]["lat"]
            req = Request(lon, lat, look_angle_min=look_angle_min,
                look_angle_max=look_angle_max, collect_duration=cdur[i], 
                reward=rewards[i], use_degrees=true, id=i)

            push!(requests, req)
            
            # Delete index to prevent resampling
            deleteat!(indexes, findfirst(x -> x == idx, indexes))
            # println(length(indexes))
        end
    end

    return requests
end

@with_kw mutable struct SolveSettings
    baseline::Bool=false
    graph::Bool=false
    milp::Bool=false
    mdpfs::Bool=false
    mdpmcts::Bool=false
    solve_gamma::Real = 1.0
    solve_depth::Int = 3
    solve_breadth::Int = 3
    solve_horizon::Real = 5400.0
    mcts_sim_iterations::Int = 10
    mcts_c::Real = 1.0
    milp_timeout::Real=900
end

function sat_sim(sc::Spacecraft,
        requests::Array{Request, 1},
        t_start::Epoch, t_end::Epoch; 
        stations::Array{GroundStation, 1}=GroundStation[],
        settings::SolveSettings)

    println("Started Executing Satellite Simulation")
    
    # Function returns
    baseline_plan, baseline_reward, baseline_time, baseline_nsr, baseline_pfs = Opportunity[], 0.0, 0.0, 0.0, 0.0
    graph_plan, graph_reward, graph_time, graph_nsr, graph_pfs = Opportunity[], 0.0, 0.0, 0.0, 0.0
    milp_plan, milp_reward, milp_time, milp_nsr, milp_pfs = Opportunity[], 0.0, 0.0, 0.0, 0.0
    mdpfs_plan, mdpfs_reward, mdpfs_time, mdpfs_nsr, mdpfs_pfs = Opportunity[], 0.0, 0.0, 0.0, 0.0
    mdpmcts_plan, mdpmcts_reward, mdpmcts_time, mdpmcts_nsr, mdpmcts_pfs = Opportunity[], 0.0, 0.0, 0.0, 0.0

    # Initialize Planning Problem
    problem = SatPlanningProblem(t_start=t_start, t_end=t_end)

    problem.solve_gamma = settings.solve_gamma
    problem.solve_depth = settings.solve_depth
    problem.solve_breadth = settings.solve_breadth
    problem.solve_horizon = settings.solve_horizon
    problem.mcts_sim_iterations = settings.mcts_sim_iterations
    problem.mcts_c = settings.mcts_c

    # Add constraints to problem
    push!(problem.constraints, constraint_agility_single_axis)

    # Add Spacecraft
    push!(problem.spacecraft, sc)

    # Add stations
    if length(stations) > 0
        add_locations(problem, stations)
    end

    # Add requests
    add_locations(problem, requests)

    # Compute Collects for problem
    pt = time()
    compute_access(problem, orbit_fraction=0.75)
    access_time = time() - pt
    @debug "Computed accesses in $access_time seconds."

    @debug "Found $(length(problem.contacts)) passes in planning horizon"
    @debug "Found $(length(problem.collects)) collects in planning horizon"
    num_accessible_requests = sum([1 for req in problem.requests if length(problem.lt_loc_opps[req.id]) > 0])
    pct_accessible_requests = floor(num_accessible_requests/length(problem.requests) * 10000)/100
    @debug "Found $num_accessible_requests out of $(length(problem.requests)) requests have accesses. $pct_accessible_requests%"

    # Baseline solve
    if settings.baseline == true
        @debug "Starting Baseline solve"
        pt = time()
        baseline_plan, baseline_reward = satellite_plan_baseline(problem)
        baseline_time = time() - pt

        baseline_feasible, baseline_realized_reward, baseline_nur, baseline_ndr, baseline_nr, baseline_nc = analyze_plan(problem, baseline_plan);
        baseline_reward = baseline_realized_reward
        baseline_nsr = baseline_nur
        baseline_pfs = baseline_nsr/num_accessible_requests

        if baseline_feasible == false
            @warn "Baseline solver produced infeasible solution for Request Count $(length(requests))"
        end
    end

    # Graph solve
    if settings.graph == true
        @debug "Starting Graph solve"
        pt = time()
        graph_plan, graph_reward = satellite_plan_graph(problem, sat_id=1)
        graph_time = time() - pt

        graph_feasible, graph_realized_reward, graph_nur, graph_ndr, graph_nr, graph_nc = analyze_plan(problem, graph_plan);
        graph_reward = graph_realized_reward
        graph_nsr = graph_nur
        graph_pfs = graph_nsr/num_accessible_requests

        if graph_feasible == false
            @warn "Graph solver produced infeasible solution for Request Count $(length(requests))"
        end
    end

    # MILP solve
    if settings.milp == true
        @debug "Starting MILP solve"
        pt = time()
        milp_plan, milp_reward = satellite_plan_milp(problem, sat_id=1, enforce_first=true, timeout=settings.milp_timeout)
        milp_time = time() - pt

        milp_feasible, milp_realized_reward, milp_nur, milp_ndr, milp_nr, milp_nc = analyze_plan(problem, milp_plan);
        milp_reward = milp_realized_reward
        milp_nsr = milp_nur
        milp_pfs = milp_nsr/num_accessible_requests

        if milp_feasible == false
            @warn "MILP solver produced infeasible solution for Request Count $(length(requests))"
        end
    end

    # MDP-FS solve
    if settings.mdpfs == true
        @debug "Starting MDP-FS solve"
        pt = time()
        mdpfs_states, mdpfs_plan, mdpfs_reward = satellite_plan_mdp_fs(problem)
        mdpfs_time = time() - pt

        mdpfs_feasible, mdpfs_realized_reward, mdpfs_nur, mdpfs_ndr, mdpfs_nr, mdpfs_nc = analyze_plan(problem, mdpfs_plan);
        mdpfs_reward = mdpfs_realized_reward
        mdpfs_nsr = mdpfs_nur
        mdpfs_pfs = mdpfs_nsr/num_accessible_requests

        if mdpfs_feasible == false
            @warn "MDP-FS solver produced infeasible solution for Request Count $(length(requests))"
        end
    end

    # MDP-MCTS solve
    if settings.mdpmcts == true
        @debug "Starting MDP-MCTS solve"
        pt = time()
        mdpmcts_states, mdpmcts_plan, mdpmcts_reward = satellite_plan_mdp_mcts(problem, parallel=false)
        mdpmcts_time = time() - pt

        mdpmcts_feasible, mdpmcts_realized_reward, mdpmcts_nur, mdpmcts_ndr, mdpmcts_nr, mdpmcts_nc = analyze_plan(problem, mdpmcts_plan);
        mdpmcts_reward = mdpmcts_realized_reward
        mdpmcts_nsr = mdpmcts_nur
        mdpmcts_pfs = mdpmcts_nsr/num_accessible_requests

        if mdpmcts_feasible == false
            @warn "MDP-MCTS produced infeasible solution for Request Count $(length(requests))"
        end
    end

    return problem, (baseline_plan, baseline_reward, baseline_time, baseline_nsr, baseline_pfs),
            (graph_plan, graph_reward, graph_time, graph_nsr, graph_pfs),
            (milp_plan, milp_reward, milp_time, milp_nsr, milp_pfs),
            (mdpfs_plan, mdpfs_reward, mdpfs_time, mdpfs_nsr, mdpfs_pfs),
            (mdpmcts_plan, mdpmcts_reward, mdpmcts_time, mdpmcts_nsr, mdpmcts_pfs)
end

function mass_sim(sc::Spacecraft,
    request_decks::Array{Array{Request, 1}, 1},
    t_start::Epoch, t_end::Epoch; 
    stations::Array{GroundStation, 1}=GroundStation[],
    settings::SolveSettings)

    # Create anonymous function to map work assignments into function all 
    fn(sc, ts, te, sta, set) = r -> sat_sim(sc, r, ts, te, stations=sta, settings=set)

    # Solve Each target deck in parallel
    results = pmap(fn(sc, t_start, t_end, stations, settings), request_decks)
 
    # Aggregate results
    baseline_rewards = [r[2][2] for r in results]
    baseline_times = [r[2][3] for r in results]
    baseline_nsr = [r[2][4] for r in results]
    baseline_pfs = [r[2][5] for r in results]

    graph_rewards = [r[3][2] for r in results]
    graph_times = [r[3][3] for r in results]
    graph_nsr = [r[3][4] for r in results]
    graph_pfs = [r[3][5] for r in results]

    milp_rewards = [r[4][2] for r in results]
    milp_times = [r[4][3] for r in results]
    milp_nsr = [r[4][4] for r in results]
    milp_pfs = [r[4][5] for r in results]

    mdpfs_rewards = [r[5][2] for r in results]
    mdpfs_times = [r[5][3] for r in results]
    mdpfs_nsr = [r[5][4] for r in results]
    mdpfs_pfs = [r[5][5] for r in results]

    mdpmcts_rewards = [r[6][2] for r in results]
    mdpmcts_times = [r[6][3] for r in results]
    mdpmcts_nsr = [r[6][4] for r in results]
    mdpmcts_pfs = [r[6][5] for r in results]


    return (baseline_rewards, baseline_times, baseline_nsr, baseline_pfs),
           (graph_rewards, graph_times, graph_nsr, graph_pfs),
           (milp_rewards, milp_times, milp_nsr, milp_pfs),
           (mdpfs_rewards, mdpfs_times, mdpfs_nsr, mdpfs_pfs),
           (mdpmcts_rewards, mdpmcts_times, mdpmcts_nsr, mdpmcts_pfs)
           
end

function sat_resource_sim()
end

function mass_resource_sim()
end