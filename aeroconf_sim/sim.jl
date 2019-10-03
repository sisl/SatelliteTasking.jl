
function run_landsat_sim(size::Integer, solve_discount::Real=0.99, 
            solve_depth::Integer=10, mcts_n_iterations::Integer=100,
            mcts_exploration_constant::Real=1.5; 
            enable_resources::Bool=false)

    ##
    ## Simulation Definition
    ##

    # Set simulation time window
    t_start = Epoch(2019, 1, 1, 0, 0, 0, tsys="UTC"); # Start of time span
    t_end = Epoch(2019, 1, 2, 0, 0, 0, tsys="UTC"); # End of simulation time span

    # Initialize Planning Problem
    problem = SatPlanningProblem(t_start=t_start, t_end=t_end);

    # Add constraints to problem
    push!(problem.constraints, constraint_agility_single_axis)

    # Create Spacecraft
    line1 = "1     1U          19  1.00000000  .00000000  00000-0  00000-0 0    05"
    line2 = "2     1  90.0000   0.0000 0010000   0.0000   0.0000 15.24162312    00"
    sc1 = Spacecraft(id=1, tle=TLE(line1, line2), slew_rate=1.0);
    push!(problem.spacecraft, sc1);

    # Compute Orital Period
    T = orbit_period(state(problem.spacecraft[1].tle, problem.spacecraft[1].tle.epoch)[1]);

    # Number of images and duration for full capacity
    DUTY_CYCLE = 10.0/100.0 # Spacecraft imaging duty cycle (as percentage)
    IMAGE_NUMBER   = 100
    IMAGE_DURATION = 30.0
    DOWNLINK_IMPROVEMENT = 10

    sc1.powergen_sunpoint =  1.0/((1.0 - DUTY_CYCLE)*T)             # Power/second
    sc1.powergen_image    = -1.0/(DUTY_CYCLE*T)                     # Power/second
    sc1.powergen_downlink = sc1.powergen_image/4.0                  # Power/second
    sc1.datagen_image     = 1.0/(IMAGE_NUMBER*IMAGE_DURATION)       # Data/second
    sc1.datagen_backorbit = sc1.datagen_image/10000.0               # Data/second
    sc1.datagen_downlink  = -DOWNLINK_IMPROVEMENT*sc1.datagen_image # Data/second

    println("Loadded $(length(problem.spacecraft)) spacecraft")

    # Set Solve Parameters
    problem.solve_discount = 0.99 # Typical values 0.99 - 0.9999
    problem.solve_depth = 3
    problem.solve_breadth = 3
    problem.solve_horizon = T

    println("Solve Parameters:\nsolve_discount: $(problem.solve_discount)\nsolve_depth: $(problem.solve_depth)\nsolve_breadth: $(problem.solve_breadth)\nsolve_horizon: $(problem.solve_horizon)")

    ##
    ## Load Deck
    ##

    
    # Load Groundstation Data
    stations = load_stations("./data/groundstations.json");
    num_stations = length(stations);

    # Load Image Data - [100, 200, 500, 1000, 2000, 5000, 10000, 12000, 15000]
    num_requests = 1000
    requests = load_requests("./data/landsat_test_$num_requests.json", id_offset=num_stations);
    @assert num_requests == length(requests);

    # Set Request Parameters
    collect_duration = 30.0
    set_collect_duration(requests, collect_duration)

    # Add Locations To Problem
    if enable_resources == true
        add_locations(problem, stations)
    end

    add_locations(problem, requests)

    ##
    ## Compute Opportunities
    ##

    # Compute Opportunities for collects and contacts

    # Remove all existing opportunities
    clear_opportunities(problem)

    # Compute problem accesses
    @time compute_access(problem, orbit_fraction=0.75)

    println("")
    println("Found $(length(problem.contacts)) passes in planning horizon")
    println("Found $(length(problem.collects)) collects in planning horizon")
    num_accessible_requests = sum([1 for req in problem.requests if length(problem.lt_loc_opps[req.id]) > 0])
    pct_accessible_requests = floor(num_accessible_requests/length(problem.requests) * 10000)/100
    println("Found $num_accessible_requests out of $(length(problem.requests)) requests have accesses. $pct_accessible_requests%")

    # MDP Planning Preparation - (Pre-compute Actions)
    @time precompute_action_space(problem, enable_resources=enable_resources);

    ##
    ## Set Solve Parameters
    ##

    problem.solve_discount = solve_discount # Typical values 0.99 - 0.9999
    problem.solve_depth = solve_depth
    problem.mcts_n_iterations = mcts_n_iterations
    problem.mcts_exploration_constant = mcts_exploration_constant

    ##
    ## Solve Planning Problem
    ##

    ts = time()
    mdp_mcts_states, mdp_mcts_plan, mdp_mcts_plan_rewards, mdp_mcts_reward = satellite_plan_mdp_mcts(problem, parallel=false)
    te = time()
    solve_time = te - ts

    # Check feasibility, and compute actual reward
    mdp_mcts_feasible, mdp_mcts_realized_reward, mdp_mcts_nur, mdp_mcts_ndr, mdp_mcts_nr, mdp_mcts_nc = analyze_plan(problem, mdp_mcts_plan);

    return (solve_discount, solve_depth, mcts_n_iterations, mcts_exploration_constant), mdp_mcts_realized_reward, solve_time
end