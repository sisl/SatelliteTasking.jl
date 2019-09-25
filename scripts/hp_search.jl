###########
# Imports #
###########

# Julia Parallel Packages
using Distributed

# Julia Packages
using Dates
using Random
using Profile
using Printf
using Logging
using Statistics
using Distributions

# Register Number of processes
# Julia must be started with `julia -p auto -e "using IJulia; jupyterlab();"`
# On Computer with 6 cores add 6-1 = 5 processes
# NOTE: This must be called _before_ @everywhere is used those functions get propery distributed    
if length(procs()) == 1
    # Remove existing process IDs so reloading works
    addprocs(length(Sys.cpu_info()));
    println("Started $(nprocs()) total processes - $(nworkers()) worker processes.")
else
    println("Found $(nprocs()) total processes - $(nworkers()) worker processes.")
end

# Satellite Dynamics Packages
@everywhere using SatelliteDynamics

# Load SatelliteTasking
@everywhere using SatelliteTasking

# Add Map Plots
using PGFMapPlots

#########################
# Simulation Definition #
#########################

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

println("Loadded $(length(problem.spacecraft)) spacecraft")

# Set Solve Parameters
problem.solve_discount = 0.99 # Typical values 0.99 - 0.9999
problem.solve_depth = 3
problem.solve_breadth = 3
problem.solve_horizon = T

println("Default Solve Parameters:\nsolve_discount: $(problem.solve_discount)\nsolve_depth: $(problem.solve_depth)\nsolve_breadth: $(problem.solve_breadth)\nsolve_horizon: $(problem.solve_horizon)")

######################
# Load Location Data #
######################

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
# add_locations(problem, stations)
add_locations(problem, requests)

#########################
# Compute Opportunities #
#########################

# Compute problem accesses
compute_access(problem, orbit_fraction=0.75);

println("")
println("Found $(length(problem.contacts)) passes in planning horizon")
println("Found $(length(problem.collects)) collects in planning horizon")
num_accessible_requests = sum([1 for req in problem.requests if length(problem.lt_loc_opps[req.id]) > 0])
pct_accessible_requests = floor(num_accessible_requests/length(problem.requests) * 10000)/100
println("Found $num_accessible_requests out of $(length(problem.requests)) requests have accesses. $pct_accessible_requests%")

# MDP Planning Preparation - (Pre-compute Actions)
println("Precomputing action space")
precompute_action_space(problem, enable_resources=false);

########################
# Define Work Function #
########################

@everywhere function run_mcts_sim(problem::SatPlanningProblem, n::Int, params::Tuple{Real, Real, Real, Real})
    println("Reporting in!")
    # Create reward storage
    hyper_time = Dict{Tuple{Real, Real, Real, Real}, AbstractVector{Real}}()
    hyper_reward = Dict{Tuple{Real, Real, Real, Real}, AbstractVector{Real}}()

    # Extract search point from params
    ec, ni, sd, ds = params

    # Copy problem
    # problem = copy(problem)
    problem.solve_discount = ds # Typical values 0.99 - 0.9999
    problem.solve_depth = sd
    problem.mcts_n_iterations = ni
    problem.mcts_exploration_constant = ec

    hyper_time[(ec, ni, sd, ds)]   = Real[]
    hyper_reward[(ec, ni, sd, ds)] = Real[]

    # Run specified number of iterations
    for idx in 1:n
        start_time = time()
        
        # Solve problem 
        println("Iteration $idx - Solving Problem with parameters:($(problem.solve_discount), $(problem.solve_depth), $(problem.mcts_n_iterations), $(problem.mcts_exploration_constant))")
        mdp_mcts_states, mdp_mcts_plan, mdp_mcts_reward = satellite_plan_mdp_mcts(problem, parallel=false)
        end_time = time()
        
        # Extract reward
        mdp_mcts_feasible, mdp_mcts_realized_reward, mdp_mcts_nur, mdp_mcts_ndr, mdp_mcts_nr, mdp_mcts_nc = analyze_plan(problem, mdp_mcts_plan);
        
        # Update serach values
        push!(hyper_time[(ec, ni, sd, ds)], end_time-start_time)
        push!(hyper_reward[(ec, ni, sd, ds)], mdp_mcts_realized_reward)
    end

    return hyper_time, hyper_reward
end

######################################
# Define Hyperparameter Search Space #
######################################

# Define parameter space
ec_space = [0.01, 0.1, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]
ni_space = [10, 20, 30, 50, 100, 200, 500, 1000]
sd_space = [2, 3, 5, 10, 20, 30, 50] # solve depth
ds_space = [0.99, 0.995, 0.999] # discount
bin_iterations = 10

ec_space = [1.5, 2.0, 3.0]
ni_space = [10, 20]
sd_space = [10, 20] # solve depth
ds_space = [0.99, 0.995] # discount
bin_iterations = 3

# ec_space = [0.01, 0.1, 0.2]
# ni_space = [5, 10]
# sd_space = [2, 3, 5]
# ds_space = [0.99]
# bin_iterations = 10

# Create reward storage
hyper_time = Dict{Tuple{Real, Real, Real, Real}, AbstractVector{Real}}()
hyper_reward = Dict{Tuple{Real, Real, Real, Real}, AbstractVector{Real}}()

println("Total number of required simulations: $(length(ec_space)*length(ni_space)*length(sd_space)*length(ds_space)*bin_iterations)")
num_jobs = length(ec_space)*length(ni_space)*length(sd_space)*length(ds_space)
println("Total number of jobs: $num_jobs")

#######################
# Run Parallel Search #
#######################

fn(p, n) = x -> run_mcts_sim(p, n, x)

# Create work assignments
nw = nworkers()        # Number of workers
lw = num_jobs # Length of work

wa = floor(Int, lw/nw) # Average work per worker
wr = lw - nw*wa        # Remaining worker

println("Have $lw items of work and $nw workers")
println("Average work: $wa, Remaining work: $wr")

# Create search points 
work = []
for ec in ec_space
    for ni in ni_space
        for sd in sd_space
            for ds in ds_space
                push!(work, (ec, ni, sd, ds))
            end
        end
    end
end

println("Work: $work")
# Construct assignments (function inputs) for Workers
# assignments = []
# println("Assignments: $assignments")

# push!(assignments, work[1:(wa+wr)])
# for i in 2:min(nw, length(work))
#     # @debug "Worker $i convering $(1+wr+(i-1)*wa):$(i*wa+wr)")
#     push!(assignments, work[(1+wr+(i-1)*wa):(i*wa+wr)])
# end

# Execute Work in parallel

println("Starting hyper parameter search")
hp_start = time()
results = pmap(fn(problem, bin_iterations), work)
# results = pmap(fn(problem, bin_iterations), assignments)
hp_end = time()
println("Parameter search finished after $(hp_end-hp_start) seconds.")

# Merge results
hyper_time = Dict{Tuple{Real, Real, Real, Real}, AbstractVector{Real}}()
hyper_reward = Dict{Tuple{Real, Real, Real, Real}, AbstractVector{Real}}()

for r in results
    ht, hr = r
    for params, time_values in ht
        hyper_time[params] = time_values
        hyper_reward[params] = hr[params]
    end
end

#################################
# Extract Best Hyper Parameters #
#################################

# Extract Maximum reward
max_max = -Inf
max_mean = -Inf
max_min = - Inf
max_hp = (-Inf, -Inf, -Inf, -Inf)

for (hparams, rvec) in hyper_reward
    minval, maxval = extrema(rvec)
    if length(rvec) > 0 && maxval > max_max
        max_max = maxval
        max_mean = mean(rvec)
        max_min = minval
        max_hp = hparams
    end
end

# Print Best Point
println("Best Hyperparameter Cominbation -- Min: $(max_min) - Mean: $(max_mean) - Max: $(max_max)")
println("Exploration Constant: $(max_hp[1])")
println("Number of Iterations: $(max_hp[2])")
println("Solve Depth: $(max_hp[3])")
println("Discount Facgtor: $(max_hp[4])")

# Save Rest to File
open("hyperparam_search_rewards_$(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")).dat") do file
    column_string = join(["Reward $idx" for idx in 1:bin_iterations], ",")
    write(file, "Exploration Constant, Number Of Iterations, Search Depth, Discount Factor, $column_string\n")
    for hparams, rewards in hyper_reward
        ec, ni, sd, ds = hparams
        reward_string = join(rewards, ",")
        write(file, "$(ec),$ni,$sd,$ds,$reward_string\n")
    end
end

open("hyperparam_search_times_$(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")).dat") do file
    column_string = join(["Runtime $idx" for idx in 1:bin_iterations], ",")
    write(file, "Exploration Constant, Number Of Iterations, Search Depth, Discount Factor, $column_string\n")
    for hparams, times in hyper_time
        ec, ni, sd, ds = hparams
        time_string = join(times, ",")
        write(file, "$(ec),$ni,$sd,$ds,$time_string\n")
    end
end