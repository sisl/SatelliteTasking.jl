import Pkg

# Julia Packages
using PGFPlots
using Profile
using StatProfilerHTML

# Allow for easy updating
using Revise

# Satellite Dynamics Packages
using SatelliteDynamics

# Load SatelliteTasking
using SatelliteTasking

# Add Map Plots
using PGFMapPlots

##
# Simulation Definition
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

println("Loadded $(length(problem.spacecraft)) spacecraft")

##
# Load Groundstation Data
##

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

##
# Compute Opportunities for collects and contacts
##

# Remove all existing opportunities
clear_opportunities(problem)

# Compute problem accesses
@time compute_access(problem, orbit_fraction=0.75)

println("Found $(length(problem.contacts)) passes in planning horizon")
println("Found $(length(problem.collects)) collects in planning horizon")
num_accessible_requests = sum([1 for req in problem.requests if length(problem.lt_loc_opps[req.id]) > 0])
pct_accessible_requests = floor(num_accessible_requests/length(problem.requests) * 10000)/100
println("Found $num_accessible_requests out of $(length(problem.requests)) requests have accesses. $pct_accessible_requests%")

##
# Baseline Solve (Single-Threaded)
##

# @time baseline_plan, baseline_reward = satellite_plan_baseline(problem, allow_repeats=false)

# # Check feasibility, and compute actual reward
# baseline_feasible, baseline_realized_reward, baseline_nur, baseline_ndr, baseline_nr, baseline_nc = analyze_plan(problem, baseline_plan);

# println("Computed plan is feasible: $(uppercase(string(baseline_feasible)))")
# println("Computed optimal schedule. Computed reward: $baseline_reward - Realized reward: $baseline_realized_reward.")
# println("Requests Satisfied: $baseline_nur/$(length(problem.requests)), Duplicate Collects: $baseline_ndr")
# println("Ground Contacts Taken: $baseline_nc")

##
# Graph Solve (Single-Threaded)
##

# @time graph_plan, graph_reward = satellite_plan_graph(problem, allow_repeats=false)

# # Check feasibility, and compute actual reward
# graph_feasible, graph_realized_reward, graph_nur, graph_ndr, graph_nr, graph_nc = analyze_plan(problem, graph_plan);

# println("Computed plan is feasible: $(uppercase(string(graph_feasible)))")
# println("Computed optimal schedule. Computed reward: $graph_reward - Realized reward: $graph_realized_reward.")
# println("Requests Satisfied: $graph_nur/$(length(problem.requests)), Duplicate Collects: $graph_ndr")
# println("Ground Contacts Taken: $graph_nc")

##
# Mixed-Integer Linear Programing Solve (Single-Threaded)
##

# @time milp_plan, milp_reward = satellite_plan_milp(problem, allow_repeats=false)

# # Check feasibility, and compute actual reward
# milp_feasible, milp_realized_reward, milp_nur, milp_ndr, milp_nr, milp_nc = analyze_plan(problem, milp_plan);

# println("")
# println("Computed plan is feasible: $(uppercase(string(milp_feasible)))")
# println("Computed optimal schedule. Computed reward: $milp_reward - Realized reward: $milp_realized_reward.")
# println("Requests Satisfied: $milp_nur/$(length(problem.requests)), Duplicate Collects: $milp_ndr")
# println("Ground Contacts Taken: $milp_nc")

##
# Markov Decision Process (No Resources)
##

# # Adjust solve parameters
# problem.solve_discount = 1.0 # Typical values 0.999 - 0.9999
# problem.solve_depth = 5
# problem.solve_breadth = 3
# problem.solve_horizon = T

# @time mdp_fs_plan, mdp_fs_reward = satellite_plan_mdp_fs(problem)

# # Check feasibility, and compute actual reward
# mdp_fs_feasible, mdp_fs_realized_reward, mdp_fs_nur, mdp_fs_ndr, mdp_fs_nr, mdp_fs_nc = analyze_plan(problem, mdp_fs_plan);

# println("Computed plan is feasible: $(uppercase(string(mdp_fs_feasible)))")
# println("Computed optimal schedule. Computed reward: $mdp_fs_reward - Realized reward: $mdp_fs_realized_reward.")
# println("Requests Satisfied: $mdp_fs_nur/$(length(problem.requests)), Duplicate Collects: $mdp_fs_ndr")
# println("Ground Contacts Taken: $mdp_fs_nc")

# # Profile MDP Planning code
# @profilehtml satellite_plan_mdp_fs(problem);

##
# Markov Decision Process - Monte Carlo Tree Search (No Resources)
##

# Adjust solve parameters
problem.solve_allow_repeats = false
problem.mdp_reward_scarcity = false
problem.solve_discount = 1.0 # Typical values 0.99 - 0.9999
problem.solve_depth = 10
problem.solve_breadth = 0
problem.solve_horizon = T
problem.mcts_sim_iterations = 30
problem.mcts_c = 3.0

mdp_mcts_plan, mdp_mcts_reward = satellite_plan_mdp_mcts(problem, parallel=false)

# Check feasibility, and compute actual reward
mdp_mcts_feasible, mdp_mcts_realized_reward, mdp_mcts_nur, mdp_mcts_ndr, mdp_mcts_nr, mdp_mcts_nc = analyze_plan(problem, mdp_mcts_plan);

println("Computed plan is feasible: $(uppercase(string(mdp_mcts_feasible)))")
println("Computed optimal schedule. Computed reward: $mdp_mcts_reward - Realized reward: $mdp_mcts_realized_reward.")
println("Requests Satisfied: $mdp_mcts_nur/$(length(problem.requests)), Duplicate Collects: $mdp_mcts_ndr")
println("Ground Contacts Taken: $mdp_mcts_nc")

# Profile MCTS 
@profilehtml satellite_plan_mdp_mcts(problem);