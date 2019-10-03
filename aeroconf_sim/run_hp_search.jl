import Pkg

# Julia Parallel Packages
using Distributed

# Julia Packages
using Random
using Profile
using Printf
using Logging
using Statistics
using Distributions
using StatProfilerHTML

## Set Worker Process Number 
NUM_WORKERS = 4
addprocs(NUM_WORKERS)
println("Started $(nprocs()) total processes - $(nworkers()) worker processes.")

# Allow for easy updating
@everywhere using Revise

# Satellite Dynamics Packages
@everywhere using SatelliteDynamics

# Load SatelliteTasking
@everywhere using SatelliteTasking

# Add Map Plots
using PGFMapPlots

include("./sim.jl")

######################################
# Define Hyperparameter Search Space #
######################################

# Define parameter space
ec_space = [0.1, 1.0, 1.5, 3.0, 5.0, 10.0]
ni_space = [25, 50, 100, 200, 500]
sd_space = [5, 10, 20, 30, 50] # solve depth
ds_space = [0.99, 0.995, 0.999] # discount
bin_iterations = 10

ec_space = [1.5, 2.0, 3.0]
ni_space = [10, 20]
sd_space = [10, 20] # solve depth
ds_space = [0.99, 0.995] # discount
bin_iterations = 3

# Create reward storage
hyper_time = Dict{Tuple{Real, Real, Real, Real}, AbstractVector{Real}}()
hyper_reward = Dict{Tuple{Real, Real, Real, Real}, AbstractVector{Real}}()

println("Total number of required simulations: $(length(ec_space)*length(ni_space)*length(sd_space)*length(ds_space)*bin_iterations)")
num_jobs = length(ec_space)*length(ni_space)*length(sd_space)*length(ds_space)
println("Total number of jobs: $num_jobs")

#######################
# Run Parallel Search #
#######################

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
                push!(work, (ds, sd, ni, ec))
            end
        end
    end
end


# Define Execution Function
fn(deck_size, er) = w -> run_landsat_sim(deck_size, w[1], w[2], w[3], w[4]; enable_resources=er)

# Execute Work in parallel

println("Starting hyper parameter search")
hp_start = time()
results = pmap(fn(1000, false), work)
hp_end = time()
println("Parameter search finished after $(hp_end-hp_start) seconds.")
println(results)

# Merge results
hyper_time = Dict{Tuple{Real, Real, Real, Real}, AbstractVector{Real}}()
hyper_reward = Dict{Tuple{Real, Real, Real, Real}, AbstractVector{Real}}()

for (solve_params, solve_reward, solve_time) in results
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