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