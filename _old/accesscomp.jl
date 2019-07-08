function compute_access(tle::TLE, locations::Array{<:Location, 1}, 
            epc_min::Union{Epoch, Nothing}, epc_max::Union{Epoch, Nothing}; 
            macro_step::Real=60.0, tol::Real=0.01, id_offset::Integer=0)
    if epc_min == nothing
        epc_min = tle.epoch
    end

    if epc_max == nothing
        epc_min += 86400.0
    end

    # Compute orbital period
    T = orbit_period(sCARTtoOSC(state(tle, epc_min), use_degrees=true)[1])

    opportunities = Array{Opportunity, 1}(undef, 0)

    next_step = macro_step

    for loc in locations
        # Perform macro adaptive stepsize search
        epc = epc_min
        while epc < epc_max

            # println("Current step: $epc")

            x_ecef = sECItoECEF(epc, state(tle, epc))

            if visible(x_ecef, loc)
                # println("Found instant of visibility - $(visible(x_ecef, loc)) - $epc")

                # Search for AOS (before initial guess epoch)
                window_open = find_access_boundary(tle, loc, epc, -macro_step, tol=tol)
                # println("Found window_open boundary: $window_open")

                # Search for LOS (after initial guess epoch)
                window_close = find_access_boundary(tle, loc, epc, macro_step, tol=tol)
                # println("Found window_close boundary: $window_close")

                teststep = (window_close - window_open)*1.1

                # If zero-doppler collection is required updated pass-times and geometry profile
                collect_duration = 0
                if typeof(loc) == Request && loc.require_zero_doppler == true
                    epc_mid = window_open + (window_close - window_open)/2.0
                    window_open  = epc_mid - loc.collect_duration/2.0
                    window_close = epc_mid + loc.collect_duration/2.0
                end

                # Create Oppourtunity
                if (window_close - window_open) > 0.0
                    if typeof(loc) == Request
                        collect = Collect(window_open, 
                            window_close,
                            spacecraft=spacecraft,
                            location=loc)

                        push!(opportunities, collect)
                    elseif typeof(loc) == GroundStation
                        contact = Contact(window_open, 
                            window_close,
                            spacecraft=spacecraft,
                            location=loc)

                        push!(opportunities, contact)
                    else
                        throw(ErrorException("Unknown type of location: $loc"))
                    end
                end

                # Step a most of an orbit because there (likely) won't be another
                # Acces until at least 1 orbit later
                # epc += T*(13/16)
                epc += T*(1/2)
                # epc += T/10
                # epc += teststep

            else
                # Take a macro step because there's a good chance there won't be an
                # opportunity until the next one
                epc += macro_step
            end
        end
    end

    # Sort opportunities in ascending order
    sort!(opportunities, by = x -> x.t_start)

    # Apply IDs for collects in ascending time order
    for (idx, opp) in enumerate(opportunities)
        # Override computed opportunity id
        opp.id = idx + id_offset
    end

    return opportunities
end

function parallel_compute_access(tle::TLE, locations::Array{<:Location, 1},
    epc_min::Union{Epoch, Nothing}, epc_max::Union{Epoch, Nothing}; 
    macro_step::Real=60.0, tol::Real=0.01, id_offset::Integer=0)

# Create anonymous function to map array inputs to 
fn(a, b, c, d, e) = x -> compute_access(a, x, b, c, macro_step=d, tol=e)

# Create work assignments
nw = nworkers()        # Number of workers
lw = length(locations) # Length of work

# @debug "Have $lw items of work and $nw workers")

wa = floor(Int, lw/nw) # Average work per worker
wr = lw - nw*wa        # Remaining worker

# Work Array
assignments = Array{<:Location, 1}[]

# @debug "Worker 1 convering $(1):$(wa+wr)")
push!(assignments, locations[1:(wa+wr)])
for i in 2:nw
# @debug "Worker $i convering $(1+wr+(i-1)*wa):$(i*wa+wr)")
push!(assignments, locations[(1+wr+(i-1)*wa):(i*wa+wr)])
end

results = pmap(fn(tle, epc_min, epc_max, macro_step, tol), assignments)

# Aggregate all results
opportunities = vcat(results...)

# Sort opportunities in ascending order
sort!(opportunities, by = x -> x.t_start)

# Aggregate results and resolve id numbering
for (idx, opp) in enumerate(opportunities)
# Override computed opportunity id
opp.id = idx + id_offset
end

return opportunities
end