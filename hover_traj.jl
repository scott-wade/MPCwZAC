" TRAJECTORIES "
" 1. Hover trajectory

- will create a reference trajectory (Xref, Uref) that is just the hexrotor hovering in one location
- will start on the ground, rise to a certain altitude and maintain hover
- At the end of the time period, it will lower back to the ground.


"

# vertical ascent
function generate_vertical_ascent_trajectory(model::NamedTuple, N, desired_altitude::Float64, ascent_duration)
    # parameters from model
    mass= model.mass
    dt= model.dt
    
    # Desired position is constant at ground level
    rx = 0.0
    ry = 0.0
    rz = min(desired_altitude * dt / ascent_duration, desired_altitude)  # Linear ascent
    # do the minimum between getting to the desired altitude and being at it
    ascent_position= [rx; ry; rz]
    
    # Desired velocities are zero
    vx = 0.0
    vy = 0.0
    vz = min(desired_altitude / ascent_duration, desired_altitude)  # Constant ascent rate
    ascent_velocity= [vx; vy; vz]
    
    # Desired attitude (roll, pitch, yaw) and angular velocities are zero
    ascent_attitude= [0.0; 0.0; 0.0]
    ascent_angvelocity= [0.0; 0.0; 0.0]

    # Construct reference trajectory vector
    # Xref= [zeros(12) for i = 1:N] 
    # for i= 1:N
    #     Xref[i]= [ascent_position; ascent_velocity; ascent_attitude; ascent_angvelocity]
    # end
    Xref= [ascent_position; ascent_velocity; ascent_attitude; ascent_angvelocity]
    
    # Linearly interpolate control thrust based on desired altitude
    hover_thrust= (9.81*mass)
    initial_thrust = 0.0
    final_thrust = hover_thrust
    current_thrust = initial_thrust + (final_thrust - initial_thrust) * min(dt / ascent_duration, 1.0)

    # Use hover thrust for ascent
    #Uref = [(current_thrust/6)*ones(6) for i = 1:(N-1)]
    Uref = [(current_thrust/6); (current_thrust/6); (current_thrust/6); (current_thrust/6); (current_thrust/6); (current_thrust/6)]
    
    return Xref, Uref
end

# maintain hover
function create_ref_hover(model, N, desired_altitude)
    mass= model.mass
    #g= model.g

    # hover state
    hover_position = [0.0; 0.0; desired_altitude]  # Hover at (0, 0, altitude)
    hover_velocity = [0.0; 0.0; 0.0]  # no velocity
    hover_attitude = [0.0; 0.0; 0.0]  # No roll, pitch, yaw
    hover_angvelocity= [0.0; 0.0; 0.0] # no angular velocity
    
    # Xref= [zeros(12) for i = 1:N]
    # for i= 1:N
    #     Xref[i]= [hover_position; hover_velocity; hover_attitude; hover_angvelocity]
    # end
    Xref= [hover_position; hover_velocity; hover_attitude; hover_angvelocity]

    hover_thrust= 9.81* mass
    #Uref = [(9.81*mass/6)*ones(6) for i = 1:(N-1)]
    Uref= [hover_thrust/6; hover_thrust/6; hover_thrust/6; hover_thrust/6; hover_thrust/6; hover_thrust/6]
    
    return Xref, Uref
end

# descend back to ground at a linear rate
function generate_vertical_descent_trajectory(model, N, desired_altitude::Float64, descent_duration)
    # parameters from model
    mass= model.mass
    dt= model.dt
    
    # Desired position during descent
    rx = 0.0
    ry = 0.0
    rz = desired_altitude - min(desired_altitude * dt / descent_duration, desired_altitude)  # Linear descent
    descent_position= [rx; ry; rz]

    # Desired velocities are zero
    vx = 0.0
    vy = 0.0
    vz = min(desired_altitude / descent_duration, desired_altitude)  # Constant descent rate
    descent_velocity= [vx; vy; vz]

    # Desired attitude (roll, pitch, yaw) and angular velocities are zero
    descent_attitude= [0.0; 0.0; 0.0]
    descent_angvelocity= [0.0; 0.0; 0.0]
    
    # Construct reference trajectory vector
    # Xref= [zeros(12) for i = 1:N] 
    # for i= 1:N
    #     Xref[i]= [descent_position; descent_velocity; descent_attitude; descent_angvelocity]
    # end
    Xref= [descent_position; descent_velocity; descent_attitude; descent_angvelocity]
    
    # Linearly interpolate control thrust based on desired altitude for descent
    hover_thrust= (9.81*mass)
    initial_thrust = hover_thrust  # Start at hover thrust
    final_thrust = 0.0  # End at zero thrust
    current_thrust = initial_thrust - (initial_thrust - final_thrust) * min(dt / descent_duration, 1.0)

    # Use hover thrust for ascent
    #Uref = [(current_thrust/6)*ones(6) for i = 1:(N-1)]
    Uref = [(current_thrust/6); (current_thrust/6); (current_thrust/6); (current_thrust/6); (current_thrust/6); (current_thrust/6)]
    
    return Xref, Uref
end


# generate total Xref, and Uref
function generate_trajectory(model, N, t_vec, total_duration::Float64, ascent_duration::Float64, hover_duration::Float64, descent_duration::Float64, desired_altitude::Float64)
    # Initialize arrays to store reference states and controls
    X_ref = []
    U_ref = []
    
    # Generate vertical ascent trajectory
    for t in t_vec
        # vertical ascent
        if t_vec[t] <= ascent_duration
            Xref, Uref = generate_vertical_ascent_trajectory(model, N, desired_altitude, ascent_duration)
            push!(X_ref, Xref)
            push!(U_ref, Uref)
        # hover
        elseif t_vec[t] > ascent_duration && t <= (ascent_duration + hover_duration)
            Xref, Uref = generate_hover_trajectory(model, N, desired_altitude)
            push!(X_ref, Xref)
            push!(U_ref, Uref)
        # vertical descent
        else t_vec[t] > (ascent_duration + hover_duration) && t <= (total_duration)
            Xref, Uref = generate_vertical_descent_trajectory(t - (ascent_duration + hover_duration), desired_altitude, hover_thrust)
            push!(X_ref, Xref)
            push!(U_ref, Uref)
        end
        
    end

    return X_ref, U_ref
end


"""
# parameters to use in the code

dt = 0.01 (from model)
total_duration = 20.0 (20 second test?)
t_vec = 0:dt:total_duration 
N = length(t_vec)
desired_altitude= 10.0
ascent_duration= 5.0 (sec)
hover_duration= 10.0 (sec)
descent_duration= 5.0 (sec)


# Generate hover trajectory for hexrotor
Xref, Uref = create_ref_hover(model, N, desired_altitude, ascent_duration)
Um = hcat(Uref...)
Xm = hcat(Xref...)

# plotting here...

"""






" 2. Planar Figure 8








"

# generate Xref vector 
function generate_reference_trajectory(t::Float64)
    # Define parameters for figure 8 trajectory
    amplitude = 1.0
    frequency = 1.0
    omega = 2 * π * frequency
    
    # Calculate desired position (x, y) for figure 8 trajectory
    x_ref = amplitude * sin(omega * t)
    y_ref = amplitude * sin(2 * omega * t)
    
    # Maintain constant altitude achieved during vertical ascent
    z_ref = 10.0
    
    # Construct reference trajectory vector
    X_ref = [x_ref, y_ref, z_ref]
    
    return X_ref
end

function generate_reference_controls(X_ref::Vector{Float64})
    # For vertical ascent and descent, use thrust required for hover
    # For figure 8 trajectory, adjust thrust to maintain constant altitude
    
    # Calculate thrust required for hover
    hover_thrust = calculate_hover_thrust(mass)
    
    # For ascent, set thrust to hover thrust
    if X_ref[3] < desired_altitude
        U_ref = [0.0, 0.0, hover_thrust]
    # For figure 8 trajectory, adjust thrust to maintain constant altitude
    elseif t < t_figure_8_end
        # Adjust thrust to maintain constant altitude
        U_ref = [thrust_x, thrust_y, adjust_thrust_for_altitude(X_ref[3])]
    # For descent, gradually reduce thrust to hover thrust
    else
        U_ref = [0.0, 0.0, hover_thrust]
    end
    
    return U_ref
end