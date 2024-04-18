" TRAJECTORIES "

###### UTILITY FUNCTIONS ######

# determines the amount of thrust required to hover
function calculate_hover_thrust(mass::Float64)
    # Acceleration due to gravity (m/s^2)
    g = 9.81
    # Calculate weight (force exerted by gravity)
    weight = mass * g
    # The thrust required for hover is equal to the weight of the drone
    hover_thrust = 1* weight
    
    return hover_thrust
end



" 1. Hover trajectory

- will create a reference trajectory (Xref, Uref) that is just the hexrotor hovering in one location
- will start on the ground, rise to a certain altitude and maintain hover
- At the end of the time period, it will lower back to the ground.


"

# vertical ascent
function generate_vertical_ascent_trajectory(model::NamedTuple, desired_altitude::Float64)
    # required thrust to hover
    mass= model.mass
    dt= model.dt
    hover_thrust= calculate_hover_thrust(mass)

    # Desired position is constant at ground level
    rx = 0.0
    ry = 0.0
    rz = min(desired_altitude * dt / ascent_duration, desired_altitude)  # Linear ascent
    # do the minimum between getting to the desired altitude and being at it
    
    # Desired velocities are zero
    vx = 0.0
    vy = 0.0
    vz = min(desired_altitude / ascent_duration, desired_altitude)  # Constant ascent rate
    
    # Desired attitude (roll, pitch, yaw) and angular velocities are zero
    p1= 0.0
    p2= 0.0
    p3= 0.0

    ω1= 0.0
    ω2= 0.0
    ω3= 0.0

    
    # Construct reference trajectory vector
    X_ref = [rx, ry, rz, vx, vy, vz, p1, p2, p3, ω1, ω2, ω3]
    
    # Use hover thrust for ascent
    U_ref = [hover_thrust, hover_thrust, hover_thrust, hover_thrust, hover_thrust, hover_thrust]
    
    return X_ref, U_ref
end

# maintain hover
function generate_hover_trajectory(t::Float64, desired_altitude::Float64, hover_duration::Float64, hover_thrust::Float64)
    # Desired position is constant at desired altitude
    x_ref = 0.0
    y_ref = 0.0
    z_ref = desired_altitude
    
    # Desired velocities and accelerations are zero
    dx_dt_ref = 0.0
    dy_dt_ref = 0.0
    dz_dt_ref = 0.0
    
    # Desired attitude (roll, pitch, yaw) and angular velocities are zero
    phi_ref = 0.0
    theta_ref = 0.0
    psi_ref = 0.0
    p_ref = 0.0
    q_ref = 0.0
    r_ref = 0.0
    
    # Construct reference trajectory vector
    if t < hover_duration
        X_ref = [x_ref, y_ref, z_ref, dx_dt_ref, dy_dt_ref, dz_dt_ref, phi_ref, theta_ref, psi_ref, p_ref, q_ref, r_ref]
    else
        X_ref = [x_ref, y_ref, z_ref, dx_dt_ref, dy_dt_ref, dz_dt_ref, phi_ref, theta_ref, psi_ref, p_ref, q_ref, r_ref]
    end
    
    # Use hover thrust for hover
    U_ref = [hover_thrust, hover_thrust, hover_thrust, hover_thrust, hover_thrust, hover_thrust]
    
    return X_ref, U_ref
end




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