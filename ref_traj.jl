function create_ref_hover(N, dt, n_inputs, mass)
    # create a hover reference (with no pitch) for whole time period

    # hover altitude
    desired_altitude= 5.0

    # hover state
    hover_position = [0.0; 0.0; desired_altitude]  # Hover at (0, 0, 5)
    hover_velocity = [0.0; 0.0; 0.0]  # no velocity
    hover_attitude = [0.0; 0.0; 0.0]  # No roll, pitch, yaw
    hover_angvelocity= [0.0; 0.0; 0.0] # no angular velocity
    
    Xref= [zeros(3) for i = 1:N]
    for i= 1:N
        Xref[i]= [hover_position; hover_velocity; hover_attitude; hover_angvelocity]
    end
    
    Uref = [(9.81*mass/n_inputs)*ones(n_inputs) for i = 1:(N-1)]
    return Xref, Uref
end