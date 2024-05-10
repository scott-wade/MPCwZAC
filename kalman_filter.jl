
# Kalman Filter 
function kalman_filter(x̂, y, wind_accel, u_k, Σ_k, A, B, wind_dynamics, W, V) # Xsim (which is augmented with the wind model), Usim, covariance Σ from previous kalman run
    # augment the state in order to include wind dynamics
    x̃_k= [x̂; wind_accel]
      
    # state matrix with wind
    M= [A wind_dynamics; zeros(3,12) I(3)]
    C_aug= [I(12) zeros(12,3)]
    B_aug= [B; zeros(3,6)]

    # prediction
    x̃_next_pred= M * x̃_k + B_aug* u_k  #wind modeled as the a random walk in sim, not here
    Σ_next_pred= M * Σ_k * M' + W 

    # measurement, innovation, and innovation covariance
    # assuming the whole state is observable
    z= y .- C_aug * x̃_next_pred # innovation
    S= C_aug * Σ_next_pred * C_aug' + V

    # kalman gain
    L= Σ_next_pred * C_aug' * inv(S)
    # println("vkstate= ", size(v_k_state))
    # println("xk= ", size(C*x̃_k))
    # println("y= ", size(y))
    # println("z= ", size(z))
    # println("L= ",size(L))

    # state update/correction
    x̃_next= x̃_next_pred + L * z

    x̃_next_states= x̃_next[1:12]
    estimated_wind= x̃_next[13:15]

    # covariance update
    Σ_next= (I- L*C_aug)* Σ_next_pred* (I-L*C_aug)' + L*V*L'

    return x̃_next_states, estimated_wind, Σ_next
end