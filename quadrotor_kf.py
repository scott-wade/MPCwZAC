import numpy as np
import matplotlib.pyplot as plt

def quadrotor_dynamics(state, inputs, dt):
    """
    Example quadrotor dynamics without wind model uncertainty.

    Inputs:
    - state: np.array, current state [x, y, z, vx, vy, vz, phi, theta, psi, p, q, r, wx, wy, wz]
    - inputs: np.array, inputs [thrust, tau_phi, tau_theta, tau_psi]
    - dt: float, timestep

    Returns:
    - new_state: np.array, next state after dt
    """
    # Unpack the state
    x, y, z, vx, vy, vz, phi, theta, psi, p, q, r, wx, wy, wz = state
    
    # Unpack the inputs
    thrust, tau_phi, tau_theta, tau_psi = inputs
    
    # Constants
    m = 1.0  # Mass of the quadrotor
    g = 9.81  # Gravitational acceleration
    Ix = 0.01  # Moment of inertia around x-axis
    Iy = 0.01  # Moment of inertia around y-axis
    Iz = 0.02  # Moment of inertia around z-axis
    
    # Translational dynamics
    ax = 0 + wx
    ay = 0 + wy
    az = thrust / m - g + wz
    
    # Rotational dynamics
    p_dot = tau_phi / Ix
    q_dot = tau_theta / Iy
    r_dot = tau_psi / Iz

    # add random walk wind noise
    wx += np.random.normal(0, 0.2, 1)
    wy += np.random.normal(0, 0.2, 1)
    wz += np.random.normal(0, 0.2, 1)
    
    # Integrate to get the new state
    new_state = np.array([
        x + vx * dt,
        y + vy * dt,
        z + vz * dt,
        vx + ax * dt,
        vy + ay * dt,
        vz + az * dt,
        phi + p * dt,
        theta + q * dt,
        psi + r * dt,
        p + p_dot * dt,
        q + q_dot * dt,
        r + r_dot * dt,
        wx,
        wy,
        wz
    ])


    
    return np.reshape(new_state, n_states)

import numpy as np

def linearize(state_eq, input_eq, dt):
    """
    Linearize the quadrotor dynamics around a given equilibrium point. #TODO do we need to ignore wind states?

    Inputs:
    - state_eq: np.array, equilibrium state
    - input_eq: np.array, equilibrium input
    - dt: float, timestep

    Returns:
    - A: np.array, system matrix around equilibrium
    - B: np.array, input matrix around equilibrium
    """
    # Number of states and inputs
    n_wind = 3
    n_states = len(state_eq)
    n_inputs = len(input_eq)
    
    # Initialize matrices
    A = np.zeros((n_states, n_states))
    B = np.zeros((n_states, n_inputs))
    
    # Approximation of derivatives
    epsilon = 1e-6
    
    # Compute A matrix
    for i in range(n_states-n_wind):
        delta = np.zeros(n_states)
        delta[i] = epsilon
        state_plus = state_eq + delta
        state_minus = state_eq - delta
        
        f_plus = quadrotor_dynamics(state_plus, input_eq, dt)
        f_minus = quadrotor_dynamics(state_minus, input_eq, dt)
        
        A[:, i] = (f_plus - f_minus) / (2 * epsilon)
    
    # Compute B matrix
    for i in range(n_inputs):
        delta = np.zeros(n_inputs)
        delta[i] = epsilon
        input_plus = input_eq + delta
        input_minus = input_eq - delta
        
        f_plus = quadrotor_dynamics(state_eq, input_plus, dt)
        f_minus = quadrotor_dynamics(state_eq, input_minus, dt)
        
        B[:, i] = (f_plus - f_minus) / (2 * epsilon)
    
    return A, B


def kalman_filter(A, B, Q, R, x_hat, P, u, y):
    """
    Kalman Filter implementation for a single time-step (has to be linear)

    Inputs:
    x_hat: Previous state estimate
    P: Previous estimate covariance
    u: Control input
    y: Measurement at current time step

    Returns:
    
    """

    # Prediction Step
    x_hat_minus = A @ x_hat + B @ u
    P_minus = A @ P @ A.T + Q

    # Update Step
    K = P_minus @ C.T @ np.linalg.inv(C @ P_minus @ C.T + R) # Kalman gain
    x_hat_plus = x_hat_minus + K @ (y - C @ x_hat_minus - D @ u)
    P_plus = (np.eye(len(P)) - K @ C) @ P_minus

    return x_hat_plus, P_plus

# Simulation parameters
dt = 0.01  # Time step
time = np.arange(0, 1, dt)  # Simulation time
n_states = 15  # Number of states
n_inputs = 4  # Number of inputs
state_eq = np.zeros(n_states)  # Equilibrium state (hovering)
input_eq = np.array([9.81, 0, 0, 0])  # Equilibrium input (counteracting gravity)
Q = np.eye(n_states) * 0.01  # Process noise covariance (simulated)
R = np.eye(n_states) * 0.1  # Measurement noise covariance (simulated)

# Initial state and covariance
x_hat = np.zeros(n_states)
P = np.eye(n_states)

# Storage for simulation results
true_states = []
estimated_states = []

for t in time:
    # Simulate the true dynamics
    if t < 0.5:
        inputs = np.array([9.81, 0.1, 0.1, 0])  # Slight control input for demonstration
    else:
        inputs = np.array([9.81, -0.1, -0.1, 0])  # Change control input
    
    true_state = quadrotor_dynamics(x_hat, inputs, dt)
    true_states.append(true_state)
    
    # Generate noisy measurements
    measurements = true_state + np.random.normal(0, 0.1, n_states)
    
    # Linearize the dynamics around the current estimated state and input
    A, B = linearize(x_hat, inputs, dt)
    
    # Kalman filter to estimate the state
    # Assuming C is an identity matrix and D is a zero matrix
    C = np.eye(n_states)
    D = np.zeros((n_states, n_inputs))
    
    x_hat, P = kalman_filter(A, B, Q, R, x_hat, P, inputs, measurements)
    estimated_states.append(x_hat)

# Convert results to numpy arrays for plotting
true_states = np.array(true_states)
estimated_states = np.array(estimated_states)

# Plotting
print("*** Plotting States ***")
plt.figure(figsize=(15, 8))

for i, label in enumerate(['x', 'y', 'z', 'vx', 'vy', 'vz', 'phi (Roll)', 'theta (Pitch)', 'psi (Yaw)', 
                           'p (Angular Vel X)', 'q (Angular Vel Y)', 'r (Angular Vel Z)', 
                           'wz (Wind x Force)', 'wy (Wind y Force)', 'wz (Wind z Force)']):
    plt.subplot(5, 3, i+1)
    plt.plot(time, true_states[:, i], label='True')
    plt.plot(time, estimated_states[:, i], label='Estimated', linestyle='--')
    plt.title(label)
    plt.legend()

plt.tight_layout()
plt.show()
