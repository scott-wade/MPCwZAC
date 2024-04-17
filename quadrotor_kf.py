import numpy as np
import matplotlib.pyplot as plt

def quadrotor_dynamics(state, inputs, dt):
    """
    Example quadrotor dynamics without wind model uncertainty.

    Inputs:
    - state: np.array, current state [x, y, z, vx, vy, vz, phi, theta, psi, p, q, r]
    - inputs: np.array, inputs [thrust, tau_phi, tau_theta, tau_psi]
    - dt: float, timestep

    Returns:
    - new_state: np.array, next state after dt
    """
    # Unpack the state
    x, y, z, vx, vy, vz, phi, theta, psi, p, q, r = state
    
    # Unpack the inputs
    thrust, tau_phi, tau_theta, tau_psi = inputs
    
    # Constants
    m = 1.0  # Mass of the quadrotor
    g = 9.81  # Gravitational acceleration
    Ix = 0.01  # Moment of inertia around x-axis
    Iy = 0.01  # Moment of inertia around y-axis
    Iz = 0.02  # Moment of inertia around z-axis
    
    # Translational dynamics
    ax = 0
    ay = 0
    az = thrust / m - g
    
    # Rotational dynamics
    p_dot = tau_phi / Ix
    q_dot = tau_theta / Iy
    r_dot = tau_psi / Iz
    
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
        r + r_dot * dt
    ])
    
    return new_state
    

def linearize(state_eq, input_eq, dt):
    """
    Linearize the quadrotor dynamics around a given equilibrium point.

    Inputs:
    - state_eq: np.array, equilibrium state
    - input_eq: np.array, equilibrium input
    - dt: float, timestep

    Returns:
    - A: np.array, system matrix around equilibrium
    - B: np.array, input matrix around equilibrium
    """
    # Number of states and inputs
    n_states = len(state_eq)
    n_inputs = len(input_eq)
    
    # Initialize matrices
    A = np.zeros((n_states, n_states))
    B = np.zeros((n_states, n_inputs))
    
    # Approximation of derivatives
    epsilon = 1e-6
    
    # Compute A matrix
    for i in range(n_states):
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


def kalman_filter(A, B, B_d, Q, R, x_hat, P, u, y, wind_mean, wind_covariance):
    """
    Kalman Filter implementation for a single time-step with disturbance (wind) model

    Inputs:
    - A: System matrix
    - B: Input matrix
    - B_d: Disturbance effect matrix
    - Q: Process noise covariance
    - R: Measurement noise covariance
    - x_hat: Previous state estimate
    - P: Previous estimate covariance
    - u: Control input
    - y: Measurement at current time step
    - wind_mean: Mean of the wind disturbance
    - wind_covariance: Covariance of the wind disturbance

    Returns:
    - x_hat_plus: Updated state estimate
    - P_plus: Updated estimate covariance
    """
    # Generate disturbance (wind) sample
    d_t = np.random.multivariate_normal(wind_mean, wind_covariance)

    # Prediction Step
    x_hat_minus = A @ x_hat + B @ u + B_d @ d_t
    P_minus = A @ P @ A.T + Q

    # Update Step
    C = np.eye(len(x_hat))  # Assuming we are measuring all states
    K = P_minus @ C.T @ np.linalg.inv(C @ P_minus @ C.T + R)  # Kalman gain
    x_hat_plus = x_hat_minus + K @ (y - C @ x_hat_minus)
    P_plus = (np.eye(len(P)) - K @ C) @ P_minus

    return x_hat_plus, P_plus


# Simulation parameters
dt = 0.01  # Time step
total_time = 10.0  # Total simulation time
time = np.arange(0, total_time, dt)  # Simulation time array
n_states = 12  # Number of states
n_inputs = 4  # Number of inputs

# Equilibrium state and input (hovering)
state_eq = np.zeros(n_states)
input_eq = np.array([9.81, 0, 0, 0])  # Input to counteract gravity

# Disturbance (wind) model parameters
wind_mean = np.zeros(n_states)  # Assume wind has zero mean
wind_covariance = np.diag([0.01] * 3 + [0] * 9)  # Variance in the translational states (x, y, z)

# System matrices (would be computed from the system's dynamics)
A, B = linearize(state_eq, input_eq, dt)
B_d = np.zeros((n_states, n_states))  # Effect of disturbance on the state
B_d[:3, :3] = np.eye(3)  # Assuming wind affects only the translational states

# Noise covariances
Q = np.eye(n_states) * 0.01  # Process noise covariance
R = np.eye(n_states) * 0.1  # Measurement noise covariance

# Initial state and covariance
x_hat = np.zeros(n_states)  # Initial estimate
P = np.eye(n_states)  # Initial estimate covariance

# Storage for simulation results
true_states = []
estimated_states = []

for t in time:
    # Simulate true dynamics with random wind disturbance
    wind_disturbance = np.random.multivariate_normal(wind_mean, wind_covariance)
    true_state = quadrotor_dynamics(x_hat, input_eq, dt) + B_d @ wind_disturbance
    true_states.append(true_state)
    
    # Generate noisy measurements
    measurements = true_state + np.random.normal(0, 0.1, n_states)
    
    # Kalman filter to estimate the state
    x_hat, P = kalman_filter(A, B, B_d, Q, R, x_hat, P, input_eq, measurements, wind_mean, wind_covariance)
    estimated_states.append(x_hat)

# Convert results to numpy arrays for plotting
true_states = np.array(true_states)
estimated_states = np.array(estimated_states)

# Plotting
print("*** Plotting States ***")
plt.figure(figsize=(12, 8))

for i, label in enumerate(['x', 'y', 'z', 'vx', 'vy', 'vz', 'phi (Roll)', 'theta (Pitch)', 'psi (Yaw)', 'p (Angular Vel X)', 'q (Angular Vel Y)', 'r (Angular Vel Z)']):
    plt.subplot(4, 3, i+1)
    plt.plot(time, true_states[:, i], label='True')
    plt.plot(time, estimated_states[:, i], label='Estimated', linestyle='--')
    plt.xlabel('Time (s)')
    plt.ylabel(label)
    plt.legend()

plt.tight_layout()
plt.show()