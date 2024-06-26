import Convex as cvx 

""" 
`u = convex_mpc(A,B,X_ref_window,xic,xg,u_min,u_max,N_mpc)`

setup and solve the above optimization problem, returning the 
first control u_1 from the solution (should be a length nu 
Vector{Float64}).  
"""
function convex_mpc(A::Matrix, # discrete dynamics matrix A
                    B::Matrix, # discrete dynamics matrix B
                    wind_dynamics::Matrix, # discrete dynamics matrix D
                    X_ref_window::Vector{Vector{Float64}}, # reference trajectory for this window 
                    U_ref_window::Vector{Vector{Float64}}, # controls reference for this window
                    xic::Vector, # current state x 
                    xg::Vector, # goal state
                    u_min::Vector, # lower bound on u 
                    u_max::Vector, # upper bound on u 
                    N_mpc::Int64,  # length of MPC window (horizon)
                    Q, #Q cost
                    R, #Input cost
                    Qn,
                    wind_accel::Vector
                    #wind_dynamics::Matrix, # discrete dynamics of wind W
                    #wind_velocity::Vector # current wind state w
                    )::Vector{Float64} # return the first control command of the solved policy 
                    
    # get our sizes for state and control
    nx,nu = size(B)
    
    # check sizes 
    @assert size(A) == (nx, nx)
    @assert length(xic) == nx 
    @assert length(xg) == nx 
    @assert length(X_ref_window) == N_mpc 

    # variables we are solving for
    X = cvx.Variable(nx,N_mpc)
    U = cvx.Variable(nu,N_mpc-1)

    # cost function
    obj = 0
    # stage cost
    for k= 1: (N_mpc-1)
        xk= X[:,k] # state at timestep k
        uk= U[:,k] # controls at timestep k
        obj += 0.5*cvx.quadform((X[:,k]-X_ref_window[k]), Q)+ 0.5 *cvx.quadform((U[:,k]- U_ref_window[k]), R)
        #obj += 0.5*cvx.quadform((X[:,k]-X_ref_window[k]), Q)+ 0.5 *cvx.quadform(U[:,k], R)
    end

    #terminal cost
    obj += 0.5* cvx.quadform((X[:,N_mpc]- X_ref_window[N_mpc]), Qn)
    ################
    
    # create problem with objective
    prob = cvx.minimize(obj)

    # TODO: add constraints with prob.constraints += 
    #################
    # initial conditions constraint
    prob.constraints += (X[:,1]== xic)
    #prob.constraints += (X[:,N_mpc]== xg)

    #D= [0 0 0; 1 0 0; 0 0 0; 0 0 0]
    #estimated_wind2= [0;0; 0; estimated_wind;0; 0; 0;0; 0; 0]
    
    for k = 1:(N_mpc-1) 
        xk= X[:,k]
        uk= U[:,k]
        # linear dynamics constraints
        prob.constraints += (X[:,k+1] == A*xk + B*uk + wind_dynamics*wind_accel)
        #prob.constraints += (X[:,k+1] == A*xk + B*uk + estimated_wind2)
        # control bound constraints 
        prob.constraints += (uk <= u_max) 
        prob.constraints += (uk >= u_min)
    end

    # solve problem 
    cvx.solve!(prob, ECOS.Optimizer; silent_solver = true)

    # get X and U solutions 
    X = X.value
    U = U.value
    
    # return first control U 
    return U[:,1]
end
       