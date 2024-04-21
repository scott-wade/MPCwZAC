using ForwardDiff

# generating a rotation matrix from MRP
function dcm_from_mrp(p)
    p1,p2,p3 = p
    den = (p1^2 + p2^2 + p3^2 + 1)^2
    a = (4*p1^2 + 4*p2^2 + 4*p3^2 - 4)
    [
    (-((8*p2^2+8*p3^2)/den-1)*den)   (8*p1*p2 + p3*a)     (8*p1*p3 - p2*a);
    (8*p1*p2 - p3*a) (-((8*p1^2 + 8*p3^2)/den - 1)*den)   (8*p2*p3 + p1*a);
    (8*p1*p3 + p2*a)  (8*p2*p3 - p1*a)  (-((8*p1^2 + 8*p2^2)/den - 1)*den)
    ]/den
end

# skew matrix formation required for dynamics
function skew(ω::Vector)
    return [0    -ω[3]  ω[2];
            ω[3]  0    -ω[1];
           -ω[2]  ω[1]  0]
end

function propeller_to_body(α, β, λ) # rotation matrix from propeller to body frame
    # x axis- rotate propeller by angle α (in radians)
    # alpha is the same for each propeller (just rotating + and -)
    Rx= [1       0      0; 
         0  cos(α) sin(α); 
         0 -sin(α) cos(α)]
    # y axis rotation- rotate propeller by angle β (in radians)
    Ry= [cos(β) 0 -sin(β);
              0 1       0; 
         sin(β) 0  cos(β)]
    # z axis rotation- rotate frame by angle λ
    Rz= [cos(λ) sin(λ) 0;
         -sin(λ) cos(λ) 0;
          0 0 1]  

    return Rx, Ry, Rz
end


function hexrotor_dynamics(model::NamedTuple,x,u)
    # hexrotor dynamics with an MRP for attitude
    # and velocity in the world frame (not body frame)

    r = x[1:3]     # position in world frame 
    v = x[4:6]     # velocity world frame
    p = x[7:9]     # n_p_b (MRP) attitude 
    ω = x[10:12]   # angular velocity 

    Q = dcm_from_mrp(p)

    mass=model.mass
    J = model.J
    gravity= model.gravity
    L= model.L
    kf=model.kf
    km=model.km
    α= model.α
    β= model.β
    #λ= model.λ

    # controls (thrust from each motor)
    # w is the spinning velocity of the ith propeller
    w1 = u[1]
    w2 = u[2]
    w3 = u[3]
    w4 = u[4]
    w5 = u[5]
    w6 = u[6]

    ### rotational dynamics
    # rotation from each rotor to the body frame
    R_propeller2body= zeros(eltype(u), 6, 3, 3)
    rotor_position= zeros(eltype(u), 6, 3)
    for i= 1:6
        α= α * ((i+1)*(-1)) # every other propeller is rotated opposite (starting with propeller one tilted +α)
        β= β * ((i+1)*(-1)) # every other propeller is rotated opposite (starting with propeller one tilted +β)
        λ= (i-1)*(π/3) ## each propeller is 60 degrees from one another
        Rx, Ry, Rz= propeller_to_body(α, β, λ)

        # total rotation from each propeller to body is the multiplication of rotation in each direction
        R_propeller2body[i, :, :] = Rz * Rx * Ry

        # can also calculate the position of each rotor/propeller in the body frame
        # each propeller is 60 degrees from each other and L distance away from the body origin in the XY frame
        rotor_position[i, :] =  Rz * [L; 0; 0]
       
    end

    # torque due to thrust in the propeller frame
    # T1_p = [0.0; 0.0; max(0, kf*w1)]
    # T2_p = [0.0; 0.0; max(0, kf*w2)]
    # T3_p = [0.0; 0.0; max(0, kf*w3)]
    # T4_p = [0.0; 0.0; max(0, kf*w4)]
    # T5_p = [0.0; 0.0; max(0, kf*w5)]
    # T6_p = [0.0; 0.0; max(0, kf*w6)]
    # # vector of vectors of the thrust of each propeller (in the propeller frame)
    # thrust_propeller= [T1_p, T2_p, T3_p, T4_p, T5_p, T6_p]
    thrust_propeller = zeros(eltype(u), 6, 3)
    for i=1:6
        thrust_propeller[i, 3] = max(0.0, kf*u[i])
    end


    # thrust expressed in the body frame
    # sum of position of each propeller crossed with the thrust rotated from propeller to body frame
    #torque_thrust_body_= [zeros(3) for _ in 1:6]
    #torque_thrust_body= [zeros(eltype(x), 3) for _ in 1:6]
    torque_thrust_body = zeros(eltype(u), 6, 3)
    for i= 1:6
        torque_thrust_body[i, :] = cross(rotor_position[i,:], vec((R_propeller2body[i,:,:] * thrust_propeller[i,:])))
        #torque_thrust_body[i]= [ForwardDiff.Dual(x, 1.0) for x in torque_thrust_body_]
    end
    # add all the vectors in torque_thrust_body together
    torque_thrust_body= sum(torque_thrust_body, dims=1)

    # torque due to drag in the propeller frame
    D1_p = [0.0; 0.0; -1*km*w1] #multiplied by -1 on odd propellers since half the propellers rotate CW and half rotate CCW
    D2_p = [0.0; 0.0; km*w2]
    D3_p = [0.0; 0.0; -1*km*w3]
    D4_p = [0.0; 0.0; km*w4]
    D5_p = [0.0; 0.0; -1*km*w5]
    D6_p = [0.0; 0.0; km*w6]
    drag_propeller= [D1_p, D2_p, D3_p, D4_p, D5_p, D6_p]
    # counter balance of drag torques at hovering
   
    # drag expressed in the body frame
    # sum of the propeller frame drag torque * the rotation matrix from propeller to body
    #torque_drag_body= [zeros(eltype(x), 3) for _ in 1:6]
    torque_drag_body = zeros(eltype(u), 6, 3)
    for i= 1:6
        torque_drag_body[i, :] = ((R_propeller2body[i,:,:]) * drag_propeller[i])
    end
    torque_drag_body= vec(sum(torque_drag_body, dims=1))

    # total torque of each propeller in the body frame is the sum of the thrust and drag torques
    τ =  torque_thrust_body' + torque_drag_body


    ### translational dynamics
    # total force produced in the world frame
    # rotation matrix from body to world (Q) * the thrust force rotated from the propeller to body frame
    thrust_force_body= zeros(eltype(u), 6, 3)
    for i= 1:6
        thrust_force_body[i,:] = R_propeller2body[i,:,:] * thrust_propeller[i,:]
    end
    thrust_force_body= sum(thrust_force_body, dims=1)

    # total force in world frame
    f = mass*gravity + Q * thrust_force_body' 

    # xdot 
    [
        v
        f/mass
        ((1+norm(p)^2)/6) *(   I + 2*(skew(p)^2 + skew(p))/(1+norm(p)^2)   )*ω
        J\(τ - cross(ω,J*ω))
    ]
end




# All controls 0
# Then rigid body (no velocity then some velocity, want to get projectile motion
# Then set one control, check if velocity lines up
# Set controls where you know it will cancel (hover) or rotate and test with that


# Get cost to go from hover (which is from k gains of IHLQR)
# For Qn

# Linearize around hover

# ^ if constraints are inactive

# If the constraints are active, horizon has to be long enough to not violate constraints before switching to LQR




function rk4(model,ode,x,u,dt)
    # rk4 for discretizing
    k1 = dt*ode(model,x, u)
    k2 = dt*ode(model,x + k1/2, u)
    k3 = dt*ode(model,x + k2/2, u)
    k4 = dt*ode(model,x + k3, u)
    result = x + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    reshape(result, length(result))
end






function vis_traj!(vis, name, X; R = 0.1, color = mc.RGBA(1.0, 0.0, 0.0, 1.0))
    # visualize a trajectory expressed with X::Vector{Vector}
    for i = 1:(length(X)-1)
        a = X[i][1:3]
        b = X[i+1][1:3]
        cyl = mc.Cylinder(mc.Point(a...), mc.Point(b...), R)
        mc.setobject!(vis[name]["p"*string(i)], cyl, mc.MeshPhongMaterial(color=color))
    end
    for i = 1:length(X)
        a = X[i][1:3]
        sph = mc.HyperSphere(mc.Point(a...), R)
        mc.setobject!(vis[name]["s"*string(i)], sph, mc.MeshPhongMaterial(color=color))
    end
end

# function animate_hexrotor(Xsim, Xref, dt)

#     # animate quadrotor, show Xref with vis_traj!, and track Xref with the green sphere
#     vis = mc.Visualizer()
#     robot_obj = mc.MeshFileGeometry(joinpath(@__DIR__,"hexrotor_assembly_notilt.obj")) 
    

#     mc.setobject!(vis[:drone][:base], robot_obj, mc.MeshPhongMaterial(color=mc.RGBA(0.6, 0.6, 1.0, 1.0)))
#     mc.settransform!(vis[:drone][:base], mc.Translation([0,0,0]) ∘ mc.LinearMap(0.002 * I))

#     vis_traj!(vis, :traj, Xref; R = 0.01, color = mc.RGBA(0.0, 0.0, 1.0, 1.0))
#     target = mc.HyperSphere(mc.Point(0,0,0.0),0.1)
#     mc.setobject!(vis[:target], target, mc.MeshPhongMaterial(color = mc.RGBA(0.0,1.0,0.0,0.4)))


#     anim = mc.Animation(floor(Int,1/dt))
#     for k = 1:length(Xsim)
#         mc.atframe(anim, k) do
#             r = Xsim[k][1:3]
#             p = Xsim[k][7:9]
#             mc.settransform!(vis[:drone], mc.compose(mc.Translation(r),mc.LinearMap((dcm_from_mrp(p)))))
#             mc.settransform!(vis[:target], mc.Translation(Xref[k][1:3]))
#         end
#     end
#     mc.setanimation!(vis, anim)

#     return (mc.render(vis))
# end

function animate_hexrotor(Xsim, dt)
    vis = mc.Visualizer()
    #robot_obj = mc.MeshFileGeometry(joinpath(@__DIR__,"hexrotor_assembly_notilt.obj"))
    robot_obj = mc.MeshFileGeometry(joinpath(@__DIR__,"utils/quadrotor.obj"))

    mc.setobject!(vis[:drone], robot_obj)

    mc.setobject!(vis[:drone][:base], robot_obj, mc.MeshPhongMaterial(color=mc.RGBA(0.6, 0.6, 1.0, 1.0)))
    mc.settransform!(vis[:drone][:base], mc.Translation([0,0,0]) ∘ mc.LinearMap(0.001 * I))

    # vis_traj!(vis, :traj, Xref[1:85]; R = 0.01, color = mc.RGBA(1.0, 0.0, 0.0, 1.0))
    # target = mc.HyperSphere(mc.Point(0,0,0.0),0.1)
    # mc.setobject!(vis[:target], target, mc.MeshPhongMaterial(color = mc.RGBA(0.0,1.0,0.0,0.4)))


    anim = mc.Animation(floor(Int,1/dt))
    for k = 1:length(Xsim)
        mc.atframe(anim, k) do
            r = Xsim[k][1:3]
            p = Xsim[k][7:9]
            mc.settransform!(vis[:drone], mc.compose(mc.Translation(r),mc.LinearMap(1.5*(dcm_from_mrp(p)))))
            #mc.settransform!(vis[:target], mc.Translation(Xref[k][1:3]))
        end
    end
    mc.setanimation!(vis, anim)

    return (mc.render(vis))
end





##### TO DO #####
#   - model parameters: coefficients, mass, inertia, length, effect
#   - reference trajectories
#   - figure out how to test dynamics model
#   - mathematically define constraints 
        # - motor torque limits (from datasheet- save reference)
        # - angle limits
        # - input torques greater than 0
        # - α and β limits? (0 to π/2 like paper?)
#   - visualization mesh cat stuff