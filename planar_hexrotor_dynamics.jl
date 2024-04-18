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

    # controls (thrust from each motor)
    w1 = u[1]
    w2 = u[2]
    w3 = u[3]
    w4 = u[4]
    w5 = u[5]
    w6 = u[6]

    F1 = max(0,kf*w1)
    F2 = max(0,kf*w2)
    F3 = max(0,kf*w3)
    F4 = max(0,kf*w4)
    F5 = max(0,kf*w5)
    F6 = max(0,kf*w6)
    F = [0., 0., F1+F2+F3+F4+F5+F6] #total rotor force in body frame

    # B matrix is mainly zeros, with a row for the motor coefficient, u is the w1-6
    M1 = km*w1
    M2 = km*w2
    M3 = km*w3
    M4 = km*w4
    M5 = km*w5
    M6 = km*w6
    
    τ = [
    L * (F2 - F5),  # Torque about x-axis: Difference of forces between rotor 2 and rotor 5
    L * ((F1 + F6) - (F3 + F4)),  # Torque about y-axis: Difference of forces between pairs of opposite rotors (1 & 6) and (3 & 4)
    L * (M1-M2+M3-M4+M5-M6)  # Torque about z-axis: Sum of forces acting on the body
    ]

    f = mass*gravity + Q*F # forces in world frame

    # xdot 
    [
        v
        f/mass
        ((1+norm(p)^2)/6) *(   I + 2*(skew(p)^2 + skew(p))/(1+norm(p)^2)   )*ω
        J\(τ - cross(ω,J*ω))
    ]
end



function rk4(model,ode,x,u,dt)
    # rk4 for discretizing
    k1 = dt*ode(model,x, u)
    k2 = dt*ode(model,x + k1/2, u)
    k3 = dt*ode(model,x + k2/2, u)
    k4 = dt*ode(model,x + k3, u)
    x + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
end





# function vis_traj!(vis, name, X; R = 0.1, color = mc.RGBA(1.0, 0.0, 0.0, 1.0))
#     # visualize a trajectory expressed with X::Vector{Vector}
#     for i = 1:(length(X)-1)
#         a = X[i][1:3]
#         b = X[i+1][1:3]
#         cyl = mc.Cylinder(mc.Point(a...), mc.Point(b...), R)
#         mc.setobject!(vis[name]["p"*string(i)], cyl, mc.MeshPhongMaterial(color=color))
#     end
#     for i = 1:length(X)
#         a = X[i][1:3]
#         sph = mc.HyperSphere(mc.Point(a...), R)
#         mc.setobject!(vis[name]["s"*string(i)], sph, mc.MeshPhongMaterial(color=color))
#     end
# end

# function animate_hexrotor(Xsim, Xref, dt)
#     # animate quadrotor, show Xref with vis_traj!, and track Xref with the green sphere
#     vis = mc.Visualizer()
#     robot_obj = mc.MeshFileGeometry(joinpath(@__DIR__,"hexrotor.obj"))
#     mc.setobject!(vis[:vic], robot_obj)

#     vis_traj!(vis, :traj, Xref; R = 0.01, color = mc.RGBA(1.0, 0.0, 0.0, 1.0))
#     target = mc.HyperSphere(mc.Point(0,0,0.0),0.1)
#     mc.setobject!(vis[:target], target, mc.MeshPhongMaterial(color = mc.RGBA(0.0,1.0,0.0,0.4)))


#     anim = mc.Animation(floor(Int,1/dt))
#     for k = 1:length(Xsim)
#         mc.atframe(anim, k) do
#             r = Xsim[k][1:3]
#             p = Xsim[k][7:9]
#             mc.settransform!(vis[:vic], mc.compose(mc.Translation(r),mc.LinearMap(1.5*(dcm_from_mrp(p)))))
#             mc.settransform!(vis[:target], mc.Translation(Xref[k][1:3]))
#         end
#     end
#     mc.setanimation!(vis, anim)

#     return (mc.render(vis))
# end



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

function animate_hexrotor(Xsim, Xref, dt)

    # animate quadrotor, show Xref with vis_traj!, and track Xref with the green sphere
    vis = mc.Visualizer()
    robot_obj = mc.MeshFileGeometry(joinpath(@__DIR__,"hexrotor_assembly_notilt.obj")) 
    

    mc.setobject!(vis[:drone][:base], robot_obj, mc.MeshPhongMaterial(color=mc.RGBA(0.6, 0.6, 1.0, 1.0)))
    mc.settransform!(vis[:drone][:base], mc.Translation([0,0,0]) ∘ mc.LinearMap(0.002 * I))

    vis_traj!(vis, :traj, Xref; R = 0.01, color = mc.RGBA(0.0, 0.0, 1.0, 1.0))
    target = mc.HyperSphere(mc.Point(0,0,0.0),0.1)
    mc.setobject!(vis[:target], target, mc.MeshPhongMaterial(color = mc.RGBA(0.0,1.0,0.0,0.4)))


    anim = mc.Animation(floor(Int,1/dt))
    for k = 1:length(Xsim)
        mc.atframe(anim, k) do
            r = Xsim[k][1:3]
            p = Xsim[k][7:9]
            mc.settransform!(vis[:drone], mc.compose(mc.Translation(r),mc.LinearMap((dcm_from_mrp(p)))))
            mc.settransform!(vis[:target], mc.Translation(Xref[k][1:3]))
        end
    end
    mc.setanimation!(vis, anim)

    return (mc.render(vis))
end