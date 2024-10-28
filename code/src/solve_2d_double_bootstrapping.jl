# common solve functions for all 1D schemes
function solve_double_bootstrapping(ode::ODE{T, LinearScalarAdvectionEquation2D, 
    Semidiscretization2D}, method::String, callback=nothing) where {T}
    if method == "1st order FD"
        return solve_1st_order_FD_double_bootstrapping(ode, callback)
    elseif method == "Lax-Wendroff"
        return solve_Lax_Wendroff_double_bootstrapping(ode, callback)
    else
        error("No such method, currently implemented: 1st order FD, Lax-Wendroff")
    end
end

# Upwind scheme
function solve_1st_order_FD_double_bootstrapping(ode::ODE{T, LinearScalarAdvectionEquation2D, 
        Semidiscretization2D, PeriodicBoundaryConditions}, callback) where {T}

    u = deepcopy(ode.u0)
    dt = ode.dt
    t = ode.time_span[1]
    N_steps = ode.N_steps
    r_x = ode.eq.a[1]*dt/ode.semi.dx
    r_y = ode.eq.a[2]*dt/ode.semi.dy
    cc = nothing
    if u.context.backend isa OpenFHEBackend
        cc = SecureArithmetic.get_crypto_context(u.context)
    end
    mult_depth_solve = 2
    mult_depth_bootstrap = 3
    mult_depth_required = mult_depth_solve + mult_depth_bootstrap
    for i in range(0, N_steps-1)
        if level(u) + mult_depth_required > ode.mult_depth
            u.data = OpenFHE.EvalBootstrap(cc, u.data, 2, 19)
        end
        # perform step
        u = u - (u - circshift(u, (1, 0)))*r_x - (u - circshift(u, (0, 1)))*r_y
        t+=dt
        if !isnothing(callback) && i%callback.frequency == 0
            callback.run(u, t, ode, callback, i)
        end
        # run garbage collector for ciphertext simulations
        if u.context.backend isa OpenFHEBackend
            GC.gc()
        end
    end
    u
end

# Lax-Wendroff scheme
function solve_Lax_Wendroff_double_bootstrapping(ode::ODE{T, LinearScalarAdvectionEquation2D, 
        Semidiscretization2D, PeriodicBoundaryConditions}, callback) where {T}

    u = deepcopy(ode.u0)
    dt = ode.dt
    t = ode.time_span[1]
    N_steps = ode.N_steps
    r_x = ode.eq.a[1]*dt/ode.semi.dx
    r_y = ode.eq.a[2]*dt/ode.semi.dy
    cc = nothing
    if u.context.backend isa OpenFHEBackend
        cc = SecureArithmetic.get_crypto_context(u.context)
    end
    mult_depth_solve = 3
    mult_depth_bootstrap = 3
    mult_depth_required = mult_depth_solve + mult_depth_bootstrap
    for i in range(0, N_steps-1)
        if level(u) + mult_depth_required > ode.mult_depth
            u.data = OpenFHE.EvalBootstrap(cc, u.data, 2, 19)
        end
        # perform step
        u = (1 - r_x*r_x - r_y*r_y)*u + (0.5*r_x*r_x - 0.5*r_x)*circshift(u, (-1, 0)) +
            (0.5*r_x*r_x + 0.5*r_x)*circshift(u, (1, 0)) + (0.5*r_y*r_y - 0.5*r_y)*circshift(u, (0, -1)) +
            (0.5*r_y*r_y + 0.5*r_y)*circshift(u, (0, 1)) + (0.25*r_x*r_y)*(circshift(u, (-1, -1)) -
            circshift(u, (1, -1)) - circshift(u, (-1, 1)) + circshift(u, (1, 1)))
        
        t+=dt
        if !isnothing(callback) && i%callback.frequency == 0
            callback.run(u, t, ode, callback, i)
        end
        # run garbage collector for ciphertext simulations
        if u.context.backend isa OpenFHEBackend
            GC.gc()
        end
    end
    u
end
