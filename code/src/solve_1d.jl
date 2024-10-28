# common solve functions for all 1D schemes
function solve(ode::ODE{T, LinearScalarAdvectionEquation1D, 
    Semidiscretization1D}, method::String, callback=nothing) where {T}
    if method == "1st order FD"
        return solve_1st_order_FD(ode, callback)
    elseif method == "Lax-Wendroff"
        return solve_Lax_Wendroff(ode, callback)
    else
        error("No such method, currently implemented: 1st order FD; Lax-Wendroff;")
    end
end

# Upwind scheme
function solve_1st_order_FD(ode::ODE{T, LinearScalarAdvectionEquation1D, 
        Semidiscretization1D, PeriodicBoundaryConditions}, callback) where {T}

    u = deepcopy(ode.u0)
    dt = ode.dt
    t = ode.time_span[1]
    N_steps = ode.N_steps
    r = ode.eq.a*dt/ode.semi.dx
    mult_depth_solve = 2
    mult_depth_bootstrap = 2
    mult_depth_required = mult_depth_solve + mult_depth_bootstrap
    for i in range(0, N_steps-1)
        if level(u) + mult_depth_required > ode.mult_depth
            u = bootstrap!(u)
        end
        # perform step
        u = u - (u - circshift(u, 1; wrap_by = :length))*r
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
function solve_Lax_Wendroff(ode::ODE{T, LinearScalarAdvectionEquation1D, 
    Semidiscretization1D, PeriodicBoundaryConditions}, callback) where {T}

    u = deepcopy(ode.u0)
    dt = ode.dt
    t = ode.time_span[1]
    N_steps = ode.N_steps
    r = ode.eq.a*dt/ode.semi.dx
    mult_depth_solve = 2
    mult_depth_bootstrap = 2
    mult_depth_required = mult_depth_solve + mult_depth_bootstrap
    for i in range(0, N_steps-1)
        if level(u) + mult_depth_required > ode.mult_depth
            u = bootstrap!(u)
        end
        # perform step
        u = u*(1-r*r) - circshift(u, -1; wrap_by = :length)*(r/2-r*r/2) + 
            circshift(u, 1; wrap_by = :length)*(r/2+r*r/2)
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
