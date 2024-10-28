# struct to collect data to solve ODEs
struct ODE{Data, Equation, Semidiscretization, BoundaryConditions <: AbstractBoundaryConditions}
    u0::Data
    eq::Equation
    semi::Semidiscretization
    time_span::Tuple
    dt::Float64
    N_steps::Int64
    mult_depth::Int64
    boundary_conditions::BoundaryConditions
    function ODE(u0::Data, eq::Equation, semi::Semidiscretization, time_span::Tuple, N_steps::Int64,
                 mult_depth::Int64, boundary_conditions::BoundaryConditions = PeriodicBoundaryConditions()
                ) where {Data, Equation, Semidiscretization, BoundaryConditions}
        dt = (time_span[2] - time_span[1])/N_steps
        new{Data, Equation, Semidiscretization, BoundaryConditions}(u0, eq, semi, time_span,
                                                                    dt, N_steps, mult_depth,
                                                                    boundary_conditions)
    end
end

# initial conditions 1D
function set_initial_conditions!(u, t, init_cond::Function, semi::Semidiscretization1D)
        for i in range(1, semi.size)
            u[i] = init_cond(semi.nodes[i], t)
        end
end

# initial conditions 2D
function set_initial_conditions!(u, t, init_cond::Function, semi::Semidiscretization2D)
        for i in range(1, semi.size_x)
            for j in range(1, semi.size_y)
                u[i, j] = init_cond([semi.nodes_x[i], semi.nodes_y[j]], t)
            end
        end
end