using CSV, DataFrames
include("../../examples/advection2d.jl")
include("../common_part.jl")
# function to measure only time in between first and second (inclusive) bootstrapping
function solve_2d_LW_measure_1_iter(ode::ODE{T, LinearScalarAdvectionEquation2D, 
    Semidiscretization2D, PeriodicBoundaryConditions}, callback) where {T}

    u = deepcopy(ode.u0)
    dt = ode.dt
    t = ode.time_span[1]
    N_steps = ode.N_steps
    r_x = ode.eq.a[1]*dt/ode.semi.dx
    r_y = ode.eq.a[2]*dt/ode.semi.dy
    mult_depth_solve = 3
    mult_depth_bootstrap = 2
    mult_depth_required = mult_depth_solve + mult_depth_bootstrap
    start_timer = 0.0
    finish_timer = 0.0
    bootstrapping_number = 1
    steps = 0
    for i in range(0, N_steps-1)
        if level(u) + mult_depth_required > ode.mult_depth
            u = bootstrap!(u)
            if bootstrapping_number == 1
                start_timer = time_ns()
                bootstrapping_number = 2
                steps = 0
            elseif bootstrapping_number == 2
                finish_timer = time_ns()
                break
            end
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
        steps += 1
    end
    finish_timer - start_timer, steps
end

function solve_2d_LW_double_bootstrapping_measure_1_iter(ode::ODE{T, LinearScalarAdvectionEquation2D, 
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
    start_timer = 0.0
    finish_timer = 0.0
    bootstrapping_number = 1
    steps = 0
    for i in range(0, N_steps-1)
        if level(u) + mult_depth_required > ode.mult_depth
            u.data = OpenFHE.EvalBootstrap(cc, u.data, 2, 19)
            if bootstrapping_number == 1
                start_timer = time_ns()
                bootstrapping_number = 2
                steps = 0
            elseif bootstrapping_number == 2
                finish_timer = time_ns()
                break
            end
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
        steps += 1
    end
    finish_timer - start_timer, steps
end

function solve_2d_FD_measure_1_iter(ode::ODE{T, LinearScalarAdvectionEquation2D, 
    Semidiscretization2D, PeriodicBoundaryConditions}, callback) where {T}

    u = deepcopy(ode.u0)
    dt = ode.dt
    t = ode.time_span[1]
    N_steps = ode.N_steps
    r_x = ode.eq.a[1]*dt/ode.semi.dx
    r_y = ode.eq.a[2]*dt/ode.semi.dy
    mult_depth_solve = 2
    mult_depth_bootstrap = 2
    mult_depth_required = mult_depth_solve + mult_depth_bootstrap
    start_timer = 0.0
    finish_timer = 0.0
    bootstrapping_number = 1
    steps = 0
    for i in range(0, N_steps-1)
        if level(u) + mult_depth_required > ode.mult_depth
            u = bootstrap!(u)
            if bootstrapping_number == 1
                start_timer = time_ns()
                bootstrapping_number = 2
                steps = 0
            elseif bootstrapping_number == 2
                finish_timer = time_ns()
                break
            end
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
        steps += 1
    end
    finish_timer - start_timer, steps
end

function configure_openfhe(a, timespan, N_steps, span, N, mult_depth, method)
    N_round = Int(2^ceil(log2(N*N)))
    cc, mult_depth =  generate_secure_cryptocontext(mult_depth, N_round)
    context = SecureContext(OpenFHEBackend(cc))
    public_key, private_key = generate_keys(context)
    init_multiplication!(context, private_key)
    init_bootstrapping!(context, private_key)
    if method == "1st order FD"
        init_matrix_rotation!(context, private_key, [(0, 1), (1, 0)], (N, N))
    elseif method == "Lax-Wendroff"
        init_matrix_rotation!(context, private_key, [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1),
                                                    (-1, 1), (1, -1), (-1, -1)], (N, N))
    else
        error("You need to initialize rotation for your method=", method)
    end
    ode = generate_setup(a, timespan, N_steps, span, N, context, 
                         mult_depth, sin2pi_sin2pi, public_key)
    return ode
end

function configure_plain(a, timespan, N_steps, span, N)
    context_unencrypted = SecureContext(Unencrypted())
    public_key, private_key = generate_keys(context_unencrypted)
    ode = generate_setup(a, timespan, N_steps, span, N, context_unencrypted, 
                         typemax(Int), sin2pi_sin2pi, public_key)
    return ode
end

function run_LW_openfhe(a, timespan, N_steps, span, N, mult_depth)
    ode = configure_openfhe(a, timespan, N_steps, span, N, mult_depth, "Lax-Wendroff")
    time, steps = solve_2d_LW_measure_1_iter(ode, nothing)
    return time * 1e-9 /steps
end

function run_LW_double_bootstrapping_openfhe(a, timespan, N_steps, span, N, mult_depth)
    ode = configure_openfhe(a, timespan, N_steps, span, N, mult_depth, "Lax-Wendroff")
    time, steps = solve_2d_LW_double_bootstrapping_measure_1_iter(ode, nothing)
    return time * 1e-9 /steps
end

function run_LW_plain(a, timespan, N_steps, span, N)
    ode = configure_plain(a, timespan, N_steps, span, N)
    start = time_ns()
    _, steps = solve_2d_LW_measure_1_iter(ode, nothing)
    finish = time_ns()
    time = finish - start
    return time * 1e-9 /steps
end

function run_FD_openfhe(a, timespan, N_steps, span, N, mult_depth)
    ode = configure_openfhe(a, timespan, N_steps, span, N, mult_depth, "1st order FD")
    time, steps = solve_2d_FD_measure_1_iter(ode, nothing)
    return time * 1e-9 /steps
end

function run_FD_plain(a, timespan, N_steps, span, N)
    ode = configure_plain(a, timespan, N_steps, span, N)
    start = time_ns()
    _, steps = solve_2d_FD_measure_1_iter(ode, nothing)
    finish = time_ns()
    time = finish - start
    return time * 1e-9 /steps
end

N = 64
a = (1.0, 1.0)
timespan = (0.0, 1.0)
x_span = (0.0, 1.0)
N_steps_plain = 1000
# doesn't play a role, as we measure only upto the second bootstrapping
N_steps_openfhe = 128

mult_depth = 25

LW_time_step_openFHE = run_LW_openfhe(a, timespan, N_steps_openfhe, x_span, N, mult_depth)
release_context_memory()
GC.gc()
GC.gc()
println("Free space: ", Sys.free_memory()/2^20, "MB")

LW_double_bootstrapping_time_step_openFHE = run_LW_double_bootstrapping_openfhe(a,
    timespan, N_steps_openfhe, x_span, N, mult_depth)
release_context_memory()
GC.gc()
GC.gc()
println("Free space: ", Sys.free_memory()/2^20, "MB")

FD_time_step_openFHE = run_FD_openfhe(a, timespan, N_steps_openfhe, x_span, N, mult_depth)
release_context_memory()
GC.gc()
GC.gc()
println("Free space: ", Sys.free_memory()/2^20, "MB")

GC.enable(false)
# plaintext first run is not relevant
_ = run_LW_plain(a, timespan, N_steps_plain, x_span, N)
LW_time_step_plain = run_LW_plain(a, timespan, N_steps_plain, x_span, N)

# plaintext first run is not relevant
_ = run_FD_plain(a, timespan, N_steps_plain, x_span, N)
FD_time_step_plain = run_FD_plain(a, timespan, N_steps_plain, x_span, N)

GC.enable(true)

# save CSV
table = DataFrame(LW_time_step_openFHE = LW_time_step_openFHE,
                  LW_double_bootstrapping_time_step_openFHE = LW_double_bootstrapping_time_step_openFHE,
                  FD_time_step_openFHE = FD_time_step_openFHE,
                  LW_time_step_plain = LW_time_step_plain,
                  FD_time_step_plain = FD_time_step_plain)

CSV.write("out/advection2d/time_per_step.csv", table)
