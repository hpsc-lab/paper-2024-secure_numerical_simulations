using CSV, DataFrames
include("../../examples/advection1d.jl")
include("../common_part.jl")
# function to measure only time in between first and second (inclusive) bootstrapping
function solve_1d_LW_measure_1_iter(ode::ODE{T, LinearScalarAdvectionEquation1D, 
    Semidiscretization1D, PeriodicBoundaryConditions}, callback) where {T}

    u = deepcopy(ode.u0)
    dt = ode.dt
    t = ode.time_span[1]
    N_steps = ode.N_steps
    r = ode.eq.a*dt/ode.semi.dx
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
        steps += 1
    end
    finish_timer - start_timer, steps
end

function solve_1d_FD_measure_1_iter(ode::ODE{T, LinearScalarAdvectionEquation1D, 
    Semidiscretization1D, PeriodicBoundaryConditions}, callback) where {T}

    u = deepcopy(ode.u0)
    dt = ode.dt
    t = ode.time_span[1]
    N_steps = ode.N_steps
    r = ode.eq.a*dt/ode.semi.dx
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
        u = u - (u - circshift(u, 1; wrap_by = :length))*r
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

function configure_openfhe(a, timespan, N_steps, xspan, N_x, mult_depth, method)
    N_round = Int(2^ceil(log2(N_x)))
    cc, mult_depth =  generate_secure_cryptocontext(mult_depth, N_round)
    context = SecureContext(OpenFHEBackend(cc))
    public_key, private_key = generate_keys(context)
    init_multiplication!(context, private_key)
    init_bootstrapping!(context, private_key)
    if method == "Lax-Wendroff"
        init_rotation!(context, private_key, [-1, 1, N_x-1, -(N_x-1)])
    elseif method == "1st order FD"
        init_rotation!(context, private_key, [1, -(N_x-1)])
    else
        error("You need to initialize rotation for your method=", method)
    end
    ode = generate_setup(a, timespan, N_steps, xspan, N_x, context, 
                            mult_depth, sin2pi, public_key)
    return ode
end

function configure_plain(a, timespan, N_steps, xspan, N_x)
    context_unencrypted = SecureContext(Unencrypted())
    public_key, private_key = generate_keys(context_unencrypted)
    ode = generate_setup(a, timespan, N_steps, xspan, N_x, context_unencrypted, 
                         typemax(Int), sin2pi, public_key)
    return ode
end

function run_LW_openfhe(a, timespan, N_steps, span, N, mult_depth)
    ode = configure_openfhe(a, timespan, N_steps, span, N, mult_depth, "Lax-Wendroff")
    time, steps = solve_1d_LW_measure_1_iter(ode, nothing)
    return time * 1e-9 /steps
end

function run_LW_plain(a, timespan, N_steps, span, N)
    ode = configure_plain(a, timespan, N_steps, span, N)
    start = time_ns()
    _, steps = solve_1d_LW_measure_1_iter(ode, nothing)
    finish = time_ns()
    time = finish - start
    return time * 1e-9 /steps
end

function run_FD_openfhe(a, timespan, N_steps, span, N, mult_depth)
    ode = configure_openfhe(a, timespan, N_steps, span, N, mult_depth, "1st order FD")
    time, steps = solve_1d_FD_measure_1_iter(ode, nothing)
    return time * 1e-9 /steps
end

function run_FD_plain(a, timespan, N_steps, span, N)
    ode = configure_plain(a, timespan, N_steps, span, N)
    start = time_ns()
    _, steps = solve_1d_FD_measure_1_iter(ode, nothing)
    finish = time_ns()
    time = finish - start
    return time * 1e-9 /steps
end

N = 64
a = 1.0
timespan = (0.0, 1.0)
x_span = (0.0, 1.0)
N_steps_plain = 10000
# doesn't play a role, as we measure only upto the second bootstrapping
N_steps_openfhe = 128

mult_depth_range = Int.(range(5, 25, 11))

LW_time_step_openFHE = zeros(length(mult_depth_range))

for i in range(1, length(mult_depth_range))
    local time_step = run_LW_openfhe(a, timespan, N_steps_openfhe, x_span, N, mult_depth_range[i])
    LW_time_step_openFHE[i] = time_step
    release_context_memory()
    GC.gc()
    GC.gc()
    println("Free space: ", Sys.free_memory()/2^20, "MB")
end

mult_depth = 25

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

# add empty places
string_vec = Vector{Any}(undef, length(mult_depth_range)-1)
string_vec .= "-"
FD_time_step_openFHE = deepcopy(append!(string_vec, FD_time_step_openFHE))
string_vec[end] = LW_time_step_plain
LW_time_step_plain = deepcopy(string_vec)
string_vec[end] = FD_time_step_plain
FD_time_step_plain = deepcopy(string_vec)

# save CSV
table = DataFrame(mult_depth_range = mult_depth_range,
                  LW_time_step_openFHE = LW_time_step_openFHE,
                  FD_time_step_openFHE = FD_time_step_openFHE,
                  LW_time_step_plain = LW_time_step_plain,
                  FD_time_step_plain = FD_time_step_plain)

CSV.write("out/advection1d/time_per_step.csv", table)

