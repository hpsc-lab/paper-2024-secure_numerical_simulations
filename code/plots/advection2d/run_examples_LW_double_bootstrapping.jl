using CSV, DataFrames, HDF5
include("../../examples/advection2d_double_bootstrapping.jl")
include("../common_part.jl")
# setup parameters
N = 64
a = (1.0, 1.0)
CFL = 0.5
mult_depth = 25
timespan = (0.0, 1.0)
frequency = 1
# cfl condition
dx = 1/N
dy = 1/N
dt = CFL/(a[1]/dx+a[2]/dy)
N_steps = Int(timespan[2]/dt)

#Lax-Wendroff
GC.gc()
GC.gc()
println("Free space: ", Sys.free_memory()/2^20, "MB")
# plaintext simulation
plain_result_LW = Array{Float64, 3}(undef, N_steps, N, N)
callback_record_run(u, t, ode, callback, i) = callback_record_2d(u, t, ode, callback, i, plain_result_LW)
callback = StandartCallback(L_inf, nothing, nothing, callback_record_run, frequency)
# disable GC for plaintext simulations
GC.enable(false)
# first run is not representative in sense of time
advection2d("PLAINTEXT", "Lax-Wendroff", N, N_steps, timespan, a, mult_depth, true, callback)
plain_result_LW = Array{Float64, 3}(undef, N_steps, N, N)
callback_record_run(u, t, ode, callback, i) = callback_record_2d(u, t, ode, callback, i, plain_result_LW)
callback = StandartCallback(L_inf, nothing, nothing, callback_record_run, frequency)
# second run
timer_plain_conf_vec_LW, _, _ = advection2d("PLAINTEXT", "Lax-Wendroff", N, N_steps, timespan,
                                            a, mult_depth, true, callback)
timer_plain_comp_vec_LW = callback.time_step .- callback.time_step[1]
# ciphertext simulation
# enable GC for ciphertext simulations
GC.enable(true)
cipher_result_LW = Array{Float64, 3}(undef, N_steps, N, N)
callback_record_run(u, t, ode, callback, i) = callback_record_2d(u, t, ode, callback, i, cipher_result_LW)
callback = StandartCallback(L_inf, nothing, nothing, callback_record_run, frequency)
timer_ciph_conf_vec_LW, _, sol_ciph_LW = advection2d("CIPHERTEXT", "Lax-Wendroff", N, N_steps,
                                                     timespan, a, mult_depth, true, callback)
timer_ciph_comp_vec_LW = callback.time_step
L_inf_cipher_plain_LW = Vector{Float64}(undef, N_steps)
for i in range(1, N_steps)
    L_inf_cipher_plain_LW[i] = L_inf(cipher_result_LW[i,:,:], plain_result_LW[i,:,:])
end

# get exact solution
sol_exact = advection2d("EXACT", "Lax-Wendroff", N, N_steps, timespan, a, mult_depth, true)
# full exact solution
sol_exact_full = Array{Float64, 3}(undef, N_steps, N, N)
for i in range(1, N_steps)
    sol_exact_full[i,:,:] = solution_linear_advection_equation2D(a, (0.0, dt*i), (0.0, 1.0), N, sin2pi_sin2pi)
end


# write out results
# solutions hdf5
h5open("out/advection2d/full_solutions_LW_double_bootstrapping.h5", "w") do file
    file["LW_plain"] = plain_result_LW
    file["LW_cipher"] = cipher_result_LW
    file["exact"] = sol_exact_full
end

steps = callback.steps
# error
table_error = DataFrame(steps = steps,
                        L_inf_cipher_plain_LW = L_inf_cipher_plain_LW)

CSV.write("out/advection2d/error_LW_double_bootstrapping.csv", table_error)
# time
timer_ciph_conf_vec_LW = timer_ciph_conf_vec_LW*ones(length(steps))
timer_plain_conf_vec_LW = timer_plain_conf_vec_LW*ones(length(steps))
table_time = DataFrame(steps = steps,
                       timer_ciph_comp_vec_LW = timer_ciph_comp_vec_LW,
                       timer_ciph_conf_vec_LW = timer_ciph_conf_vec_LW,
                       timer_plain_comp_vec_LW = timer_plain_comp_vec_LW,
                       timer_plain_conf_vec_LW = timer_plain_conf_vec_LW)

CSV.write("out/advection2d/time_LW_double_bootstrapping.csv", table_time)
# visualization
table_visualization = DataFrame(length = N*ones(N*N),
                                sol_ciph_LW = vec(sol_ciph_LW),
                                sol_exact = vec(sol_exact))

CSV.write("out/advection2d/visualization_LW_double_bootstrapping.csv", table_visualization)
