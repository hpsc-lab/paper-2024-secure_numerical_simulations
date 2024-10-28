using CSV, DataFrames, HDF5
include("../../examples/advection2d.jl")
include("../common_part.jl")
# setup parameters
N = 64
a = (1.0, 1.0)
CFL = 0.5
mult_depth_range = Int.(range(5, 25, 11))
timespan = (0.0, 1.0)
frequency = 1
# cfl condition
dx = 1/N
dy = 1/N
dt = CFL/(a[1]/dx+a[2]/dy)
N_steps = Int(timespan[2]/dt)

#Lax-Wendroff
# disable GC for plaintext simulations
GC.enable(false)
_, _, plain_sol = advection2d("PLAINTEXT", "Lax-Wendroff", N, N_steps, timespan, a, mult_depth_range[1], true)
# ciphertext simulation
# enable GC for ciphertext simulations
GC.enable(true)
cipher_result_LW_vec = Array{Float64, 4}(undef, length(mult_depth_range), N_steps, N, N)
LW_conf_time_vec = Vector{Float64}(undef, length(mult_depth_range))
LW_comp_time_vec = Array{Float64, 2}(undef, length(mult_depth_range), N_steps)
LW_time_per_step_vec = Vector{Float64}(undef, length(mult_depth_range))
LW_error_vec = Vector{Float64}(undef, length(mult_depth_range))
LW_error_per_step_vec = Vector{Float64}(undef, length(mult_depth_range))
for i in range(1, length(mult_depth_range))
    GC.gc()
    GC.gc()
    println("Free space: ", Sys.free_memory()/2^20, "MB")
    local cipher_result_LW = Array{Float64, 3}(undef, N_steps, N, N)
    callback_record_run(u, t, ode, callback, i) = callback_record_2d(u, t, ode, callback, i, cipher_result_LW)
    cipher_result_LW_vec[i,:,:,:] .= cipher_result_LW
    local callback = StandartCallback(L_inf, nothing, nothing, callback_record_run, frequency)
    local LW_conf_time, LW_comp_time, sol_ciph_LW = advection2d("CIPHERTEXT", "Lax-Wendroff", N, N_steps,
                                                                timespan, a, mult_depth_range[i], true, callback)
    LW_conf_time_vec[i] = LW_conf_time
    LW_comp_time_vec[i,:] .= callback.time_step
    LW_time_per_step_vec[i] = LW_comp_time/N_steps
    LW_error_vec[i] = L_inf(plain_sol, sol_ciph_LW)
    LW_error_per_step_vec[i] = L_inf(plain_sol, sol_ciph_LW)/N_steps
end


# write out results
# solutions hdf5
h5open("out/advection2d/LW_diff_mult_depth.h5", "w") do file
    file["mult_depth_range"] = mult_depth_range
    file["cipher_result_LW_vec"] = cipher_result_LW_vec
    file["LW_conf_time_vec"] = LW_conf_time_vec
    file["LW_comp_time_vec"] = LW_comp_time_vec
    file["LW_time_per_step_vec"] = LW_time_per_step_vec
    file["LW_error_vec"] = LW_error_vec
    file["LW_error_per_step_vec"] = LW_error_per_step_vec
end
# CSV file to plot
table = DataFrame(mult_depth_range = mult_depth_range,
                  LW_error_vec = LW_error_vec,
                  LW_comp_time_vec = LW_comp_time_vec[:, end])

CSV.write("out/advection2d/LW_diff_mult_depth.csv", table)