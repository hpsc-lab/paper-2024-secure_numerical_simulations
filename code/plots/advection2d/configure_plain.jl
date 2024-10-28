using CSV, DataFrames
include("../../examples/advection2d.jl")
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
# plaintext simulation with callback to write out solution at every step
# disable GC for plaintext simulations
GC.enable(false)
# first run is not representative in sense of time
advection2d("PLAINTEXT", "Lax-Wendroff", N, N_steps, timespan, a, mult_depth, true)
# second run
timer_conf_LW_vec = []
for i in range(1, 100)
    local timer_conf, _, _ = advection2d("PLAINTEXT", "Lax-Wendroff", N, N_steps, timespan,
                                         a, mult_depth, true)
    append!(timer_conf_LW_vec, timer_conf)
end
# compute mean value
append!(timer_conf_LW_vec, sum(timer_conf_LW_vec)/length(timer_conf_LW_vec))
GC.enable(true)


#1st order finite differences
GC.gc()
GC.gc()
println("Free space: ", Sys.free_memory()/2^20, "MB")
# plaintext simulation with callback to write out solution at every step
# disable GC for plaintext simulations
GC.enable(false)
# first run is not representative in sense of time
advection2d("PLAINTEXT", "1st order FD", N, N_steps, timespan, a, mult_depth, true)
# second run
timer_conf_FD_vec = []
for i in range(1, 100)
    local timer_conf, _, _ = advection2d("PLAINTEXT", "1st order FD", N, N_steps, timespan,
                                         a, mult_depth, true)
    append!(timer_conf_FD_vec, timer_conf)
end
# compute mean value
append!(timer_conf_FD_vec, sum(timer_conf_FD_vec)/length(timer_conf_FD_vec))
GC.enable(true)

# write out results
table_time = DataFrame(timer_conf_LW_vec = timer_conf_LW_vec,
                       timer_conf_FD_vec = timer_conf_FD_vec)

CSV.write("out/advection2d/conf_time_plain.csv", table_time)
