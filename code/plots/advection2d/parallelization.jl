using CSV, DataFrames
include("../../examples/advection2d.jl")
include("../common_part.jl")

# setup parameters
N = 64
a = (1.0, 1.0)
CFL = 0.5
mult_depth = 25
timespan = (0.0, 1.0)
# cfl condition
dx = 1/N
dy = 1/N
dt = CFL/(a[1]/dx+a[2]/dy)
N_steps = Int(timespan[2]/dt)

#Lax-Wendroff
# ciphertext simulations
timer_ciph_conf_LW, timer_ciph_comp_LW, _ = advection2d("CIPHERTEXT", "Lax-Wendroff", N, N_steps,
                                                        timespan, a, mult_depth, true)

table = DataFrame(OMP_THREADS = [ENV["OMP_NUM_THREADS"]],
                  conf_time = [timer_ciph_conf_LW],
                  comp_time = [timer_ciph_comp_LW])

if !isfile("out/advection2d/parallelization.csv")
    CSV.write("out/advection2d/parallelization.csv", table)
else
    CSV.write("out/advection2d/parallelization.csv", table, append=true)
end
