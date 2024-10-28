using Plots, CSV, DataFrames
include("../../examples/advection2d.jl")
include("../common_part.jl")
# setup parameters
a = (1.0, 1.0)
timespan = (0.0, 0.5)
mult_depth = 25
N_range = [32, 64, 128, 256]
CFL = 0.5

#Lax-Wendroff
L2_mean_LW_cfl = []
for N in N_range
    local dx = 1/N
    local dy = 1/N
    # cfl condition
    local dt = CFL/(a[1]/dx+a[2]/dy)
    local N_steps = Int(timespan[2]/dt)
    # ciphertext simulation
    GC.gc()
    GC.gc()
    println("Free space: ", Sys.free_memory()/2^20, "MB")
    local _, _, sol_ciph = advection2d("CIPHERTEXT", "Lax-Wendroff", N, N_steps, timespan, a, mult_depth, true)
    # exact solution
    local exact = advection2d("EXACT", "Lax-Wendroff", N, N_steps, timespan, a, mult_depth, true)
    #compute ad save errors
    push!(L2_mean_LW_cfl, L_2(sol_ciph, exact))
end
# experimental order of convergence
EOC_LW_cfl = []
for i in range(2, length(L2_mean_LW_cfl))
    push!(EOC_LW_cfl, log2(L2_mean_LW_cfl[i-1]/L2_mean_LW_cfl[i]))
end

# save in csv
push!(EOC_LW_cfl, sum(EOC_LW_cfl)/(length(N_range)-1))
table_LW = DataFrame(N = vcat(N_range, "mean"),
                    L2_mean = vcat(L2_mean_LW_cfl, "-"),
                    EOC = vcat("-", EOC_LW_cfl))

rename!(table_LW,["\$N_x=N_y\$","\$L2_{mean}\$", "\$EOC\$"])
CSV.write("out/advection2d/LW_convergence.csv", table_LW)
