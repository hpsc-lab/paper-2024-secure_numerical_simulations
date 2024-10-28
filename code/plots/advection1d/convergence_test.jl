using Plots, CSV, DataFrames
include("../../examples/advection1d.jl")
include("../common_part.jl")
# setup parameters
a = 1.0
mult_depth = 25
timespan = (0, 0.5)
N_range = [32, 64, 128, 256]
CFL = 0.5

#Lax-Wendroff
L2_mean_LW_cfl = []
for N in N_range
    local dx = 1/N
    # cfl condition
    local dt = CFL/(a[1]/dx)
    local N_steps = Int(timespan[2]/dt)
    # ciphertext simulation
    GC.gc()
    GC.gc()
    println("Free space: ", Sys.free_memory()/2^20, "MB")
    local _, _, sol_ciph = advection1d("CIPHERTEXT", "Lax-Wendroff", N, N_steps, timespan, a, mult_depth, true)
    # exact
    local exact = advection1d("EXACT", "Lax-Wendroff", N, N_steps, timespan, a, mult_depth, true)
    #compute and save errors
    push!(L2_mean_LW_cfl, L_2(sol_ciph, exact))
    GC.gc()
end
# experimental order of convergence
EOC_LW_cfl = []
for i in range(2, length(L2_mean_LW_cfl))
    push!(EOC_LW_cfl, log2(L2_mean_LW_cfl[i-1]/L2_mean_LW_cfl[i]))
end

#1st order finite differences
L2_mean_FD_cfl = []
for N in N_range
    local dx = 1/N
    # cfl condition
    local dt = CFL/(a[1]/dx)
    local N_steps = Int(timespan[2]/dt)
    # ciphertext simulation
    GC.gc()
    GC.gc()
    println("Free space: ", Sys.free_memory()/2^20, "MB")
    local _, _, sol_ciph = advection1d("CIPHERTEXT", "1st order FD", N, N_steps, timespan, a, mult_depth, true)
    # exact
    local exact = advection1d("EXACT", "1st order FD", N, N_steps, timespan, a, mult_depth, true)
    #compute and save errors
    push!(L2_mean_FD_cfl, L_2(sol_ciph, exact))
end
# experimental order of convergence
EOC_FD_cfl = []
for i in range(2, length(L2_mean_FD_cfl))
    push!(EOC_FD_cfl, log2(L2_mean_FD_cfl[i-1]/L2_mean_FD_cfl[i]))
end

# save data in csv
push!(EOC_LW_cfl, sum(EOC_LW_cfl)/(length(N_range)-1))
table_LW = DataFrame(N = vcat(N_range, "mean"),
                    Error = vcat(L2_mean_LW_cfl, "-"),
                    EOC = vcat("-", EOC_LW_cfl))

rename!(table_LW,["\$N_x\$","\$L2_{mean}\$", "\$EOC\$"])
CSV.write("out/advection1d/LW_convergence.csv", table_LW)
push!(EOC_FD_cfl, sum(EOC_FD_cfl)/(length(N_range)-1))
table_FD = DataFrame(N = vcat(N_range, "mean"),
                     Error = vcat(L2_mean_FD_cfl, "-"),
                     EOC = vcat("-", EOC_FD_cfl))

rename!(table_FD,["\$N_x\$","\$L2_{mean}\$", "\$EOC\$"])
CSV.write("out/advection1d/FD_convergence.csv", table_FD)
