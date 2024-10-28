using CSV, DataFrames
include("../../src/secure_simulations.jl")
include("../common_part.jl")

GC.gc()
GC.gc()
println("Free space: ", Sys.free_memory()/2^20, "MB")

N = 64
mult_depth = 15
private_key, public_key, context, mult_depth_resulted = init_setup(mult_depth, N)
ciph_sin, _, _, u_sin = init_data(context, public_key, N, mult_depth_resulted-1)

# error investigation
res = deepcopy(ciph_sin)
error_bootstrap = zeros(mult_depth)
for i in range(1, mult_depth)
    bootstrap!(res)
    error_bootstrap[i] = maximum(abs.(collect(decrypt(res, private_key)).-u_sin))
end
release_context_memory()
GC.gc()


# runtime investigation
mult_depth_range = Int.(range(5, 25, 5))
time_bootstrap = zeros(length(mult_depth_range))
for i in range(1, length(mult_depth_range))
    local mult_depth = mult_depth_range[i]
    local private_key, public_key, context, mult_depth_resulted = init_setup(mult_depth, N)
    mult_depth_range[i] = mult_depth_resulted
    local ciph_sin, _, _, _ = init_data(context, public_key, N, mult_depth_resulted-1)
    # Bootstrapping
    time_records = zeros(6)
    for i in 1:6
        start = time_ns()
        bootstrap!(ciph_sin)
        finish = time_ns()
        time_records[i] = (finish-start)*10^-9
        GC.gc()
    end
    time_bootstrap[i] = sum(time_records[2:end])/5 # first run to precompile
    release_context_memory()
    GC.gc()
end

table = DataFrame(N_boots = range(1, mult_depth),
                  error = error_bootstrap,
                  precision_bits = -log2.(error_bootstrap),
                  mult_depth_range = vcat(mult_depth_range, zeros(mult_depth-length(mult_depth_range))),
                  time = vcat(time_bootstrap, zeros(mult_depth-length(mult_depth_range))),
                  N_time = ones(mult_depth)*length(mult_depth_range))

CSV.write("out/efficiency/bootstrapping.csv", table)
