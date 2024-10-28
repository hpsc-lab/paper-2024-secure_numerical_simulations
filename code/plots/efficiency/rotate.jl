using CSV, DataFrames
include("../../src/secure_simulations.jl")
include("../common_part.jl")

GC.gc()
GC.gc()
println("Free space: ", Sys.free_memory()/2^20, "MB")

# error investigation
N = 64
mult_depth = 15
private_key, public_key, context, _ = init_setup(mult_depth, N)
ciph_sin, _, _, u_sin = init_data(context, public_key, N)

# rotate ciphertext
res_1 = deepcopy(ciph_sin)
res_5 = deepcopy(ciph_sin)
res_25 = deepcopy(ciph_sin)
exact_1 = deepcopy(u_sin)
exact_5 = deepcopy(u_sin)
exact_25 = deepcopy(u_sin)
error_rotate_1 = zeros(mult_depth)
error_rotate_5 = zeros(mult_depth)
error_rotate_25 = zeros(mult_depth)
for i in range(1, mult_depth)
    global res_1 = circshift(res_1, 1)
    global exact_1 = circshift(exact_1, 1)
    error_rotate_1[i] = maximum(abs.(exact_1.-collect(decrypt(res_1, private_key))))
    global res_5 = circshift(res_5, -5)
    global exact_5 = circshift(exact_5, -5)
    error_rotate_5[i] = maximum(abs.(exact_5.-collect(decrypt(res_5, private_key))))
    global res_25 = circshift(res_25, 25)
    global exact_25 = circshift(exact_25, 25)
    error_rotate_25[i] = maximum(abs.(exact_25.-collect(decrypt(res_25, private_key))))
    GC.gc()
end
release_context_memory()
GC.gc()


# runtime investigation
mult_depth_range = Int.(range(5, 25, 5))
time_rotate_1 = zeros(length(mult_depth_range))
time_rotate_5 = zeros(length(mult_depth_range))
time_rotate_25 = zeros(length(mult_depth_range))
for i in range(1, length(mult_depth_range))
    local mult_depth = mult_depth_range[i]
    local private_key, public_key, context, mult_depth_resulted = init_setup(mult_depth, N)
    mult_depth_range[i] = mult_depth_resulted
    local ciph_sin, _, _, _ = init_data(context, public_key, N)
    # rotate ciphertext
    time_records = zeros(6)
    for i in 1:6
        start = time_ns()
        circshift(ciph_sin, 1)
        finish = time_ns()
        time_records[i] = (finish-start)*10^-9
        GC.gc()
    end
    time_rotate_1[i] = sum(time_records[2:end])/5 # first run to precompile
    
    for i in 1:6
        start = time_ns()
        circshift(ciph_sin, -5)
        finish = time_ns()
        time_records[i] = (finish-start)*10^-9
        GC.gc()
    end
    time_rotate_5[i] = sum(time_records[2:end])/5 # first run to precompile
    
    for i in 1:6
        start = time_ns()
        circshift(ciph_sin, 25)
        finish = time_ns()
        time_records[i] = (finish-start)*10^-9
        GC.gc()
    end
    time_rotate_25[i] = sum(time_records[2:end])/5 # first run to precompile
    
    release_context_memory()
    GC.gc()
end

table = DataFrame(N_rotates = range(1, mult_depth),
                  error_rotate_1 = error_rotate_1,
                  error_rotate_5 = error_rotate_5,
                  error_rotate_25 = error_rotate_25,
                  mult_depth_range = vcat(mult_depth_range, zeros(mult_depth-length(mult_depth_range))),
                  time_rotate_1 = vcat(time_rotate_1, zeros(mult_depth-length(mult_depth_range))),
                  time_rotate_5 = vcat(time_rotate_5, zeros(mult_depth-length(mult_depth_range))),
                  time_rotate_25 = vcat(time_rotate_25, zeros(mult_depth-length(mult_depth_range))),
                  N_time = ones(mult_depth)*length(mult_depth_range))

CSV.write("out/efficiency/rotate.csv", table)
