using CSV, DataFrames
include("../../src/secure_simulations.jl")
include("../common_part.jl")

GC.gc()
GC.gc()
println("Free space: ", Sys.free_memory()/2^20, "MB")

N = 64
mult_depth = 15
private_key, public_key, context, _ = init_setup(mult_depth, N)
ciph_sin, _, _, u_sin = init_data(context, public_key, N)

# error investigation
exact = deepcopy(u_sin)
error_encrypt_decrypt = zeros(mult_depth)
for i in range(1, mult_depth)
    ciph = encrypt(u_sin, public_key, context)
    global u_sin = collect(decrypt(ciph, private_key))
    error_encrypt_decrypt[i] = maximum(abs.(exact.-u_sin))
end
release_context_memory()
GC.gc()

# runtime investigation
mult_depth_range = Int.(range(5, 25, 5))
time_encrypt = zeros(length(mult_depth_range))
time_decrypt = zeros(length(mult_depth_range))
for i in range(1, length(mult_depth_range))
    local mult_depth = mult_depth_range[i]
    local private_key, public_key, context, mult_depth_resulted = init_setup(mult_depth, N)
    mult_depth_range[i] = mult_depth_resulted
    local ciph_sin, _, _, _ = init_data(context, public_key, N)
    # decrypt
    time_records = zeros(6)
    for i in 1:6
        start = time_ns()
        decrypt(ciph_sin, private_key)
        finish = time_ns()
        time_records[i] = (finish-start)*10^-9
        GC.gc()
    end
    time_decrypt[i] = sum(time_records[2:end])/5 # first run to precompile
    ptxt = decrypt(ciph_sin, private_key)
    # encrypt
    for i in 1:6
        start = time_ns()
        encrypt(collect(ptxt), public_key, ptxt.context)
        finish = time_ns()
        time_records[i] = (finish-start)*10^-9
        GC.gc()
    end
    time_encrypt[i] = sum(time_records[2:end])/5 # first run to precompile
    release_context_memory()
    GC.gc()
end

table = DataFrame(N_enc_dec = range(1, mult_depth),
                  error = error_encrypt_decrypt,
                  mult_depth_range = vcat(mult_depth_range, zeros(mult_depth-length(mult_depth_range))),
                  time_encrypt = vcat(time_encrypt, zeros(mult_depth-length(mult_depth_range))),
                  time_decrypt = vcat(time_decrypt, zeros(mult_depth-length(mult_depth_range))),
                  N_time = ones(mult_depth)*length(mult_depth_range))

CSV.write("out/efficiency/encryption-decryption.csv", table)
