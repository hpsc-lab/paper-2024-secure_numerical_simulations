using CSV, DataFrames, OpenFHE
include("../../src/secure_simulations.jl")
include("../common_part.jl")

GC.gc()
GC.gc()
println("Free space: ", Sys.free_memory()/2^20, "MB")

# error investigation
N = 64
mult_depth = 15
private_key, public_key, context, _ = init_setup(mult_depth, N)
ciph_sin, ciph_factor, ptxt_factor, u_sin = init_data(context, public_key, N)

# multiply ciphertext by ciphertext
# uncorrellated
res = deepcopy(ciph_sin)
exact = deepcopy(u_sin)
error_ciph_by_ciph_uncorrelated = zeros(mult_depth)
for i in range(1, mult_depth)
    global res = res*ciph_factor
    # create new ciphertext with the same level
    global ciph_factor = encrypt(PlainVector(ones(N).*(1 + pi/30), context, level(res)+1), public_key)
    global exact = exact*(1 + pi/30)
    error_ciph_by_ciph_uncorrelated[i] = maximum(abs.(exact.-collect(decrypt(res, private_key))))
end
# correllated
res = deepcopy(ciph_sin)
exact = deepcopy(u_sin)
error_ciph_by_ciph_correlated = zeros(mult_depth)
cc = SecureArithmetic.get_crypto_context(context)
ciph_factor = encrypt(ptxt_factor, public_key)
for i in range(1, mult_depth)
    global res = res*ciph_factor
    # need to drop level
    OpenFHE.ModReduceInPlace(cc, ciph_factor.data)
    global exact = exact*(1 + pi/30)
    error_ciph_by_ciph_correlated[i] = maximum(abs.(exact.-collect(decrypt(res, private_key))))
end
# multiply ciphertext by plaintext
res = deepcopy(ciph_sin)
exact = deepcopy(u_sin)
error_ciph_by_plain = zeros(mult_depth)
for i in range(1, mult_depth)
    global res = res*ptxt_factor
    global exact = exact*(1 + pi/30)
    error_ciph_by_plain[i] = maximum(abs.(exact.-collect(decrypt(res, private_key))))
end
# multiply ciphertext by scalar
res = deepcopy(ciph_sin)
exact = deepcopy(u_sin)
error_ciph_by_scalar = zeros(mult_depth)
for i in range(1, mult_depth)
    global res = res*(1 + pi/30)
    global exact = exact*(1 + pi/30)
    error_ciph_by_scalar[i] = maximum(abs.(exact.-collect(decrypt(res, private_key))))
end
release_context_memory()
GC.gc()


# runtime investigations
mult_depth_range = Int.(range(5, 25, 5))
time_ciph_by_ciph = zeros(length(mult_depth_range))
time_ciph_by_plain = zeros(length(mult_depth_range))
time_ciph_by_scalar = zeros(length(mult_depth_range))
for i in range(1, length(mult_depth_range))
    local mult_depth = mult_depth_range[i]
    local private_key, public_key, context, mult_depth_resulted = init_setup(mult_depth, N)
    mult_depth_range[i] = mult_depth_resulted
    local ciph_sin, ciph_factor, ptxt_factor, _ = init_data(context, public_key, N)
    time_records = zeros(6)
    # multiply ciphertext by ciphertext
    for i in 1:6
        start = time_ns()
        ciph_sin = ciph_sin*ciph_factor
        finish = time_ns()
        time_records[i] = (finish-start)*10^-9
        GC.gc()
    end
    time_ciph_by_ciph[i] = sum(time_records[2:end])/5 # first run to precompile
    # multiply ciphertext by plaintext
    for i in 1:6
        start = time_ns()
        ciph_sin = ciph_sin*ptxt_factor
        finish = time_ns()
        time_records[i] = (finish-start)*10^-9
        GC.gc()
    end
    time_ciph_by_plain[i] = sum(time_records[2:end])/5 # first run to precompile
    # multiply ciphertext by scalar
    for i in 1:6
        start = time_ns()
        ciph_sin = ciph_sin*(1+pi/30)
        finish = time_ns()
        time_records[i] = (finish-start)*10^-9
        GC.gc()
    end
    time_ciph_by_scalar[i] = sum(time_records[2:end])/5 # first run to precompile
    release_context_memory()
    GC.gc()
end

table = DataFrame(N_mults = range(1, mult_depth),
                  error_ciph_by_ciph_uncorrelated = error_ciph_by_ciph_uncorrelated,
                  error_ciph_by_ciph_correlated = error_ciph_by_ciph_correlated,
                  error_ciph_by_plain = error_ciph_by_plain,
                  error_ciph_by_scalar = error_ciph_by_scalar,
                  mult_depth_range = vcat(mult_depth_range, zeros(mult_depth-length(mult_depth_range))),
                  time_ciph_by_ciph = vcat(time_ciph_by_ciph, zeros(mult_depth-length(mult_depth_range))),
                  time_ciph_by_plain = vcat(time_ciph_by_plain, zeros(mult_depth-length(mult_depth_range))),
                  time_ciph_by_scalar = vcat(time_ciph_by_scalar, zeros(mult_depth-length(mult_depth_range))),
                  N_time = ones(mult_depth)*length(mult_depth_range))

CSV.write("out/efficiency/multiplication.csv", table)
