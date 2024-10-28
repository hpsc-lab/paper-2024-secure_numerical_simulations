using CSV, DataFrames
include("../../src/secure_simulations.jl")
include("../common_part.jl")

GC.gc()
GC.gc()
println("Free space: ", Sys.free_memory()/2^20, "MB")

N = 64

# runtime investigation
mult_depth_range = Int.(range(5, 25, 5))
time_encode = zeros(length(mult_depth_range))
for i in range(1, length(mult_depth_range))
    local mult_depth = mult_depth_range[i]
    local _, public_key, context, mult_depth_resulted = init_setup(mult_depth, N)
    mult_depth_range[i] = mult_depth_resulted
    local _, _, _, u_sin = init_data(context, public_key, N)
    # encode
    time_records = zeros(6)
    for i in 1:6
        start = time_ns()
        ptxt = PlainVector(u_sin, context)
        OpenFHE.Encode(ptxt.data[])
        finish = time_ns()
        time_records[i] = (finish-start)*10^-9
        GC.gc()
    end
    time_encode[i] = sum(time_records[2:end])/5 # first run to precompile
    release_context_memory()
    GC.gc()
end

table = DataFrame(mult_depth_range = mult_depth_range,
                  time_encode = time_encode)

CSV.write("out/efficiency/encode-decode.csv", table)
