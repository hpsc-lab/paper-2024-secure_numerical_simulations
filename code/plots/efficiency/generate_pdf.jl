using Plots, LaTeXStrings, Measures, CSV, DataFrames
include("../common_part.jl")

# superscript
function sup(i::Int)
    if i < 0
        c = [Char(0x207B)]
    else
        c = []
    end
    for d in reverse(digits(abs(i)))
        if d == 0 push!(c, Char(0x2070)) end
        if d == 1 push!(c, Char(0x00B9)) end
        if d == 2 push!(c, Char(0x00B2)) end
        if d == 3 push!(c, Char(0x00B3)) end
        if d > 3 push!(c, Char(0x2070+d)) end
    end
    return join(c)
end

# Global options
my_size = (540,360)
my_palette = :seaborn_colorblind
my_linewidth = 3
my_tickfont = font(10)
my_guidefont = font(12)
my_legendfont = font(12)
my_gridlinewidth = 2
my_linestyle = [:solid :dot :dashdot :dashdotdot :dash]
# my_grid = (:xy, :gray, :solid, 0.5, 0.5)

# addition
table = CSV.read("out/efficiency/addition.csv", DataFrame)
# plot error
ticks = [2.0e-13, 4.0e-13, 6.0e-13, 8.0e-13]
ticks_string = ["2"*"×10"*sup(-13),"4"*"×10"*sup(-13),"6"*"×10"*sup(-13),"8"*"×10"*sup(-13)]
p1 = plot(table.N_adds, [table.error_ciph_by_ciph_correlated table.error_ciph_by_ciph_uncorrelated table.error_ciph_by_plain table.error_ciph_by_scalar],
          label=["Ciphertext + ciphertext (correlated)" "Ciphertext + ciphertext (uncorrelated)" "Ciphertext + plaintext" "Ciphertext + scalar"],
          legend_position = :topleft,
          xaxis="Number of addition operations",
          yaxis=L"L^\infty"*" error",
          ylims=(0.0, 8.5e-13),
          yticks=ticks,
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
yticks!(ticks, ticks_string)
savefig(p1, "out/efficiency/add_error.pdf")
# plot time
N_time = Int(table.N_time[1])
p2 = plot(table.mult_depth_range[1:N_time], [table.time_ciph_by_ciph[1:N_time] table.time_ciph_by_plain[1:N_time] table.time_ciph_by_scalar[1:N_time]],
          label=["Ciphertext + ciphertext" "Ciphertext + plaintext" "Ciphertext + scalar"],
          legend_position = :topleft,
          xaxis="Multiplicative depth " * L"$l_\mathrm{max}$",
          yaxis="Runtime [s]",
          ylims=(0.0, 0.042),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p2, "out/efficiency/add_time.pdf")

# bootstrapping
table = CSV.read("out/efficiency/bootstrapping.csv", DataFrame)
# double bootstrapping
table_double = CSV.read("out/efficiency/double_bootstrapping.csv", DataFrame)
# error
ticks = [1.0e-10,1.0e-9,1.0e-8,1.0e-7,1.0e-6,1.0e-5]
ticks_Latex = ["10"*sup(-10),"10"*sup(-9),"10"*sup(-8),"10"*sup(-7),"10"*sup(-6),"10"*sup(-5)]
ticks_bits = [2^(-32),2^(-28),2^(-24),2^(-20),2^(-16)]
p1 = plot(table.N_boots, [table.error table_double.error], label = ["Standard bootstrapping" "Iterative bootstrapping"],
          xaxis="Number of bootstrapping operations",
          yaxis=:log2,
          ylims=(2^(-35), 2^(-14)),
          yticks = ticks,
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
ylabel!(L"L^\infty"*" error")
yticks!(ticks, ticks_Latex)
plot!(twinx(), table.N_boots, fill(NaN, length(table.error)), label="",
      ylims=(2^(-35), 2^(-14)),
      yaxis=(:log2, "Precision [bits]"),
      yticks=(ticks_bits, string.(Int.(round.(-log2.(ticks_bits))))),
      size = my_size,
      linewidth = my_linewidth,
      tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
      gridlinewidth = my_gridlinewidth,
      gridcolor = :white,
      palette = my_palette,
      linestyle = my_linestyle,)
ygrid!(true)
savefig(p1, "out/efficiency/bootstrapping_error.pdf")
# time
N_time = Int(table.N_time[1])
p2 = plot(table.mult_depth_range[1:N_time], [table.time[1:N_time], table_double.time[1:N_time]], label=["Standard bootstrapping" "Iterative bootstrapping"],
          xaxis="Multiplicative depth " * L"$l_\mathrm{max}$",
          yaxis="Runtime [s]",
          ylims=(0.0, 180.0),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p2, "out/efficiency/bootstrapping_time.pdf")

# encode-decode
table = CSV.read("out/efficiency/encode-decode.csv", DataFrame)
# time
p2 = plot(table.mult_depth_range, table.time_encode,
          label="Encode",
          xaxis="Multiplicative depth " * L"$l_\mathrm{max}$",
          yaxis="Runtime [s]",
          ylims=(0.0, 0.13),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p2, "out/efficiency/encode-decode_time.pdf")

# encryption-decryption
table = CSV.read("out/efficiency/encryption-decryption.csv", DataFrame)
# error
p1 = plot(table.N_enc_dec, table.error, label="",
          xaxis="Number of encryption/decryption operations",
          yaxis=L"L^\infty"*" error",
          ylims=(0.0, 2.15e-13),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p1, "out/efficiency/encryption-decryption_error.pdf")
# time
N_time = Int(table.N_time[1])
p2 = plot(table.mult_depth_range[1:N_time], [table.time_encrypt[1:N_time], table.time_decrypt[1:N_time]],
          label=["Encrypt" "Decrypt"],
          xaxis="Multiplicative depth " * L"$l_\mathrm{max}$",
          yaxis="Runtime [s]",
          ylims=(0.0, 1.25),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p2, "out/efficiency/encryption-decryption_time.pdf")

# multiplication
table = CSV.read("out/efficiency/multiplication.csv", DataFrame)
# error
ticks = [5.0e-13, 1.0e-12, 1.5e-12, 2.0e-12]
ticks_string = ["5.0"*"×10"*sup(-13),"1.0"*"×10"*sup(-12),"1.5"*"×10"*sup(-12),"2.0"*"×10"*sup(-12)]
p1 = plot(table.N_mults, [table.error_ciph_by_ciph_correlated table.error_ciph_by_ciph_uncorrelated table.error_ciph_by_plain table.error_ciph_by_scalar],
          label=["Ciphertext " * L"*" * " ciphertext (correlated)" "Ciphertext " * L"*" * " ciphertext (uncorrelated)" "Ciphertext " * L"*" * " plaintext" "Ciphertext " * L"*" * " scalar"],
          legend_position = :topleft,
          xaxis="Number of multiplication operations",
          yaxis=L"L^\infty"*" error",
          ylims=(0.0, 2.4e-12),
          yticks=ticks,
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
yticks!(ticks, ticks_string)
savefig(p1, "out/efficiency/mult_error.pdf")
# plot time
N_time = Int(table.N_time[1])
p2 = plot(table.mult_depth_range[1:N_time], [table.time_ciph_by_ciph[1:N_time] table.time_ciph_by_plain[1:N_time] table.time_ciph_by_scalar[1:N_time]],
          label=["Ciphertext " * L"*" * " ciphertext" "Ciphertext " * L"*" * " plaintext" "Ciphertext " * L"*" * " scalar"],
          legend_position = :topleft,
          xaxis="Multiplicative depth " * L"$l_\mathrm{max}$", 
          yaxis="Runtime [s]",
          ylims=(0.0, 2.55),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p2, "out/efficiency/mult_time.pdf")
# log plot time
N_time = Int(table.N_time[1])
ticks = [0.01, 0.1, 1.0, 10.0]
p2 = plot(table.mult_depth_range[1:N_time], [table.time_ciph_by_ciph[1:N_time] table.time_ciph_by_plain[1:N_time] table.time_ciph_by_scalar[1:N_time]],
          label=["Ciphertext " * L"*" * " ciphertext" "Ciphertext " * L"*" * " plaintext" "Ciphertext " * L"*" * " scalar"],
          legend_position = :topleft,
          xaxis="Multiplicative depth " * L"$l_\mathrm{max}$", 
          yaxis=:log,
          yticks=ticks,
          ylims=(0.008, 17),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
yticks!(ticks, string.(ticks))
ylabel!("Runtime [s]")
savefig(p2, "out/efficiency/mult_time_log.pdf")

# rotate
table = CSV.read("out/efficiency/rotate.csv", DataFrame)
# plot error
p1 = plot(table.N_rotates, [table.error_rotate_1, table.error_rotate_5, table.error_rotate_25],
          label=[L"i = -1" L"i = 5" L"i = -25"],
          legend=:topleft,
          xaxis="Number of rotation operations",
          yaxis=L"L^\infty"*" error",
          ylims=(0.0, 6.5e-14),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p1, "out/efficiency/rotate_error.pdf")
# plot time
N_time = Int(table.N_time[1])
p2 = plot(table.mult_depth_range[1:N_time], [table.time_rotate_1[1:N_time], table.time_rotate_5[1:N_time], table.time_rotate_25[1:N_time]],
          label=[L"i = -1" L"i = 5" L"i = -25"],
          xaxis="Multiplicative depth " * L"$l_\mathrm{max}$", 
          yaxis="Runtime [s]",
          ylims=(0.0, 2.05),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p2, "out/efficiency/rotate_time.pdf")
