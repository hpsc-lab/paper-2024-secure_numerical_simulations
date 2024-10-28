using Plots, LaTeXStrings, Measures, CSV, DataFrames
include("../common_part.jl")

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

# error
# read data 
table_error = CSV.read("out/advection2d/error.csv", DataFrame)
# plot
p1 = plot(table_error.steps, [table_error.L_inf_cipher_plain_FD table_error.L_inf_cipher_plain_LW],
          label=["Upwind" "Lax-Wendroff"],
          xaxis="Time step",
          yaxis=L"L^\infty"*" error",
          ylims=(0.0, 8.1e-6),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p1, "out/advection2d/ciph_vs_plain_error.pdf")

# double bootstrapping
table_error = CSV.read("out/advection2d/error_LW_double_bootstrapping.csv", DataFrame)
# plot
p1 = plot(table_error.steps, table_error.L_inf_cipher_plain_LW,
          label=nothing,
          xaxis="Time step",
          yaxis=L"L^\infty"*" error",
          ylims=(0.0, 8.1e-9),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p1, "out/advection2d/ciph_vs_plain_error_double_bootstrapping.pdf")

# time
# read data
table_time = CSV.read("out/advection2d/time.csv", DataFrame)
# plot
p1 = plot(table_time.steps,
          [table_time.timer_ciph_comp_vec_FD table_time.timer_ciph_comp_vec_LW],
          label=["Upwind" "Lax-Wendroff"],
          xaxis="Time step",
          yaxis="Runtime [s]",
          ylims=(0.0, 6100.0),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p1, "out/advection2d/execution_time_cipher.pdf")
       
p2 = plot(table_time.steps, 
          [table_time.timer_plain_comp_vec_FD table_time.timer_plain_comp_vec_LW],
          label=["Upwind" "Lax-Wendroff"],
          xaxis="Time step",
          yaxis="Runtime [s]",
          yformatter = :scientific,
          ylims=(0.0, 6.1e-2),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p2, "out/advection2d/execution_time_plain.pdf")

# double bootstrapping
table_time = CSV.read("out/advection2d/time_LW_double_bootstrapping.csv", DataFrame)
# plot
p1 = plot(table_time.steps,
          table_time.timer_ciph_comp_vec_LW,
          label=nothing,
          xaxis="Time step",
          yaxis="Runtime [s]",
          ylims=(0.0, 8100.0),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p1, "out/advection2d/execution_time_cipher_double_bootstrapping.pdf")

# visualization
# read data
table_visualization = CSV.read("out/advection2d/visualization.csv", DataFrame)
# plot
# errors
Linf_FD = round(L_inf(table_visualization.sol_exact, table_visualization.sol_ciph_FD); sigdigits=2)
Linf_LW = round(L_inf(table_visualization.sol_exact, table_visualization.sol_ciph_LW); sigdigits=2)
N = Int(table_visualization.length[1])
x = range(0.0, 1.0, N)
y = x

p1 = plot(x, y, reshape(table_visualization.sol_ciph_FD, (N,N)); st=:surface, colorbar=false,
          xaxis="x", yaxis="y", zaxis="u",
          xlims=(0, 1.0), ylims=(0, 1.0), zlims=(-1.0, 1.0),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p1, "out/advection2d/visualization_cipher_FD.pdf")
p2 = plot(x, y, reshape(table_visualization.sol_ciph_LW, (N,N)); st=:surface, colorbar=false,
          xaxis="x", yaxis="y", zaxis="u",
          xlims=(0, 1.0), ylims=(0, 1.0), zlims=(-1.0, 1.0),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p2, "out/advection2d/visualization_cipher_LW.pdf")
p3 = plot(x, y, reshape(table_visualization.sol_exact, (N,N)); st=:surface, colorbar=false,
          xaxis="x", yaxis="y", zaxis="u",
          xlims=(0, 1.0), ylims=(0, 1.0), zlims=(-1.0, 1.0),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p3, "out/advection2d/visualization_exact.pdf")

# parallelization
# read data
table_parallelization = CSV.read("out/advection2d/parallelization.csv", DataFrame)
# plot
p1 = plot(table_parallelization.OMP_THREADS,
          [table_parallelization.comp_time table_parallelization.conf_time],
          label=["Compute" "Configure"],
          xaxis=:log2,
          yaxis="Runtime [s]",
          ylims=(0.0, 6100.0),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
xlabel!("Number of OpenMP threads")
xticks!([1,2,4,8,16,32,64], string.([1,2,4,8,16,32,64]))
savefig(p1, "out/advection2d/parallelization_time.pdf")

speedup_conf = table_parallelization.conf_time[1] ./ table_parallelization.conf_time
speedup_comp = table_parallelization.comp_time[1] ./ table_parallelization.comp_time
p2 = plot(table_parallelization.OMP_THREADS, 
          [speedup_comp speedup_conf],
          xaxis=:log2,
          label=["Compute" "Configure"],
          yaxis="Speedup",
          ylims=(0.0, 12.7),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
xlabel!("Number of OpenMP threads")
xticks!([1,2,4,8,16,32,64], string.([1,2,4,8,16,32,64]))
savefig(p2, "out/advection2d/parallelization_speed_up.pdf")

efficiency_conf = 100*speedup_conf./table_parallelization.OMP_THREADS 
efficiency_comp = 100*speedup_comp./table_parallelization.OMP_THREADS
p3 = plot(table_parallelization.OMP_THREADS,
          [efficiency_comp efficiency_conf],
          xaxis=:log2,
          label=["Compute" "Configure"],
          yaxis="Efficiency [%]",
          ylims=(0.0, 102.0),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
xlabel!("Number of OpenMP threads")
xticks!([1,2,4,8,16,32,64], string.([1,2,4,8,16,32,64]))
savefig(p3, "out/advection2d/parallelization_efficiency.pdf")

# time_per_step and error per step
table= CSV.read("out/advection2d/LW_diff_mult_depth.csv", DataFrame)
# plot
plot(table.mult_depth_range,
     table.LW_error_vec,
     xaxis=L"$l_\mathrm{refresh}$",
     yaxis=L"L^\infty"*" error",
     ylims=(0, 4.1e-5),
     size = my_size,
     linewidth = my_linewidth,
     tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
     gridlinewidth = my_gridlinewidth,
     palette = my_palette,
     linestyle = my_linestyle,

    )

p1 = plot!(twinx(),
           table.mult_depth_range,
           [fill(NaN, length(table.LW_error_vec)) table.LW_comp_time_vec],
           label=["Error" "Runtime"],
           yaxis="Runtime [s]",
           ylims=(0, 12300),
           yticks=[0, 3000, 6000, 9000, 12000],
           grid=true,
           yformatter=:plain,
           size = my_size,
           linewidth = my_linewidth,
           tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
           gridlinewidth = my_gridlinewidth,
           palette = my_palette,
           linestyle = my_linestyle,
          )

savefig(p1, "out/advection2d/time_error_per_step.pdf")