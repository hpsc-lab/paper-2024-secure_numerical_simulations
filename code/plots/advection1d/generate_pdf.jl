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
table_error = CSV.read("out/advection1d/error.csv", DataFrame)
# plot
p1 = plot(table_error.steps, [table_error.L_inf_cipher_plain_FD table_error.L_inf_cipher_plain_LW],
          label=["Upwind" "Lax-Wendroff"],
          xaxis="Time step",
          yaxis=L"L^\infty"*" error",
          ylims=(0.0, 3.06e-6),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p1, "out/advection1d/ciph_vs_plain_error.pdf")

# time
# read data
table_time = CSV.read("out/advection1d/time.csv", DataFrame)
# plot
p1 = plot(table_time.steps,
          [table_time.timer_ciph_comp_vec_FD table_time.timer_ciph_comp_vec_LW],
          label=["Upwind" "Lax-Wendroff"],
          xaxis="Time step",
          yaxis="Runtime [s]",
          ylims=(0.0, 810.0),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p1, "out/advection1d/execution_times_cipher.pdf")

p2 = plot(table_time.steps, 
          [table_time.timer_plain_comp_vec_FD table_time.timer_plain_comp_vec_LW],
          label=["Upwind" "Lax-Wendroff"],
          legend=:topleft,
          xaxis="Time step",
          yaxis="Runtime [s]",
          yformatter = :scientific,
          ylims=(0.0, 4.1e-4),
          size = my_size,
          linewidth = my_linewidth,
          tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
          gridlinewidth = my_gridlinewidth,
          palette = my_palette,
          linestyle = my_linestyle,
         )
savefig(p2, "out/advection1d/execution_times_plain.pdf")

# visualization
# read data
table_visualization = CSV.read("out/advection1d/visualization.csv", DataFrame)
# plot
p1 = plot(table_visualization.steps, [table_visualization.sol_ciph_FD table_visualization.sol_ciph_LW table_visualization.sol_exact], 
    label=["Upwind" "Lax-Wendroff" "Exact"],
    xaxis="x",
    yaxis="u",
    ylims=(-1.05, 1.05),
    size = my_size,
    linewidth = my_linewidth,
    tickfont = my_tickfont, guidefont = my_guidefont, legendfont = my_legendfont,
    gridlinewidth = my_gridlinewidth,
    palette = my_palette,
    linestyle = my_linestyle,
   )
savefig(p1, "out/advection1d/visualization.pdf")
