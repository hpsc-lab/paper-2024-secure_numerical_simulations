# Numerical experiments
In this folder and its subfolders, all code to reproduce the numerical experiments in the paper is
located.  Furthermore, it also contains the raw data for the results reported in the paper.


## Set up Julia
Install Julia following the instructions at https://julialang.org/downloads/. We have used Julia
v1.10.4 for all our results.

To install all necessary Julia packages, execute the following statement from within the folder that
contains the `README.md` file you are currently reading:
```shell
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```
This will recreate the exact Julia environment we used to obtain our results for full
reproducibility.


## Run the experiments
To re-create the results, run one of the following commands.

*Note:* Running all scripts may take more than 30 hours, depending on your system. Our tests were
conducted on a machine with an AMD Ryzen Threadripper 3990X 64 core processor and 256 GB RAM. Since
many of the scripts use a lot of memory, the execution may be aborted if your machine runs out of
memory.

### Accuracy and performance analysis of the CKKS scheme (Sec. 3)
```bash
OMP_NUM_THREADS=1 julia --project=. ./plots/efficiency/runall.jl
```

### Convergence tests (Sec. 5)
```bash
OMP_NUM_THREADS=8 julia --project=. ./plots/advection1d/convergence_test.jl
OMP_NUM_THREADS=8 julia --project=. ./plots/advection2d/convergence_test.jl
```

### Accuracy and performance analysis of the secure numerical simulations (Sec. 5)
```bash
OMP_NUM_THREADS=1 julia --project=. ./plots/advection1d/configure_plain.jl
OMP_NUM_THREADS=1 julia --project=. ./plots/advection2d/configure_plain.jl
OMP_NUM_THREADS=1 julia --project=. ./plots/advection1d/run_examples.jl
OMP_NUM_THREADS=1 julia --project=. ./plots/advection2d/run_examples.jl
OMP_NUM_THREADS=1 julia --project=. ./plots/advection1d/time_per_step.jl
OMP_NUM_THREADS=1 julia --project=. ./plots/advection2d/time_per_step.jl
OMP_NUM_THREADS=1 julia --project=. ./plots/advection2d/LW_diff_mult_depth.jl
OMP_NUM_THREADS=1 julia --project=. ./plots/advection2d/run_examples_LW_double_bootstrapping.jl
```

### Parallelization results (Sec. 5)
```bash
OMP_NUM_THREADS=1 julia --project=. plots/advection2d/parallelization.jl
OMP_NUM_THREADS=2 julia --project=. plots/advection2d/parallelization.jl
OMP_NUM_THREADS=4 julia --project=. plots/advection2d/parallelization.jl
OMP_NUM_THREADS=8 julia --project=. plots/advection2d/parallelization.jl
OMP_NUM_THREADS=16 julia --project=. plots/advection2d/parallelization.jl
OMP_NUM_THREADS=32 julia --project=. plots/advection2d/parallelization.jl
OMP_NUM_THREADS=64 julia --project=. plots/advection2d/parallelization.jl
```


## Plot generation
To generate the plots after running the experiments, execute the following lines:
```bash
julia --project=. ./plots/efficiency/generate_pdf.jl
julia --project=. ./plots/advection1d/generate_pdf.jl
julia --project=. ./plots/advection2d/generate_pdf.jl
```
All plots will be generated in the corresponding subfolder within the `out/` directory, i.e.,
`out/efficiency`, `out/advection1d`, or `out/advection2d`.


## Results
The results we obtained from running the experiments on our machine can be found in the folder
[`out/`](out/).
