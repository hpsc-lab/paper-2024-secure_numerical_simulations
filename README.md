# Secure numerical simulations using fully homomorphic encryption 

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14003844.svg)](https://doi.org/10.5281/zenodo.14003844)
[![arXiv:2410.21824](https://img.shields.io/badge/arXiv-2410.21824-yellow)](https://arxiv.org/abs/2410.21824)

This repository contains information and code to reproduce the results presented in the
article
```bibtex
@online{kholod2024secure,
  title={{S}ecure numerical simulations using fully homomorphic encryption},
  author={Kholod, Arseniy and Polyakov, Yuriy and Schlottke-Lakemper, Michael},
  year={2024},
  month={10},
  doi={10.48550/arXiv.2410.21824},
  eprint={2410.21824},
  eprinttype={arxiv},
  eprintclass={math.NA}
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{kholod2024secureRepro,
  title={Reproducibility repository for
         "{S}ecure numerical simulations using fully homomorphic encryption"},
  author={Kholod, Arseniy and Polyakov, Yuriy and Schlottke-Lakemper, Michael},
  year= {2024},
  howpublished={\url{https://github.com/hpsc-lab/paper-2024-secure_numerical_simulations}},
  doi={10.5281/zenodo.14003844}
}
```


## Abstract

Data privacy is a significant concern when using numerical simulations for sensitive
information such as medical, financial, or engineering data---especially in untrusted
environments like public cloud infrastructures. Fully homomorphic encryption (FHE) offers a
promising solution for achieving data privacy by enabling secure computations directly on
encrypted data. Aimed at computational scientists, this work explores the viability of
FHE-based, privacy-preserving numerical simulations of partial differential equations. The
presented approach utilizes the CKKS scheme, a widely used FHE method for approximate
arithmetic on real numbers. Two Julia packages are introduced, OpenFHE.jl and
SecureArithmetic.jl, which wrap the OpenFHE C++ library to provide a convenient interface
for secure arithmetic operations. With these tools, the accuracy and performance of key FHE
operations in OpenFHE are evaluated, and implementations of finite difference
schemes for solving the linear advection equation with encrypted data are demonstrated. The
results show that cryptographically secure numerical simulations are possible, but that
careful consideration must be given to the computational overhead and the numerical errors
introduced by using FHE. An analysis of the algorithmic restrictions imposed by FHE
highlights potential challenges and solutions for extending the approach to other models and
methods. While it remains uncertain how broadly the approach can be generalized to more complex
algorithms due to CKKS limitations, these findings lay the groundwork for further research on
privacy-preserving scientific computing.


## Numerical experiments

The numerical experiments presented in the paper use
[SecureArithmetic.jl](https://github.com/hpsc-lab/SecureArithmetic.jl).
To reproduce the numerical experiments using SecureArithmetic.jl, you need to install
[Julia](https://julialang.org/).

The subfolder `code` of this repository contains a `README.md` file with
instructions to reproduce the numerical experiments.
Both subfolders also include the result data and scripts for postprocessing.

All numerical experiments were carried out using Julia v1.10.4.


## Authors

- Arseniy Kholod (University of Augsburg, Germany)
- [Yuriy Polyakov](https://ypolyakov.gitlab.io/) (Duality Technologies, Inc., United States)
- [Michael Schlottke-Lakemper](https://lakemper.eu) (University of Augsburg, Germany)


## License
The contents of this repository are available under the [MIT license](LICENSE.md). If you reuse our
code or data, please also cite us (see above).


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
