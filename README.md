# Secure numerical simulations using fully homomorphic encryption 

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14003844.svg)](https://doi.org/10.5281/zenodo.14003844)

This repository contains information and code to reproduce the results presented in the
article
```bibtex
@online{kholod2024secure,
  title={{S}ecure numerical simulations using fully homomorphic encryption},
  author={Kholod, Arseniy and Polyakov, Yuriy and Schlottke-Lakemper, Michael},
  year={2023},
  month={11},
  doi={10.48550/arXiv.TODO},
  eprint={TODO},
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

Data privacy is a significant concern in many environments today.
This is particularly true if sensitive information, e.g.,  engineering, medical, or
financial data, is to be processed on potentially insecure systems, as it is often the
case in cloud computing. Fully homomorphic encryption (FHE) offers a potential solution
to this problem, as it allows for secure computations on encrypted data. In this paper,
we investigate the viability of using FHE for privacy-preserving numerical simulations of partial
differential equations. We first give an overview of the CKKS scheme, a popular FHE
method for computations with real numbers. This is followed by an introduction of our
Julia packages OpenFHE.jl and SecureArithmetic.jl, which provide a Julia wrapper for the
C++ library OpenFHE and offer a user-friendly interface for secure arithmetic
operations. We then present a performance analysis of the CKKS scheme within OpenFHE,
focusing on the error and efficiency of different FHE operations. Finally, we
demonstrate the application of FHE to secure numerical simulations by implementing two
finite difference schemes for the linear advection equation using the
SecureArithmetic.jl package. Our results show that FHE can be used to perform
cryptographically secure numerical simulations, but that the error and efficiency of FHE
operations must be carefully considered when designing applications.


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
