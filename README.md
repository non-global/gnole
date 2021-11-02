[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5637033.svg)](https://doi.org/10.5281/zenodo.5637033)

Gnole (Generating functionals for NOn-globaL Evolution)
=================================================

This repository contains the code used in [arXiv:2111.XXXXX](https://arxiv.org/abs/2111.XXXXX "MC non-global paper").

## About

Gnole is a Monte Carlo program to resum non-global logarithms to next-to-leading logarithmic accuracy.

## Installation

create a build directory and change to it
```
mkdir build; cd build
```

For a compilation without neural network options, you can then simply compile with
```
cmake ..
make -j
```

to compile with Torch, use instead
```
cmake -D Torch_DIR="/path/to/torch/share/cmake/Torch" -DNNET=ON ..
make -j
```

executables will be placed in the bin/ directory.

## Authors

* Andrea Banfi
* Frédéric A. Dreyer
* Pier Francesco Monni

## References

This code is based on the work published in the following two articles:

* A. Banfi, F. A. Dreyer and P.F. Monni, "Next-to-leading non-global logarithms in QCD,"
  [arXiv:2104.06416](https://arxiv.org/abs/2104.06416 "theory non-global paper")

* A. Banfi, F. A. Dreyer and P.F. Monni, "Higher-order non-global logarithms from jet calculus,"
  [arXiv:2111.XXXXX](https://arxiv.org/abs/2111.XXXXX "MC non-global paper")

Gnole was built on previous work by Mrinal Dasgupta and Gavin Salam,
and incorporates algorithms and code originally presented in the
following articles:

* M. Dasgupta and G. Salam, "Resummation of nonglobal QCD observables,"
  Phys. Lett. B **512** (2001), 323-330
  [arXiv:hep-ph/0104277 [hep-ph]](https://arxiv.org/abs/hep-ph/0104277 "Original MC for non-global resummation")

* M. Dasgupta and G. Salam, "Accounting for coherence in interjet E(t) flow: A Case study,"
  JHEP **03** (2002), 017
  [arXiv:hep-ph/0203009 [hep-ph]](https://arxiv.org/abs/hep-ph/0203009 "Original resummation for Et")
