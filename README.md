[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

Gnole (Generating functionals for NOn-globaL Evolution)
============

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
* Frederic Dreyer
* Pier Monni

## References

* A. Banfi, F. A. Dreyer and P.F. Monni, "Next-to-leading non-global logarithms in QCD,"
  [arXiv:2104.06416](https://arxiv.org/abs/2104.06416 "theory non-global paper")

* A. Banfi, F. A. Dreyer and P.F. Monni, "Higher-order non-global logarithms from jet calculus,"
  [arXiv:2111.XXXXX](https://arxiv.org/abs/2111.XXXXX "MC non-global paper")
