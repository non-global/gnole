[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5637033.svg)](https://doi.org/10.5281/zenodo.5637033)

Gnole (Generating functionals for NOn-globaL Evolution)
=================================================

This repository contains the code used in [arXiv:2111.02413](https://arxiv.org/abs/2111.02413 "MC non-global paper").

## About

The [complete user guide](USER_GUIDE.md) documents all options in this
edition, including matching with the SL observable and observable correction.

Gnole is a Monte Carlo program to resum non-global logarithms to next-to-leading logarithmic accuracy.

## Installation

create a build directory and change to it
```
mkdir build; cd build
```

Compile with
```
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
```

executables will be placed in the bin/ directory.
Change back to the source directory (`cd ..`) before using the commands below.

## Running

This edition supports LL and expanded NLL. The obsolete unexpanded NLL implementation
and its inverse-grid/neural interface, present in v1.0, have been removed. `-order 1` selects
expanded NLL directly; the old `-expand-nll` option is no longer accepted. Omitting `-order` 
retains the existing LL default.

```
./bin/gnole -order 0 -nev 1e5 -seed 2 -out ll.dat
./bin/gnole -order 1 -nev 1e5 -seed 2 -out nll.dat
```

Renormalisation and resummation scale factors remain `-xmur` and `-xQ`.
With `-order 1`, the component selectors `-hard-only`,
`-single-insertion-only` and `-double-insertion-only` singles out a given NLL 
correction. Matching to NLO is performed with `-matching` and requires order one 
and the hard contribution (full NLL or `-hard-only`). It can also be combined with
`-sl-obs`.

The option `-obs-correction-only` runs the LL shower and extracts the leading 
measurement shift between physical transverse energy and the SL scale. It implies
`-sl-obs`, requires `-order 0` and currently `-p -1`, and rejects matching,
NLL component selectors and adaptive insertion.

```
./bin/gnole -order 0 -obs-correction-only -mc-errors -nev 1e6 \
  -drap 2 -nbins 40 -lnktmax 10 -seed 1 -out obs_correction.dat
```

## Authors

* Andrea Banfi
* Frédéric A. Dreyer
* Pier Francesco Monni

## References

This code is based on the work published in the following two articles:

* A. Banfi, F. A. Dreyer and P.F. Monni, "Next-to-leading non-global logarithms in QCD,"
  [arXiv:2104.06416](https://arxiv.org/abs/2104.06416 "theory non-global paper")

* A. Banfi, F. A. Dreyer and P.F. Monni, "Higher-order non-global logarithms from jet calculus,"
  [arXiv:2111.02413](https://arxiv.org/abs/2111.02413 "MC non-global paper")

Gnole was built on previous work by Mrinal Dasgupta and Gavin Salam,
and incorporates algorithms and code originally presented in the
following articles:

* M. Dasgupta and G. Salam, "Resummation of nonglobal QCD observables,"
  Phys. Lett. B **512** (2001), 323-330
  [arXiv:hep-ph/0104277 [hep-ph]](https://arxiv.org/abs/hep-ph/0104277 "Original MC for non-global resummation")

* M. Dasgupta and G. Salam, "Accounting for coherence in interjet E(t) flow: A Case study,"
  JHEP **03** (2002), 017
  [arXiv:hep-ph/0203009 [hep-ph]](https://arxiv.org/abs/hep-ph/0203009 "Original resummation for Et")
