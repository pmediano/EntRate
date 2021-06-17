EntRate -- Entropy rate estimators for neuroscience
===================================================

This library contains Octave/Matlab (>R2016a) functions to compute entropy
rate- and Lempel-Ziv complexity-related functions in continuous and discrete
data. An example application of LZ complexity to neuroimaging data with a
relevant discussion on entropy rate can be found here:

* Mediano, P., Rosas, F., Timmermann, C., Roseman, L., Nutt, ...  Bor, D. &
  Carhart-Harris, R. L. (2020). Effects of external stimulation on psychedelic
  state neurodynamics.
  [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.11.01.356071v1).


**Note**: this repository contains a copy of the
[MVGC2](https://github.com/SacklerCentre/MVGC2) library, lightly modified to be
Octave-compatible. Full credit for the MVGC2 library goes to the original
authors, and interested users are referred to the [original
repository](https://github.com/SacklerCentre/MVGC2) and [associated
papers](https://doi.org/10.1016/j.jneumeth.2013.10.018) for details. It also
contains a precompiled (unmodified) jar of the VMM code by Ron Begleiter (see
[paper](https://arxiv.org/abs/1107.0051) for details).


Download and installation
-------------------------

If you intend to use the LZ76 entropy rate estimator, you need to compile it
first by running:

```octave
mex COPTIMFLAGS="-O3" LZ76.c
```

Once this is done, you can run `LZ76()` and call `help LZ76` as usual.

Note, however, that if you use Octave you will need to install the `statistics`
package. You can do this simply by running `pkg install -forge io statistics`.

Tests are provided in the `tests/` subfolder. To run them in Matlab, run
`runtests('tests/')` from this repository's root folder.


Usage
-----

This library implements a few entropy rate estimators (via state-space entropy
[CSER], context tree-weighted predictor [CTW], and Lempel-Ziv compression)
applicable to discrete or continuous data. For example, to compute the entropy
rate of a maximum-entropy random sequence of bits, run:

```octave
X = 1*(rand([1, 1000]) < 0.5);
H = CTWEntropyRate(X);
```

Simple example scripts comparing LZ, CTW and CSER can be found in the
`examples/` folder.

Feature requests and bug reports are warmly welcome. Email Pedro Mediano (see
email in the paper above) for any questions or comments.


Use in Python
-------------

The `CTWEntropyRate` and `StateSpaceEntropyRate` functions in this repository
are Octave-friendly, which means they can be easily called from Python through
the wholesome [oct2py](https://oct2py.readthedocs.io/) package, as long as you
have a functional Octave installation.

You can simply run (from the root folder of the repo):

```python
import numpy.random as rn
from oct2py import Oct2Py

oc = Oct2Py()
oc.StateSpaceEntropyRate(rn.randn(1, 1000))
```

Note that `StateSpaceEntropyRate` requires the `statistics` Octave package
(which you can install running `pkg install -forge io statistics` from Octave),
and `CTWEntropyRate` requires an Octave installation with Java support.


Licence
-------

This software is distributed under the GNU General Public Licence, version 3.


Further reading
---------------

* Barnett, L., & Seth, A. (2015). Granger causality for state-space models.
  Physical Review E, 91(4), 040101.

* Begleiter, R., El-Yaniv, R., & Yona, G. (2004). On prediction using variable
  order Markov models. Journal of Artificial Intelligence Research, 22,
  385-421.

* Kaspar, F., & Schuster, H. G. (1987). Easily calculable measure for the
  complexity of spatiotemporal patterns. Physical Review A, 36(2), 842.

* Mediano, P., Rosas, F., Barrett, A., & Bor, D. (2020). Decomposing
  spectral and phasic differences in non-linear features between datasets.
  arXiv:2009.10015.


(C) Pedro Mediano, Fernando Rosas and Andrea Luppi, 2019-21

