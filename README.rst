cfme
----

This repository started as a simple command line utility to perform Fourier
Motzkin elimination (FME) for my master thesis. I later added a few handy
related components but then finally moved to python when I had to use more
advanced numerics to implement CHM - a geometric algorithm for performing
variable elimination (aka polyhedral projection). In this sense, the successor
of this project is pystif_.  However, the tools in this repository may still
be useful and are much faster than the pystif equivalents.

.. _pystif: https://github.com/coldfix/pystif


Building
~~~~~~~~

Dependencies:

- boost_
- glpk_

.. _boost: http://www.boost.org/
.. _glpk: https://www.gnu.org/software/glpk/

In a standard linux environment with modern g++ installed the project can be
built by typing::

    make

The binaries are built in the ``bin/`` subfolder.


File format
~~~~~~~~~~~

When reading or writing systems of linear inequalities (SLI) the input/output
has the following simple text format::

    0 -1  0  1
    0  0 -1  1
    0  1  1 -1

Each row holds the coefficients of one inequality. The above example shows the
Shannon cone of two random variables.


Usage
~~~~~

The most general purpose binaries are

- ``check_equivalence`` check if two SLI are equivalent
- ``diff_systems`` print the differing parts of two SLI
- ``elemental-inequalities`` print the Shannon cone for given number of variables
- ``eliminate`` eliminate all but the first few columns of a SLI
- ``minimize_system`` remove all redundant constraints from a SLI

There are a few other binaries which should not be expected to be useful or
even finished. I myself have already forgotten most of their purposes by now.

These two I remember and use for my master thesis:

- ``init-cca`` initialize a cyclic CCA of the specified width using i.i.d.
  initial states
- ``next-layer`` do the same thing using arbitrary specified constraints of
  the initial states. Using this in succession with the ``eliminate`` utility
  the time evolution of a CCA can be computed.
