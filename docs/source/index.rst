.. genominterv documentation master file, created by
   sphinx-quickstart on Sat May 18 10:24:07 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to genominterv
============================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

:py:mod:`genomeinterv` provides support for working with intervals on genomes. A
genomic interval is specified as a chromosome, start, and end. It is half-open
so that a value ``x`` is in an interval ``(start, end)`` included in the inrval
if ``start <= x and x < end``. All functions take `pandas.DataFrames
<http://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`_ as
arguments. These data frames must include ``chrom``, ``start``, and ``end``
columns.

See the `library reference <code.html>`_ for detailed documentation of each
function and decorator.

Set operations
=================

The three functions:

- :any:`interval_diff`
- :any:`interval_intersect`
- :any:`interval_union`

do the standard difference, intersection and union set operations on two sets of
genomic intervals. The intervals returned from all three functions are collapsed
to produce non-overlapping intervals. The genomic intervals in each set must be
non-overlapping. This can be achieved using function:

- :any:`interval_collapse`

which produces the union of genomic intervals in a single set genomic of
intervals.

Genomic decorator
======================

To make it easy to create other interval functions that work across chromosomes,
the module provides a :any:`genomic` decorator that can be applied to functions
that operate lists of ``(start, end)`` tuples. Applying the decorator changes
the signature of a function to make it operate on DataFrames that include
``chrom``, ``start``, and ``end`` columns. Here is an example function that
shifts intervals by 1000bp::

    @genomic
    def inverval_shift(tuples):
        return [(x+1000, y+1000) for (x, y) in tuples]

    df = pandas.DataFrame()

    shifted = inverval_shift(df)

Remapping functions
======================

The function :any:`interval_distance` onverts coordinates of one set of
genomic intervals into distances to the closest interval in a second set.
:any:`interval_relative_distance` does the same but returns
relative distances.

Two-set statistics
=====================

The module also provides two statistics for relations between sets:
:any:`jaccard` computes the `Jaccard index <https://en.wikipedia.org/wiki/Jaccard_index>`_
statistic for two sets of genomic intervals.


Bootstrap decorator
======================

The module provides a :any:`bootstrap` decorator that turns a function producing
a statistic into one that also produces a p-value. The bootstrapping resamples
the intervals of the second argument for each chromosome independently. Only
required argument to bootstrap is the name of the chromosome assembly used.

This example does this for the provided :any:`jaccard` satistic::

    @bootstrap('hg19', samples=1000)
    def jaccard_test(query, annot):
        return jaccard(query, annot)

    jaccard_stat, p_value = jaccard_test(intervals, other_intervals)

The decorator works on any function that takes two sets of intervals.


Ready-made tests
==================

:any:`proximity_test` computes tests if intervals in one set is significantly
proximal to intervals in another
set.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

