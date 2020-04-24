import pandas
import numpy

from functools import wraps, reduce, partial
import bisect, math, random, sys
from itertools import chain
from statsmodels.distributions.empirical_distribution import ECDF

from collections import namedtuple

from genominterv.chrom_sizes import chrom_sizes

def by_chrom(func):
    """
    Decorator that converting a function that operates on a data frame with only
    one chromosome to one operating on a data frame with many chromosomes.
    """
    @wraps(func)
    def wrapper(*args):

        # make a local copy with reset indexes
        data_frames = [df.reset_index() for df in args]

        # get all chromosoems in arguments
        chromosomes = set()
        for df in data_frames:
            chromosomes.update(df['chrom'].unique())
        chromosomes = sorted(chromosomes)

        # get indexes (possibly none) for each chromosome in each frame
        idx = list()
        for df in data_frames:
            d = dict((chrom, []) for chrom in chromosomes)
            gr = df.groupby('chrom').groups
            d.update(gr)
            idx.append(d)

        # call func on subsets of each argument matching a chromosome
        results = list()
        for chrom in chromosomes:
            func_args = list()
            for i, df in enumerate(data_frames):
                func_args.append(df.loc[idx[i][chrom]])
            results.append(func(*func_args))

        return pandas.concat(results).reset_index()
    return wrapper


def with_chrom(func):
    """
    Decorator for converting a function operating on (start, end) tuples to one
    that takes data frames with chrom, start, end columns. Also sorts intervals.
    """
    @wraps(func)
    def wrapper(*args):
        chrom_set = set()
        tps_list = list()
        for df in args:            
            chrom_set.update(df['chrom'])
            tps = sorted(zip(df['start'], df['end']))
            tps_list.append(tps)
        assert len(chrom_set) == 1
        chrom = chrom_set.pop()
        res_df = pandas.DataFrame.from_records(func(*tps_list), columns = ['start', 'end'])
        res_df['chrom'] = chrom
        return res_df
    return wrapper


def genomic(func):
    """
    Decorator for converting a function operating on (start, end) tuples to one
    that takes data frames with chrom, start, end columns and executes on each
    chromosome individually.
    """
    @wraps(func)
    @by_chrom
    @with_chrom
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper

# In all of the following, the list of intervals must be sorted and 
# non-overlapping. We also assume that the intervals are half-open, so
# that x is in tp(start, end) iff start <= x and x < end.

def flatten(list_of_tps):
    """
    Convert a list of sorted intervals to a list of endpoints.

    :param query: Sorted list of (start, end) tuples.
    :type query: list
    :returns: A list of interval ends
    :rtype: list
    """
    return reduce(lambda ls, ival: ls + list(ival), list_of_tps, [])


def unflatten(list_of_endpoints):
    """
    Convert a list of sorted endpoints into a list of intervals.

    :param query: Sorted list of ends.
    :type query: list
    :returns: A list of intervals.
    :rtype: list
    """
    return [ [list_of_endpoints[i], list_of_endpoints[i + 1]]
          for i in range(0, len(list_of_endpoints) - 1, 2)]


def merge(query, annot, op):
    """
    Merge two lists of sorted intervals according to the boolean function op.

    :param query: List of (start, end) tuples.
    :type query: list
    :param query: List of (start, end) tuples.
    :type query: list
    :param op: Boolean function taking two a bolean arguments.
    :type op: function
    :returns: A list of interval merged according to op
    :rtype: list
    """
    a_endpoints = flatten(query)
    b_endpoints = flatten(annot)

    assert a_endpoints == sorted(a_endpoints), "not sorted or non-overlaping"
    assert b_endpoints == sorted(b_endpoints), "not sorted or non-overlaping"


    sentinel = max(a_endpoints[-1], b_endpoints[-1]) + 1
    a_endpoints += [sentinel]
    b_endpoints += [sentinel]

    a_index = 0
    b_index = 0

    res = []

    scan = min(a_endpoints[0], b_endpoints[0])
    while scan < sentinel:
        in_a = not ((scan < a_endpoints[a_index]) ^ (a_index % 2))
        in_b = not ((scan < b_endpoints[b_index]) ^ (b_index % 2))
        in_res = op(in_a, in_b)

        if in_res ^ (len(res) % 2):
            res += [scan]
        if scan == a_endpoints[a_index]: 
            a_index += 1
        if scan == b_endpoints[b_index]: 
            b_index += 1
        scan = min(a_endpoints[a_index], b_endpoints[b_index])

    return unflatten(res)

def diff(a, b):
    if not (a and b):
        return a and a or b
    return merge(a, b, lambda in_a, in_b: in_a and not in_b)

def union(a, b):
    if not (a and b):
        return []
    return merge(a, b, lambda in_a, in_b: in_a or in_b)

def intersect(a, b):
    if not (a and b):
        return []
    return merge(a, b, lambda in_a, in_b: in_a and in_b)
    
@genomic
def interval_diff(query, annot):
    """
    This function computes the difference between two sets of genomic intervals.
    The genomic intervals in each set must be non-overlapping. This can be
    achieved using :any:interval_collapse.

    :param query: data frame with query intervals.
    :type query: pandas.DataFrame
    :param annot: data frame with annotation intervals.
    :type annot: pandas.DataFrame
    :returns: A data frame with chr, start, end columns representing the difference.
    :rtype: pandas.DataFrame
    """ 
    return diff(query, annot)


@genomic
def interval_union(query, annot):
    """
    This function computes the union of two sets of genomic intervals. The genomic intervals
    in each set must be non-overlapping. This can be achieved using :any:interval_collapse.

    :param query: data frame with query intervals.
    :type query: pandas.DataFrame
    :param annot: data frame with annotation intervals.
    :type annot: pandas.DataFrame
    :returns: A data frame with chr, start, end columns representing the union.
    :rtype: pandas.DataFrame
    """ 
    return union(query, annot)


@genomic
def interval_intersect(query, annot):
    """
    # This function computes the intersection of two sets of genomic intervals. The genomic intervals
    in each set must be non-overlapping. This can be achieved using :any:interval_collapse.

    :param query: data frame with query intervals.
    :type query: pandas.DataFrame
    :param annot: data frame with annotation intervals.
    :type annot: pandas.DataFrame
    :returns: A data frame with chr, start, end columns representing the intersection.
    :rtype: pandas.DataFrame
    """ 
    return intersect(query, annot)


@genomic
def interval_collapse(interv):
    """
    This function computes the union of intervals in a single set.

    :param interv: data frame with intervals.
    :type interv: pandas.DataFrame
    :returns: A data frame with chr, start, end columns representing the union.
    :rtype: pandas.DataFrame
    """
    interv_union = [list(interv[0])]
    for i in range(1, len(interv)):
        x = interv[i]
        if interv_union[-1][1] < x[0]:
            interv_union.append(list(x))
        else:
            interv_union[-1][1] = x[1]
    return interv_union


def remap(query, annot, relative=False, include_prox_coord=False):
    """
    Remap the coordinates of a single interval in ``query`` to the distance from
    the closet interval in ``annot``. Returns empty set if annot is empty for
    the chromosome. Intervals from ``query`` that overlap intervals in ``annot``
    are discarded.

    :param query: A ``(start, end)`` tuple.
    :type query: tuple or list
    :param annot: List of ``(start, end)`` tuples on the same chromosome.
    :type annot: list
    :returns: List of ``(start, end)`` tuples
    :rtype: list
    """

    if not (query and annot):
        return []

    query_start, query_end = query
    annot_starts, annot_ends = zip(*annot)

    # find interval betweent two annotations
    idx = bisect.bisect_right(annot_ends, query_start)

    interval_start = idx != 0 and annot_ends[idx-1] or None
    interval_end = idx != len(annot) and annot_starts[idx] or None
    assert interval_start is not None or interval_end is not None, "maybe your query annd annot overlaps..."

    if interval_start is None and interval_end is None:
        # spans accross all annotation
        return []
    if not ((interval_start is None or interval_start < query_start) and (interval_end is None or interval_end >= query_end)):
        # spans at least one annotation interval
        return []

    if interval_start is None:
        interval_mid = - float('inf')
    elif interval_end is None:
        interval_mid = float('inf')
    else:
        # assert query_end <= interval_end, "query and annotations coordinates overlap"
        interval_mid = int(interval_end + math.ceil((interval_start - interval_end) / 2.0))

    if interval_mid < query_start:
        if include_prox_coord:
            remapped = [(query_end - interval_end, query_start - interval_end,
                        idx == 0 and numpy.nan or annot_starts[idx], 
                        idx == 0 and numpy.nan or annot_ends[idx])]
        else:
            remapped = [(query_end - interval_end, query_start - interval_end)]
    elif interval_mid >= query_end:
        if include_prox_coord:
            remapped = [(query_start - interval_start, query_end - interval_start,
                        idx == len(annot_starts) and numpy.nan or annot_starts[idx-1], 
                        idx == len(annot_ends) and numpy.nan or annot_ends[idx-1])]
        else:
            remapped = [(query_start - interval_start, query_end - interval_start)]            
    else:
        if include_prox_coord:
            remapped = [(query_start - interval_start, interval_mid - interval_start,
                        idx == 0 and numpy.nan or annot_starts[idx-1], 
                        idx == 0 and numpy.nan or annot_ends[idx-1]),
                        (query_end - interval_end, interval_mid - interval_end,
                        idx == len(annot_starts) and numpy.nan or annot_starts[idx], 
                        idx == len(annot_ends) and numpy.nan or annot_ends[idx])]
        else:
            remapped = [(query_start - interval_start, interval_mid - interval_start),
                        (query_end - interval_end, interval_mid - interval_end)]

    # FIXME: Test above change and update documentations

    # FIXME: The stuff below should be done before the include_prox_coord stuff is done...

    # compute remapped distance relative to the interval length (so that is is max 0.5)
    if relative:
        if interval_start is None or interval_end is None:
            remapped = [(numpy.nan, numpy.nan)]
        else:
            interval_size = float(interval_end - interval_start)
            remapped = [(s/interval_size, e/interval_size) for (s, e) in remapped]

    return remapped

@genomic
def interval_distance(query, annot):
    """
    Computes the distance from each query interval to the closest interval in
    annot. If a query interval overlaps the midpoint between two annot intarvals
    it is split into two intervals proximal to each annot interval.    Intervals
    from ``query`` that overlap intervals in ``annot`` are discarded.

    :param query: Data frame with intervals.
    :type query: pandas.DataFrame
    :param annot: Data frame with intervals.
    :type annot: pandas.DataFrame
    :returns: A data frame with remapped intervals.
    :rtype: pandas.DataFrame

    If you want to retain the origianl columns in ``query``, use :any:remap_interval_data.
    """
    return list(chain.from_iterable(remap(q, annot) for q in query))


@genomic
def interval_relative_distance(query, annot):
    """
    Same as :any:interval_relative_distance, but computes the *relative* distance.
    I.e. distances between 0 and 0.5.

    :param query: Data frame with intervals.
    :type query: pandas.DataFrame
    :param annot: Data frame with intervals.
    :type annot: pandas.DataFrame
    :returns: A data frame with remapped intervals.
    :rtype: pandas.DataFrame
    """
    return list(chain.from_iterable(remap(q, annot, relative=True) for q in query))


def remap_interval_data(query, annot):
    """
    Computes the distance from each query interval to the closest interval
    in annot. Original coordinates are preserved as ``'orig_start'`` and
    ``'orig_end'`` columns. If a query interval overlaps the midpoint between two
    annot intarvals it is split into two intervals proximal to each
    annot interval, thus contributing two rows to the returned data frame.
    Intervals from ``query`` that overlap intervals in ``annot`` are discarded.

    :param query: Data frame with intervals.
    :type query: pandas.DataFrame
    :param annot: Data frame with intervals.
    :type annot: pandas.DataFrame
    :returns: A data frame with remapped intervals.
    :rtype: pandas.DataFrame
    """

    annot_grouped = annot.groupby('chrom')

    df_list = list()
    column_names = tuple(query.columns.values)
    for chrom, group in query.groupby('chrom'):

        chrom_annot = annot_grouped.get_group(chrom)
        annot_tups = [tuple(t) for t in chrom_annot[['start', 'end']].itertuples(index=False)]

        remapped = list()
        # for index, row in group.iterrows():            
            # start, end = (row['start'], row['end'])
        for row in group.itertuples(index=False):            
            start, end = row.start, row.end
            for remapped_start, remapped_end, start_prox, end_prox in remap((start, end), annot_tups, include_prox_coord=True):
                remapped.append((remapped_start, remapped_end, start_prox, end_prox) + tuple(row))

        df = pandas.DataFrame().from_records(remapped, 
                columns=('start_remap', 'end_remap', 'start_prox', 'end_prox') + column_names)
            # df = pandas.DataFrame().from_records(remapped, columns=['idx', 'start_remap', 'end_remap']).set_index('idx')
            # df = pandas.merge(group, df, right_index=True, left_index=True)

        df_list.append(df)

    df = (pandas.concat(df_list)
            .reset_index(drop=True)
            .rename(columns={'start': 'start_orig',
                            'end': 'end_orig'})
            .rename(columns={'start_remap': 'start',
                            'end_remap': 'end'})
            )

    df['start_orig'] = df['start_orig'].astype('Int64')
    df['end_orig'] = df['end_orig'].astype('Int64')
    df['start_prox'] = df['start_prox'].astype('Int64')
    df['end_prox'] = df['end_prox'].astype('Int64')

    return df



# def ovl_interval_data(query, annot):

#     query_grouped = query.groupby('chrom')
#     annot_grouped = annot.groupby('chrom')

#     query_df_list = list()
#     annot_df_list = list()

#     for chrom, query_group in query_grouped:
#         annot_group = annot_grouped.get_group(chrom)

#         starts = query_group.start.tolist()
#         ends = query_group.end.tolist()

#         idx_list = list()    
#         annot_idx_list = list()
#         for tup in annot_group.itertuples():

#             start_idx = bisect.bisect_right(starts, tup.start) - 1
#             end_idx = bisect.bisect_right(starts, tup.end) - 1

#             if start_idx > -1 and tup.start < ends[start_idx]:
#                 idx_list.append(start_idx)
#                 annot_idx_list.append(tup.Index)
#             elif start_idx > -1 and tup.end < ends[start_idx]:
#                 idx_list.append(end_idx)
#                 annot_idx_list.append(tup.Index)

#         query_df_list.append(query_group.iloc[idx_list])
#         annot_df_list.append(annot_group.loc[annot_idx_list, ['start', 'end']])

#     query_data_overlap = (pandas.concat(query_df_list)
#                             .reset_index(drop=True)
#                          )
#     annot_intervals = (pandas.concat(annot_df_list)
#                         .reset_index(drop=True)
#                         .rename(columns={'start': 'ovl_start', 'end': 'ovl_end'})
#                         )
#     return pandas.concat([query_data_overlap, annot_intervals], axis=1)


def interval_permute(df, chromosome_sizes):
    """
    Permute intervals not preserving size of gaps.
    """

    group_list = list()
    for chrom, group in df.groupby('chrom'):

        assert group.end.max() <= chromosome_sizes[chrom]

        segment_lengths = group.end - group.start
        total_gap = numpy.sum(group.start - group.end.shift())
        if numpy.isnan(total_gap): # in case there are no internal gaps (one segment)
            total_gap = 0
        else:
            total_gap = int(total_gap)
        if group.start.iloc[0] != 0:
            total_gap += group.start.iloc[0]
        if group.end.iloc[-1] != chromosome_sizes[chrom] + 1:
            total_gap += chromosome_sizes[chrom] + 1 - group.end.iloc[-1]

        assert total_gap >= len(segment_lengths)+1, (total_gap, len(segment_lengths)+1)
        idx = pandas.Series(sorted(random.sample(range(total_gap), len(segment_lengths)+1)))
        gap_lengths = (idx - idx.shift()).dropna()

        borders = numpy.cumsum([j for i in zip(gap_lengths, segment_lengths) for j in i])
        starts, ends = borders[::2], borders[1::2]

        new_df = pandas.DataFrame({'chrom': chrom, 'start': starts, 'end': ends})
        group_list.append(new_df)

    return pandas.concat(group_list)


def bootstrap(chromosome_sizes, samples=1000, smaller=False, return_boot=False):
    """
    Decorator that turns a function producing a statistic into one that also
    produces a p-value from bootstrapping. The bootstrapping resamples the
    intervals of the second argument for each chromosome independently. Only
    required argument to bootstrap is the name of the genome assembly used.

    :param chromosome_sizes: The name of the genome assembly used, or a dictionary mapping chromosome names to chromosome lengths.
    :type chromosome_sizes: str or dict
    :param samples: The number of bootstrap samples to use. Default is 1000. 
    :type samples: int
    :param smaller: Whether to consider extremely small values rather than more extreme large ones.
    :type smaller: bool
    :param return_boot: Whether to return the bootstrap samples too.
    :type return_boot: bool

    :returns: statistic, p-value With ``return_boot=True`` a list of bootstrap values are returned as well.
    :rtype: float, float, [list]

    :raises: KeyError if ``chromosome_sizes`` string is not an assembly known to the module. In that case specify as dictionary.
    """
    if type(chromosome_sizes) is str:
        chromosome_sizes = chrom_sizes[chromosome_sizes]

    def decorator(func):
        @wraps(func)
        def wrapper(query, annot, **kwargs):

            stat = func(query, annot, **kwargs)

            boot = list()
            for i in range(samples):
                perm = interval_permute(query, chromosome_sizes)
                boot.append(func(perm, annot, **kwargs))
            boot.sort()
            if smaller:
                p_value = bisect.bisect_left(boot, stat) / len(boot)
            else:
                p_value = (len(boot) - bisect.bisect_left(boot, stat)) / len(boot)
            if p_value == 0:
                sys.stderr.write('p-value is zero smaller than {}. Increase nr samples to get actual p-value.\n'.format(1/samples))

            if return_boot:
                return stat, p_value, boot
            else:
                return stat, p_value

        return wrapper
    return decorator


def proximity_test(query, annot, samples=10000, npoints=1000):
    """This function does something.

    :param name: The name to use.
    :type name: str.
    :param state: Current state to be in.
    :type state: bool.
    :returns:  int -- the return code.
    :raises: AttributeError, KeyError

    """
    remapped_df = interval_relative_distance(query, annot)
    distances = remapped_df.start

    def _stat(distances, npoints):
        obs_ecdf = ECDF(distances)
        points = numpy.linspace(0, 0.5, num=npoints)    
        test_stat = sum(obs_ecdf(points) - 2 * points) * 2 / npoints
        return test_stat
    
    test_stat = _stat(distances, npoints)
    
    null_distr = list()
    for i in range(samples):
        sampled_distances = numpy.random.uniform(0, 0.5, len(distances))
        # we compute absolute values for the null distribution
        null_distr.append(abs(_stat(sampled_distances, npoints)))
    null_distr.sort()
    # we compare abs value of test_stat to null distr to get pvalue for extremity in both directions
    p_value = (len(null_distr) - bisect.bisect_left(null_distr, abs(test_stat))) / len(null_distr)

    TestResult = namedtuple('TestResult', ['statistic', 'pvalue'])
    return TestResult(test_stat, p_value)


def jaccard(a, b):
    """
    Compute Jaccard overlap statistic

    :param name: The name to use.
    :type name: str.
    :param state: Current state to be in.
    :type state: bool.
    :returns:  int -- the return code.
    :raises: AttributeError, KeyError
    """
    inter = interval_intersect(a, b)
    union = interval_union(a, b)

    return sum(inter.end - inter.start) / sum(union.end - union.start)


if __name__ == "__main__":


    print(remap((300, 400), [(0, 100), (1000, 1100)], include_prox_coord=True))

    assert 0

    query = pandas.DataFrame(dict(chrom='X', start=[3, 5], end=[15, 7], extra=['foo', 'bar']))
    print(query)
    annot = pandas.DataFrame(dict(chrom='X', start=[1, 20], end=[2, 25]))
    print(annot)
    print(remap_interval_data(query, annot))

    assert 0

    query = pandas.DataFrame(dict(chrom='X', start=[3, 5], end=[22, 7], extra=['foo', 'bar']))
    print(query)
    annot = pandas.DataFrame(dict(chrom='X', start=[1, 20], end=[2, 25]))
    print(annot)

    print(ovl_interval_data(query, annot))

    assert 0

    # annotation
    tp = [('chr1', 1, 3), ('chr1', 4, 10), ('chr1', 25, 30), ('chr1', 20, 27), ('chr2', 1, 10), ('chr2', 1, 3)]
    annot = pandas.DataFrame.from_records(tp, columns=['chrom', 'start', 'end'])
    print("annot\n", annot)

    # query
    tp = [('chr1', 8, 22), ('chr8', 14, 15)]
    query = pandas.DataFrame.from_records(tp, columns=['chrom', 'start', 'end'])
    print("query\n", query)

    annot_collapsed = interval_collapse(annot)
    print("annot_collapsed\n", annot_collapsed)

    non_ovl_query = interval_diff(query, annot_collapsed)
    print("non_ovl_query\n", non_ovl_query)

    distances = interval_distance(non_ovl_query, annot_collapsed)
    print("distances\n", distances)

    print("distance test\n", proximity_test(non_ovl_query, annot_collapsed))


    sys.exit()
    # ##################################################################

    print('jaccard:')
    annot = pandas.DataFrame({'chrom': 'chr1', 'start': range(0, 1000000, 1000), 'end': range(100, 1000100, 1000)})
    query = pandas.DataFrame({'chrom': 'chr1', 'start': range(50, 1000050, 1000), 'end': range(150, 1000150, 1000)})

    print(annot)
    print(query)

    # print(interval_jaccard(query, annot, samples=10, chromosome_sizes={'chr1': 1500000, 'chr2': 1500000}))


    annot = pandas.DataFrame({'chrom': 'chr1', 'start': [500, 2000], 'end': [1000, 2500]})
    query = pandas.concat([pandas.DataFrame({'chrom': 'chr1', 'start': [1100, 1800], 'end': [1200, 1900]}),
                           pandas.DataFrame({'chrom': 'chr2', 'start': [1100, 1800], 'end': [1200, 1900]})
                           ])

    print(annot)
    print(query)

    chromosome_sizes={'chr1': 1500000, 'chr2': 1500000}
    @bootstrap(chromosome_sizes, samples=10)
    def my_stat(a, b):
        return jaccard(a, b)

    print(my_stat(query, annot))

    @genomic
    def interval_computation(a, b):
        # stupid function that always returns:
        return [(1, 1)]

    print(interval_computation(query, annot))

    @bootstrap('hg19', samples=10)
    def statistic(a, b):
        return 42

    print(statistic(query, annot))

    # def overlap_count(query, annot):
    #       center = annot.start + (annot.end-annot.start)/2
    #       b = numpy.equal(numpy.searchsorted(annot.start, center) - 1,
    #                       numpy.searchsorted(annot.end, center, side='left')) & \
    #             (query.chrom == annot.chrom)

    #       return b.sum()

    # @bootstrap(chromosome_sizes, samples=1000)
    # def my_stat(a, b):
    #     return overlap_count(a, b)

    # print(my_stat(query, annot))

    ####################################################

    # annot_collapsed = (annot.groupby('chrom')
    #                    .apply(interval_collapse)
    #                    .reset_index(drop=True)
    #                    )
    # print("annot_collapsed\n", annot_collapsed)

    # non_ovl_query = (DataFrameList(query, annot_collapsed)
    #                  .groupby('chrom')
    #                  .apply(interval_diff)
    #                  .reset_index(drop=True)
    #                  )
    # print("non_ovl_query\n", non_ovl_query)

    # distances = (DataFrameList(non_ovl_query, annot_collapsed)
    #        .groupby('chrom')
    #        .apply(interval_distance)
    #        .reset_index(drop=True)
    #        )
    # print("distances\n", distances)

