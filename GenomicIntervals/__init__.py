import pandas
import numpy

from functools import wraps, reduce, partial
import bisect, math, random, sys
from itertools import chain

# def with_chrom(strip=True, addback=True):
#     """decorator for converting a function operating on (start, end) tuples 
#     to one that takes data frames with chrom, start, end columns"""
#     def decorator(func):
#         @wraps(func)
#         def wrapper(*args):
#             chrom_set = set()
#             tps_list = list()
#             for df in args:
#                 chrom_set.update(df['chrom'])
#                 if strip:
#                     tps = sorted(zip(df['start'], df['end']))
#                 else:
#                     tps = sorted(zip(df['chrom'], df['start'], df['end']))    
#                 tps_list.append(tps)
#             assert len(chrom_set) == 1
#             chrom = chrom_set.pop()
#             res_df = pandas.DataFrame.from_records(func(*tps_list), columns = ['start', 'end'])
#             if addback:
#                 res_df['chrom'] = chrom
#             return res_df
#         return wrapper
#     return decorator    


def by_chrom(func):
    """decorator to converting a function that operates on a data frame
    with only one chromosome to one operating on a data frame with many chromosomes."""
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
    """decorator for converting a function operating on (start, end) tuples 
    to one that takes data frames with chrom, start, end columns"""
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


# In all of the following, the list of intervals must be sorted and 
# non-overlapping. We also assume that the intervals are half-open, so
# that x is in tp(start, end) iff start <= x and x < end.

def flatten(list_of_tps):
    """Convert a list of intervals to a list of endpoints"""
    return reduce(lambda ls, ival: ls + list(ival), list_of_tps, [])


def unflatten(list_of_endpoints):
    """Convert a list of endpoints, with an optional terminating sentinel,
    into a list of intervals"""
    return [ [list_of_endpoints[i], list_of_endpoints[i + 1]]
          for i in range(0, len(list_of_endpoints) - 1, 2)]


def merge(a_tps, b_tps, op):
    """Merge two lists of intervals according to the boolean function op"""
    a_endpoints = flatten(a_tps)
    b_endpoints = flatten(b_tps)

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


@by_chrom
@with_chrom
def interval_diff(a, b):
    return merge(a, b, lambda in_a, in_b: in_a and not in_b)


@by_chrom
@with_chrom
def interval_union(a, b):
    return merge(a, b, lambda in_a, in_b: in_a or in_b)


@by_chrom
@with_chrom
def interval_intersect(a, b):
    return merge(a, b, lambda in_a, in_b: in_a and in_b)


@by_chrom
@with_chrom
def interval_collapse(a):
    a_un = [list(a[0])]
    for i in range(1, len(a)):
        x = a[i]
        if a_un[-1][1] < x[0]:
            a_un.append(list(x))
        else:
            a_un[-1][1] = x[1]
    return a_un


def remap(query, annot):
    query_start, query_end = query
    annot_starts, annot_ends = zip(*annot)

    # find interval betweent two annotations
    idx = bisect.bisect_right(annot_ends, query_start)
    interval_start = idx != 0 and annot_ends[idx-1] or None
    interval_end = idx != len(annot) and annot_starts[idx] or None
    assert interval_start is not None or interval_end is not None
    
    if interval_start is None:
        interval_mid = - float('inf')
    elif interval_end is None:
        interval_mid = float('inf')
    else:
        assert query_end <= interval_end, "query and annotations coordinates overlap"
        interval_mid = int(interval_end + math.ceil((interval_start - interval_end) / 2.0))

    if interval_mid < query_start:
        remapped = [(interval_end - query_end, interval_end - query_start)]
    elif interval_mid >= query_end:
        remapped = [(query_start - interval_start, query_end - interval_start)]
    else:
        remapped = [(query_start - interval_start, interval_mid - interval_start),
                    (interval_end - query_end, interval_end - interval_mid)]
    return remapped


@by_chrom
@with_chrom
def interval_distance(query, annot):
    """
    Convertes absolute coordinates of each query intervals so that start and end
    to distances relative to the most proximal annotation.  
    If an interval overlaps the midpoint between two annotations it is split 
    into two intervals proximal to each annotation.
    It is assumed that the query intervals do not overlap the 
    """
    return list(chain.from_iterable(remap(q, annot) for q in query))



def jaccard_stat(a, b, chromosome_sizes, permute=False):
    """
    Compute Jaccard overlap statistic, optionally 
    preceeded by permuting intervals in first argument.
    """

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

    if permute:
        a = interval_permute(a, chromosome_sizes)

    inter = interval_intersect(a, b)
    union = interval_union(a, b)

    return sum(inter.end - inter.start) / float(sum(union.end - union.start))



def interval_jaccard(query, annot, samples=1000, chromosome_sizes={}, dview=None):
    """
    Compute jaccard test statistic and p-value.
    """

    # compute actual jaccard stat for query and annot
    test_stat = jaccard_stat(query, annot)

    # partial function for jaccard stat with permuted query intervals
    jaccard_stat_perm = partial(jaccard_stat, chromosome_sizes=chromosome_sizes, permute=True)

    # sampling of jaccard stats with permuted query intervals
    if dview is not None:
        assert 0, "this does not work..., don't use dview..."
        null_distr = list(dview.map_sync(jaccard_stat_perm, [query] * samples, [annot] * samples))
    else:
        null_distr = list(map(jaccard_stat_perm, [query] * samples, [annot] * samples))

    # get p value
    null_distr.sort()
    p_value = (len(null_distr) - bisect.bisect_left(null_distr, test_stat)) / float(len(null_distr))
    if p_value == 0:
        sys.stderr.write('p-value is zero smaller than {}. Increase nr samples to get actual p-value.\n'.format(1/float(samples)))

    return test_stat, p_value


class DataFrameList(object):
    
    def __init__(self, *args):
        self.frames = args
    
    def groupby(self, by):
        return _MultiGroupBy(self, by)
    

class _MultiGroupBy(object):

    def __init__(self, dflist, by):
        self.dflist = dflist
        self.groupings = [df.groupby(by).groups for df in dflist.frames]
    
    def apply(self, fun):
        group_results = list()
        for by in self.groupings[0].keys():
            subframes = [df.loc[gr[by]] for df, gr in zip(self.dflist.frames, self.groupings)]
            group_results.append(fun(*subframes))
        return pandas.concat(group_results)



if __name__ == "__main__":
    # annotation
    tp = [('chr1', 1, 3), ('chr1', 4, 10), ('chr1', 25, 30), ('chr1', 20, 27), ('chr2', 1, 10), ('chr2', 1, 3)]
    annot = pandas.DataFrame.from_records(tp, columns=['chrom', 'start', 'end'])
    print("annot\n", annot)

    # query
    tp = [('chr1', 8, 22), ('chr2', 14, 15)]
    query = pandas.DataFrame.from_records(tp, columns=['chrom', 'start', 'end'])
    print("query\n", query)

    annot_collapsed = interval_collapse(annot)
    print("annot_collapsed\n", annot_collapsed)

    non_ovl_query = interval_diff(query, annot_collapsed)
    print("non_ovl_query\n", non_ovl_query)

    distances = interval_distance(non_ovl_query, annot_collapsed)
    print("distances\n", distances)

    print('orig:')
    print(query)

    # print('permuted:')
    # print(interval_permute(query, {'chr1': 50, 'chr2': 50}))

    print('jaccard:')
    annot = pandas.DataFrame({'chrom': 'chr1', 'start': range(0, 1000000, 1000), 'end': range(100, 1000100, 1000)})
    query = pandas.DataFrame({'chrom': 'chr1', 'start': range(50, 1000050, 1000), 'end': range(150, 1000150, 1000)})

    print(annot.head())
    print(query.head())

    print(interval_jaccard(query, annot, samples=10, chromosome_sizes={'chr1': 1500000, 'chr2': 1500000}))

    # FIXME: it does not work if query has chromosoems not in annotation...





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

