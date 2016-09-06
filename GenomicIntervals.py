import pandas as pd

from functools import wraps
import bisect, math
from itertools import chain, izip


def chromosomal(strip=True, addback=True):
    """decorator for converting a function operating on (start, end) tuples 
    to one that takes data frames with chrom, start, end columns"""
    def decorator(func):
        @wraps(func)
        def wrapper(*args):
            chrom_set = set()
            tps_list = list()
            for df in args:
                chrom_set.update(df['chrom'])
                if strip:
                    tps = sorted(zip(df['start'], df['end']))
                else:
                    tps = sorted(zip(df['chrom'], df['start'], df['end']))    
                tps_list.append(tps)
            assert len(chrom_set) == 1
            chrom = chrom_set.pop()
            res_df = pd.DataFrame.from_records(func(*tps_list), columns = ['start', 'end'])
            if addback:
                res_df['chrom'] = chrom
            return res_df
        return wrapper
    return decorator    


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

@chromosomal()
def interval_diff(a, b):
    return merge(a, b, lambda in_a, in_b: in_a and not in_b)

@chromosomal()
def interval_union(a, b):
    return merge(a, b, lambda in_a, in_b: in_a or in_b)

@chromosomal()
def interval_intersect(a, b):
    return merge(a, b, lambda in_a, in_b: in_a and in_b)

@chromosomal()
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

@chromosomal()
def interval_distance(query, annot):
    """
    Convertes absolute coordinates of each query intervals so that start and end
    to distances relative to the most proximal annotation.  
    If an interval overlaps the midpoint between two annotations it is split 
    into two intervals proximal to each annotation.
    It is assumed that the query intervals do not overlap the 
    """
    return list(chain.from_iterable(remap(q, annot) for q in query))



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
            subframes = [df.loc[gr[by]] for df, gr in izip(self.dflist.frames, self.groupings)]
            group_results.append(fun(*subframes))
        return pd.concat(group_results)


if __name__ == "__main__":

	# In all of the following, the list of intervals must be sorted and 
	# non-overlapping. We also assume that the intervals are half-open, so
	# that x is in tp(start, end) iff start <= x and x < end.

	# annotation
	tp = [('chr1', 1, 3), ('chr1', 4, 10), ('chr1', 25, 30), ('chr1', 20, 27), ('chr2', 1, 10), ('chr2', 1, 3)]
	annot = pd.DataFrame.from_records(tp, columns=['chrom', 'start', 'end'])
	print "annot\n", annot

	# query
	tp = [('chr1', 8, 22), ('chr2', 14, 15)]
	query = pd.DataFrame.from_records(tp, columns=['chrom', 'start', 'end'])
	print "query\n", query

	annot_collapsed = (annot.groupby('chrom')
	                   .apply(interval_collapse)
	                   .reset_index(drop=True)
	                   )
	print "annot_collapsed\n", annot_collapsed

	non_ovl_query = (DataFrameList(query, annot_collapsed)
	                 .groupby('chrom')
	                 .apply(interval_diff)
	                 .reset_index(drop=True)
	                 )
	print "non_ovl_query\n", non_ovl_query

	distances = (DataFrameList(non_ovl_query, annot_collapsed)
	       .groupby('chrom')
	       .apply(interval_distance)
	       .reset_index(drop=True)
	       )
	print "distances\n", distances





