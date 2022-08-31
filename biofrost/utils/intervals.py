"""Functions for processing intervals"""

__all__ = [
    "cluster_positions",
    "merge_intervals",
]


def cluster_positions(x, gap=5):
    """Cluster positions into intervals

    Args:
        x (list): List of genomic positions, does not required to be sorted
        gap (int, optional): Minimum distance to consider as gap. Defaults to 5.

    Returns:
        list: [(st, en), ], list of clustered intervals
    """
    if len(x) == 0:
        return []
    x_sorted = sorted(x)
    clustered = [[x_sorted[0], x_sorted[0]], ]
    for i in x_sorted[1:]:
        if i-clustered[-1][1] > gap:
            clustered.append([i,i])
        else:
            clustered[-1][1] = i
    return clustered


def merge_intervals(blocks, gap=5):
    """Merge small blocks into large intervals

    Args:
        blocks (list): Input list of (start, end), does not required to be sorted
        gap (int, optional): Minimum gap to be considered as different interval. Defaults to 5.

    Returns:
        [(start, end)]: Merged intervals
    """
    from operator import itemgetter
    tmp = sorted(blocks, key=itemgetter(0, 1))
    merged = []
    last_st, last_en = tmp[0][0], tmp[0][1]
    for x in tmp[1:]:
        st, en = x[0], x[1]
        if st <= last_en + gap:
            last_en = max(en, last_en)
            last_st = min(st, last_st)
        else:
            merged.append([last_st, last_en])
            last_st, last_en = st, en
    merged.append([last_st, last_en])
    return merged