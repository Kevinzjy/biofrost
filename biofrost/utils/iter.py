"""Other useful functions"""
import itertools

__all__ = [
    "flatten",
    "grouper",
    "empty_iter",
    "pairwise",
]


def flatten(x):
    """Flatten list of lists"""
    flatted_list = list(itertools.chain(*x))
    return flatted_list


def grouper(iterable, n, fillvalue=None):
    """
    Collect data info fixed-length chunks or blocks
    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    """
    from itertools import zip_longest
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=None)


def empty_iter(iterable):
    """
    Return None if the iterable object is empty
    """
    try:
        first = next(iterable)
    except StopIteration:
        return None
    return itertools.chain([first], iterable)


def pairwise(iterable):
    """"s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    from itertools import tee
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)
