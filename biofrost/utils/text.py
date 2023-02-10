"""Functions for dealing with string / bytes"""
import os
from pathlib import Path
# import sys
# import re
# import itertools
# import time

__all__ = [
    "to_str", "to_bytes", "generate_random_key", "mk_random_dir",
]


def to_str(bytes_or_str):
    """Return Instance of str"""
    if isinstance(bytes_or_str, bytes):
        value = bytes_or_str.decode('utf-8')
    else:
        value = bytes_or_str
    return value


def to_bytes(bytes_or_str):
    """Return Instance of bytes"""
    if isinstance(bytes_or_str, str):
        value = bytes_or_str.encode('utf-8')
    else:
        value = bytes_or_str
    return value


def generate_random_key(length):
    """Generate random key of specific length

    Args:
        length (int): length

    Returns:
        str: random string
    """
    import random
    import string
    return ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(length))


def mk_random_dir(tmp_dir=None):
    """Create random directory

    Args:
        tmp_dir (str, optional): path of temporary directory. Defaults to /tmp or TMP_DIR.

    Returns:
        pathlib.Path: Path object
    """
    tmp_dir = None
    if tmp_dir is None and 'TMP_DIR' in os.environ:
        tmp_dir = os.environ['TMP_DIR']
    else:
        tmp_dir = "/tmp"

    root = Path(tmp_dir) / f"pyblast_{generate_random_key(6).upper()}"
    while root.exists():
        root = Path(tmp_dir) / f"pyblast_{generate_random_key(6).upper()}"
    root.mkdir()

    return root


# def ranking(ranks, names, order=1):
#     """Ranking a list using another list as key

#     Args:
#         ranks (list): list of keys for ranking
#         names (list): list of labels corresponding to each key in ranks
#         order  (list): 1 or -1, if order equals -1, the order is reversed

#     Returns:
#         dict: A dict of rank

#     """
#     import numpy as np
#     from sklearn.preprocessing import MinMaxScaler
#     minmax = MinMaxScaler()
#     ranks = minmax.fit_transform(order * np.array([ranks]).T).T[0]
#     ranks = map(lambda x: round(x, 2), ranks)
#     return dict(zip(names, ranks))




# def tree():
#     """Tree structure"""
#     from collections import defaultdict
#     return defaultdict(tree)




# def hash_str(f):
#     """Return sha256 of input string"""
#     import hashlib
#     return hashlib.sha256(str(f).encode()).hexdigest()


# def sorted_iters(iters, key, reverse=False):
#     """
#     Return iters with minimum values of given keys
#     Parameters
#     ----------
#     iters : iterable object
#         Iterable object
#     key : key
#         key to sort
#     reverse : boolean
#         if True, return maximum values
#     Returns
#     -------
#     str
#         Absolute path of directory
#     """
#     from operator import itemgetter
#     x = sorted(iters, key=itemgetter(key), reverse=reverse)
#     return [i for i in x if i[key] == x[0][key]]


# def download(url, outfile):
#     """Download file"""
#     import requests
#     from contextlib import closing
#     from .logger import ProgressBar

#     with closing(requests.get(url, stream=True)) as response, open(outfile, 'wb') as out:
#         chunk_size = 1024
#         try:
#             content_size = int(response.headers['content-length'])
#             sys.stderr.write('Downloading file to {}, total size: {}\n'.format(outfile, content_size))
#             prog = ProgressBar()
#             cnt = 0
#             for data in response.iter_content(chunk_size=chunk_size):
#                 out.write(data)
#                 cnt += len(data)
#                 prog.update(100 * cnt / content_size)
#         except KeyError:
#             out.write(response.content)
#     return 1


