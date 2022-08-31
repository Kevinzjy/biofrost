from .io import is_gzipped, check_dir, check_file
from .text import to_bytes, to_str
from .iter import flatten, grouper, empty_iter, pairwise
from .intervals import cluster_positions, merge_intervals

__all__ = [
    "is_gzipped",
    "check_file",
    "check_dir",

    "to_bytes",
    "to_str",

    "flatten",
    "grouper",
    "empty_iter",
    "pairwise",

    "cluster_positions",
    "merge_intervals",
]