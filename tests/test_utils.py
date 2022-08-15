import sys
from biofrost import utils


def test_ranking():
    print(utils.ranking([1, 2, 3], ['A', 'B', 'C']))
