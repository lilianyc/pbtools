
"""Tests for pbtools."""

import networkx as nx
import numpy as np
import pandas as pd
import pbxplore as pbx
import pytest

import pbtools as pbt

# Precision of float approximations.
REL = 1e-6

#TODO: add more tests on networks + on pbxplore.

def test_mutual_information():
    assert pbt.mutual_information(pd.Series(["a", "b"]), pd.Series(["a", "c"])) \
           == pytest.approx(0.25, REL)

    with pytest.raises(AssertionError):
        pbt.mutual_information(1, 2)
        pbt.mutual_information([], [])
        pbt.mutual_information(pd.Series([]), pd.Series([1]))


@pytest.fixture(scope="session")
def get_MI_matrix():
    return [(["aaa", "cab"], np.array([[0, 0, 0.25],
                                      [0, 0, 0],
                                      [0, 0, 0]])),
            (["zzaaazz", "zzcabzz"], np.array([[0, 0, 0, 0, 0, 0, 0],
                                              [0, 0, 0, 0, 0, 0, 0],
                                              [0, 0, 0, 0, 0.25, 0, 0],
                                              [0, 0, 0, 0, 0, 0, 0],
                                              [0, 0, 0, 0, 0, 0, 0],
                                              [0, 0, 0, 0, 0, 0, 0],
                                              [0, 0, 0, 0, 0, 0, 0]])),
            ]


def test_mutual_information_matrix(get_MI_matrix):

    for data in get_MI_matrix:
        sequences = data[0]
        expected = data[1]
        assert pbt.mutual_information_matrix(sequences) == pytest.approx(expected)


def test_network_creation():

    with pytest.raises(AssertionError):
        pbt.interaction_graph(np.array([[0, 0, 1],
                                        [0, 0, 0]]))
