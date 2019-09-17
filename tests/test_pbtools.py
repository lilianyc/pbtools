
"""Tests for pbtools."""

import numpy as np
import pandas as pd
import pbxplore as pbx
import pytest

import pbtools as pbt


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


def test_mutual_information():

    for data in get_MI_matrix():
        sequences = data[0]
        expected = data[1]
        assert pbt.mutual_information_matrix(sequences) == pytest.approx(expected)