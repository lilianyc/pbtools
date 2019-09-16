
"""Main module of pbtools

"""

import argparse
import itertools
import logging
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pbxplore as pbx

#print(globals().get("__file__"))
#DATA_DIR = Path(__file__).parent.resolve().joinpath("data")
#DATA_DIR = Path().resolve().parent.joinpath("data")

PB_NAMES = 'abcdefghijklmnop'

# Set custom logger.
log = logging.getLogger(__name__)
log_handler = logging.StreamHandler()
log_formatter = logging.Formatter("[%(asctime)s] %(levelname)-8s: %(message)s",
                                  "%Y-%m-%d, %H:%M:%S")
log_handler.setFormatter(log_formatter)
log.addHandler(log_handler)
log.setLevel(logging.DEBUG)



trajectory = "md.trr"
topology = "md.gro"

sequences = []
try:
    for chain_name, chain in pbx.chains_from_trajectory(trajectory, topology):
        dihedrals = chain.get_phi_psi_angles()
        pb_seq = pbx.assign(dihedrals)
        sequences.append(pb_seq)
except Exception as exc:
    print("Error", exc)

# Count matrix with one row per sequence and one column per PB.
# The readthedocs documentation has the 2 confused.
# !!!: Duplicate sequences seem grouped together.
count_matrix = pbx.analysis.count_matrix(sequences)

count_matrix.shape

df = pd.DataFrame(count_matrix)
# Probability matrix.
proba_df = df / len(sequences)
# See if all probabilities are < 1.
assert (proba_df < 1).all().all()

small_seq = sequences[:10]


def mutual_information_matrix(sequences):
    """
    """
    df_seq = pd.DataFrame((list(seq) for seq in sequences))
    positions = len(sequences[0])
    MI_matrix = np.zeros((positions, positions))
#    positions = range(len(sequences[0]))
    # Get all permutations of positions.
    for pos1, pos2 in itertools.permutations(range(positions), 2):
        # for col in range(len(sequences))
        # sequences[pos1] seq[pos2]
#        MI = 0 # mutual_information(pos1, pos2)
#        for PB1, PB2 in itertools.combinations(PB_NAMES, 2):
#            pass
        log.debug(f"{pos1}, {pos2}")
        MI = mutual_information(df_seq[pos1], df_seq[pos2])
        MI_matrix[pos1, pos2] = MI

    return MI_matrix


def mutual_information(pos1, pos2):
    """
    """

    assert len(pos1) == len(pos2), "Series have different lengths"
    MI = 0

    for PB1, PB2 in itertools.combinations(PB_NAMES, 2):
        PB1_count = (pos1 == PB1).sum()
        PB2_count = (pos2 == PB2).sum()
        joint_prob = ((pos1 == PB1) &
              (pos2 == PB2)).sum() / PB1_count
        # Denominator is 0
        if not (PB1_count and PB2_count and joint_prob):
            continue
        PB1_prob = PB1_count / len(pos1)
        PB2_prob = PB2_count / len(pos1)

        log.debug(f"{joint_prob}, {PB1_prob}, {PB2_prob}")

        MI += joint_prob + math.log((joint_prob/(PB1_prob * PB2_prob)),
                                    len(PB_NAMES))
    return MI

c = mutual_information_matrix(small_seq)

for i in itertools.combinations(PB_NAMES, 2):
    print(i, end=" ") 


df_seq = pd.DataFrame((list(seq) for seq in small_seq))


b = mutual_information(df_seq[2], df_seq[3])





