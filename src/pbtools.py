
"""

"""

from pathlib import Path

import pandas as pd
import pbxplore as pbx

#print(globals().get("__file__"))
#DATA_DIR = Path(__file__).parent.resolve().joinpath("data")
#DATA_DIR = Path().resolve().parent.joinpath("data")

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
