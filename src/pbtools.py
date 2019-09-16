
"""

"""

from pathlib import Path

import pbxplore as pbx

#print(globals().get("__file__"))
#DATA_DIR = Path(__file__).parent.resolve().joinpath("data")
#DATA_DIR = Path().resolve().parent.joinpath("data")

trajectory = "md.trr"
topology = "md.gro"

try:
    for chain_name, chain in pbx.chains_from_trajectory(trajectory, topology):
        dihedrals = chain.get_phi_psi_angles()
        pb_seq = pbx.assign(dihedrals)
    
    print(chain_name)
    print(pb_seq)
except Exception as exc:
    print("Error", exc)




