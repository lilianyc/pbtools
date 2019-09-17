
"""Main module of pbtools

"""

import argparse
import itertools
import logging
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pbxplore as pbx


#print(globals().get("__file__"))
#DATA_DIR = Path(__file__).parent.resolve().joinpath("data")
#DATA_DIR = Path().resolve().parent.joinpath("data")

PB_NAMES = 'abcdefghijklmnop'
DEBUG = True

# Set custom logger.
log = logging.getLogger(__name__)
log_handler = logging.StreamHandler()
log_formatter = logging.Formatter("[%(asctime)s] %(levelname)-8s: %(message)s",
                                  "%Y-%m-%d, %H:%M:%S")
log_handler.setFormatter(log_formatter)
log.addHandler(log_handler)
log.setLevel(logging.DEBUG)



trajectory = "psi_md_traj.xtc"
topology = "psi_md_traj.gro"

'''

def user_inputs():
    """Parse and handle the submitted command line.

    Parameters
    ----------
    args : list of str
        List of arguments received from the CLI.

    Returns
    -------
    argparse.Namespace
        Object containing the arguments parsed from the CLI.

    Raises
    ------
    SystemExit
        If the file provided is not found.

    Notes
    -----
    | Heavily based upon PBXplore's user_input model at:
    | https://github.com/pierrepo/PBxplore/blob/master/pbxplore/scripts/
      PBassign.py

    """
    parser = argparse.ArgumentParser(
        description='Read PDB structures and assign protein blocs (PBs).')

    # arguments
    parser.add_argument("-p", action="append",
                        help=("name of a pdb file "
                              "or name of a directory containing pdb files"))
    parser.add_argument("-o", action="store", required=True,
                        help="name for results")
    # arguments for MDanalysis
    group = parser.add_argument_group(
        title='other options to handle molecular dynamics trajectories')
    group.add_argument("-x", action="store", metavar='TRAJECTORY',
                       help="name of the trajectory file")
    group.add_argument("-g", action="store", metavar='TOPOLOGY',
                       help="name of the topology file")

    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {}'.format(__version__))
    # get all arguments
    options = parser.parse_args()

    # check options
    if not options.p:
        if not options.x:
            parser.print_help()
            parser.error("use at least option -p or -x")
        elif not options.g:
            parser.print_help()
            parser.error("option -g is mandatory, with use of option -x")

    # check files
    pdb_name_lst = []
    if options.p:
        for name in options.p:
            # input is a file: store file name
            if Path(name).isfile():
                pdb_name_lst.append(name)
            # input is a directory: list and store all PDB and PDBx/mmCIF files
            elif Path(name).isdir():
#                pdb_name_lst = list(Path().glob("sequence*.pdb"))
#                for extension in (pbx.structure.PDB_EXTENSIONS + pbx.structure.PDB.PDBx_EXTENSIONS):
#                    pdb_name_lst += glob.glob(os.path.join(name, "*" + extension))
                pass
            # input is neither a file nor a directory.
            elif (not Path(name).isfile() or not Path(name).isdir()):
                parser.error("{0}: not a valid file or directory".format(name))
    else:
        if not Path(options.x).isfile():
            sys.exit("{0}: not a valid file".format(options.x))
        elif not Path(options.g).isfile():
            sys.exit("{0}: not a valid file".format(options.g))

    return options, pdb_name_lst


def cli(args=None):
    """Entry point for seq_to_first_iso's CLI.

    Parameters
    ----------
    args : list of str, optional
        CLI arguments, args are used for testing (default is None for CLI).
    
    Returns
    -------
    None
        Writes a tsv file.

    Raises
    ------
    SystemExit
        If no sequences were found on the file.

    Notes
    -----
    Main function of the script, for use with CLI.

    """
    if not args:
        args = sys.argv[1:]

    options, pdb_name_lst = user_inputs()
    if options.p:
        if pdb_name_lst:
            print("{} PDB file(s) to process".format(len(pdb_name_lst)))
        else:
            print('Nothing to do. Good bye.')
            return
        # PB assignement of PDB structures
        chains = pbx.chains_from_files(pdb_name_lst)
    else:
        # PB assignement of a Gromacs trajectory
        chains = pbx.chains_from_trajectory(options.x, options.g)

    all_comments = []
    all_sequences = []

    for comment, chain in chains:
        try:
            dihedrals = chain.get_phi_psi_angles()
            sequence = pbx.assign(dihedrals)
            all_comments.append(comment)
            all_sequences.append(sequence)
        except FloatingPointError:
            log.error("The computation of angles produced NaN. "
                      "This typically means there are issues with "
                      "some residues coordinates. "
                      f"Check your input file ({comment})")

'''

def mutual_information(pos1, pos2):
    """Computes the MI from 2 Series.

    Parameters
    ----------
    pos1 : pandas.Series
        Protein Blocks of a given sequence position.
    pos2 : pandas.Series
        Protein Blocks of a given sequence position.

    Returns
    -------
    MI : float
        Mutual Information of the 2 positions

    Raises
    ------
    AssertionError
        If the Series have different lengths

    """

    assert len(pos1) == len(pos2), "Series have different lengths"
    MI = 0

    for PB1, PB2 in itertools.permutations(PB_NAMES, 2):
        PB1_count = (pos1 == PB1).sum()
        PB2_count = (pos2 == PB2).sum()
        joint_prob = ((pos1 == PB1) &
              (pos2 == PB2)).sum() / PB1_count
        # Denominator or joint probability are 0.
        if not (PB1_count and PB2_count and joint_prob):
            continue
        PB1_prob = PB1_count / len(pos1)
        PB2_prob = PB2_count / len(pos1)

        log.debug(f"{joint_prob}, {PB1_prob}, {PB2_prob}")

        MI += joint_prob + math.log((joint_prob/(PB1_prob * PB2_prob)),
                                    len(PB_NAMES))
    return MI


def mutual_information_matrix(sequences):
    """Return the matrix of Mutual Informations for all positions.

    Parameters
    ----------
    sequences : list of str
        Sequences obtained from the file(s)

    Returns
    -------
    MI_matrix : numpy.ndarray
        Matrix of Mutual informations.

    """
    # Create dataframe from a generator.
    df_seq = pd.DataFrame((list(seq) for seq in sequences))
    positions = len(sequences[0])
    MI_matrix = np.zeros((positions, positions))
    # Get all combinations of positions.
    for pos1, pos2 in itertools.combinations(range(positions), 2):
        log.debug(f"{pos1}, {pos2}")
        MI = mutual_information(df_seq[pos1], df_seq[pos2])
        MI_matrix[pos1, pos2] = MI

    return MI_matrix



if __name__ == "__main__":
    if DEBUG:
        sequences = []
        try:
            for chain_name, chain in pbx.chains_from_trajectory(trajectory,
                                                                topology):
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
    
        
        c = mutual_information_matrix(small_seq)
        
        for i in itertools.combinations(PB_NAMES, 2):
            print(i, end=" ") 
        
    
        df_seq = pd.DataFrame((list(seq) for seq in small_seq))
        
        
        b = mutual_information(df_seq[2], df_seq[3])


