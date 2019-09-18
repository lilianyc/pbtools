
"""Perform analysis on protein sequences encoded as Protein Blocks

Provide functions to compute a Mutual Information matrix from sequences.

"""

__authors__ = "Lilian Yang-crosson"
__license__ = "BSD 3-Clause License"
__version__ = "0.1.0"
__maintainer__ = "Lilian Yang-crosson"


from .pbtools import (PB_NAMES,
                      MISSING_BLOCK,
                      interaction_graph,
                      mutual_information,
                      mutual_information_matrix)

__all__ = ["PB_NAMES",
           "MISSING_BLOCK",
           "interaction_graph",
           "mutual_information",
           "mutual_information_matrix"]
