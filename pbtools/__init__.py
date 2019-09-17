
"""Perform analysis on protein sequences encoded as Protein Blocks

Provide functions to compute a matrix of Mutual Information.

"""

__authors__ = "Lilian Yang-crosson"
__license__ = "BSD 3-Clause License"
__version__ = "0.1.0"
__maintainer__ = "Lilian Yang-crosson"


from .pbtools import (mutual_information, mutual_information_matrix)

__all__ = ["mutual_information",
           "mutual_information_matrix"]
