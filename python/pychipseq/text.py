# -*- coding: utf-8 -*-
"""
Text parsing functions

Created on Sat Jan 31 16:46:42 2015

@author: antony
"""

from typing import List
import numpy as np

NA = "n/a"
EMPTY_ARRAY = np.array([])


def get_header(f):
    return f.readline().strip().split("\t")


def find_index(tokens: List[str], text: str, offset: int = 0):
    """
    Find the first heading in list that matches some text.
    """

    idx = find_indices(tokens, text, offset)

    if idx.size > 0:
        return idx[0]
    else:
        return -1


def find_indices(tokens: List[str], text: str, offset: int = 0) -> int:
    """
    Find all the headings in list that matches some text.
    """

    ltokens = [x.lower() for x in tokens]
    lt = text.lower()

    idx = np.where([lt in x for x in ltokens])[0]

    if idx.size > 0:
        return idx
    else:
        return EMPTY_ARRAY


def empty_line(n: int) -> str:
    """
    Produce an empty line of n tab separated n/as
    Args:
        n:  number of columns
    Return
        Tab separated string of n n/as
    """

    return '\t'.join([NA] * n)
