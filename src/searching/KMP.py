#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Implementations for the Knuth-Morris-Pratt (KMP) algorithm used for searching a word in a text.
These implementations are derived from the lecture 'Algorithmic Bioinformatics I' at LMU Munich summer term 2022

Created: July 28 15:50
@author: mpjw
"""

import numpy as np


def naive_text_search(text: str, word: str) -> bool:
    """


    See Naive pseudo code in https://www.bio.ifi.lmu.de/mitarbeiter/volker-heun/notes/ab8.pdf#subsection.2.2.1
    :param text: text string to search in (char t[] in script)
    :param word: word to search in text as string (char s[] in script)

    The parameter int n and int m from the script correspond to text and word length.
    """

    n = len(text)
    m = len(word)

    i = j = 0

    while i <= n - m:
        while text[i + j] == word[j]:
            j += 1
            if j == m:
                return True
        i += 1
        j = 0
    return False


def compute_borders(border_table: np.ndarray, word: str):
    """
    In place calculation of border table.
    See https://www.bio.ifi.lmu.de/mitarbeiter/volker-heun/notes/ab8.pdf#subsection.2.2.6

    :param border_table: Numpy array to be filled with border lengths (int[] border in script)
    :param word: word to be searched in text (char[] s in script)
    Parameter int m in script corresponds to length of s aka word.

    :return:
    """

    border_table[0] = -1
    border_table[1] = 0

    i = border_table[1]
    for j in range(2, len(word)):
        # Note that we have: i = border_table[j - 1]

        while i >= 0 and word[i] != word[j - 1]:
            i = border_table[i]

        i += 1
        border_table[j] = i


def knuth_morris_pratt(text: str, word: str, border_table=None) -> bool:
    """

    See https://www.bio.ifi.lmu.de/mitarbeiter/volker-heun/notes/ab8.pdf#subsection.2.2.4
    :param text: text string to search in (char t[] in script)
    :param word: word to search in text as string (char s[] in script)
    :param border_table: precomputed border table

    The parameter int n and int m from the script correspond to text and word length.
    """

    # initialize and fill border table
    if border_table is None:
        border_table = np.zeros(len(word) + 1, dtype=int)
        compute_borders(border_table=border_table, word=word)

    n = len(text)
    m = len(word)

    i = j = 0
    while i <= n - m:
        while text[i + j] == word[j]:
            j += 1
            if j == m:
                return True

        i = i + (j - border_table[j])  # It holds that j − border[j] > 0
        j = max(0, border_table[j])

    return False
