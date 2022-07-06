#!/bin/python
# -*- coding: utf-8 -*-

"""
Algorithmic Bioinfromatics 1: Maximum scoring subsequence algorithms.
Covered in exercise session 1, exercise 3

Created: May 3 08:54
@author: mpjw
"""
from typing import List


def MSS_Naive(a: List[int], n: int):
    """
    Naive maximum scoring subsequence algorithm from Script (p. 46).
    :param a: Array of integers (aka our sequence)
    :param n: length of Array a
    :return: A maximal scoring subsequence as interval (i, j)
    """

    # always is maxscore = σ(l, r)
    maxscore, l, r = 0, 0, 0

    for i in range(0, n):
        for j in range(i, n):

            # compute s = σ(i, j)
            s = 0
            for k in range(i, j):
                s = s + a[k]
            if s > maxscore:
                maxscore, l, r = s, i, j

    return maxscore, l + 1, r


def test_MSS_Naive(**kwargs):
    """
    Test the MSS_Naive algorithm with default test data.
    Prints solutions as (maxscore, l, r) with l, r denoting interval boundaries of the subsequence
    """

    # default test data
    a = kwargs.get("a") if "a" in kwargs.keys() else [1, -1, 3, -1, 3, -2, 2, -3, -2, 4, -2, 3, -2, 2]

    print("Test: MSS_Naive")
    print("Input:", a)
    print("Solution:", MSS_Naive(a, len(a)))


def MSS_Short_Naive(a: List[int], n: int):
    """
    Naive maximum scoring subsequence algorithm for shortest solution from Tutor-session 1.
    :param a: Array of integers (aka our sequence)
    :param n: length of Array a
    :return: A shortest maximal scoring subsequence as interval (i, j)
    """

    # always is maxscore = σ(l, r)
    maxscore, l, r = 0, 0, 0

    for i in range(0, n):
        for j in range(i, n):

            # compute s = σ(i, j)
            s = 0
            for k in range(i, j):
                s = s + a[k]
            if s > maxscore:
                maxscore, l, r = s, i, j

            # check subsequence length
            elif s == maxscore and j - i + 1 < r - l + 1:
                l, r = i, j

    return maxscore, l + 1, r


def test_MSS_Short_Naive(**kwargs):
    """
    Test the MSS_Short_Naive algorithm with default test data.
    Prints solutions as (maxscore, l, r) with l, r denoting interval boundaries of the subsequence
    """

    # default test data
    a = kwargs.get("a") if "a" in kwargs.keys() else [1, -1, 3, -1, 3, -2, 2, -3, -2, 4, -2, 3, -2, 2]

    print("Test: MSS_Short_Naive")
    print("Input:", a)
    print("Solution:", MSS_Short_Naive(a, len(a)))


def MSS_Short_Inverse_Naive(a: List[int], n: int):
    """
    Naive maximum scoring subsequence algorithm for shortest solution from Tutor-session 1.
    :param a: Array of integers (aka our sequence)
    :param n: length of Array a
    :return: A shortest maximal scoring subsequence as interval (i, j)
    """

    # always is maxscore = σ(l, r)
    maxscore, l, r = 0, 0, 0

    for i in range(n, 0, -1):
        for j in range(i, n):

            # compute s = σ(i, j)
            s = 0
            for k in range(i, j):
                s = s + a[k]
            if s > maxscore:
                maxscore, l, r = s, i, j

    return maxscore, l + 1, r


def test_MSS_Short_Inverse_Naive(**kwargs):
    """
    Test the MSS_Short_Inverse_Naive algorithm with default test data.
    Prints solutions as (maxscore, l, r) with l, r denoting interval boundaries of the subsequence
    """

    # default test data
    a = kwargs.get("a") if "a" in kwargs.keys() else [1, -1, 3, -1, 3, -2, 2, -3, -2, 4, -2, 3, -2, 2]

    print("Test: MSS_Short_Inverse_Naive")
    print("Input:", a)
    print("Solution:", MSS_Short_Inverse_Naive(a, len(a)))


# non naive algorithm: https://www.aaai.org/Papers/ISMB/1999/ISMB99-027.pdf


if __name__ == '__main__':
    test_MSS_Naive()
    test_MSS_Short_Naive()
    test_MSS_Short_Inverse_Naive()
