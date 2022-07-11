#!/bin/python
# -*- coding: utf-8 -*-

"""
Algorithmic Bioinfromatics 1: Maximum scoring subsequence algorithms.
Covered in exercise session 1, exercise 3

Created: May 3 08:54
@author: mpjw
"""
import random
from typing import List

import numpy as np


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


def MSS_recursive(a: List[int], n: int):
    """
    Recursive algorithm for finding minimal Scoring subsequences
    :param a: Array of integers
    :param n: Length of array a
    :return: Maximum scoring subsequence as Triple (max_score, left, right)
    """

    def sigma(a: List[int], i: int, j: int):
        if i > j:
            return 0
        elif i == j:
            return a[i]
        else:
            k = random.randint(i, j - 1)
            return sigma(a, i, k) + sigma(a, k + 1, j)  # for some k in [i:j-1]

    maxscore = 0
    l = 1
    r = 0

    for i in range(n):
        for j in range(i, n):
            s = sigma(a, i, j)  # compute s = sigma(a, i, j)
            if s > maxscore:
                maxscore = s
                l = i
                r = j

    return maxscore, l + 1, r


def MSS_dp(a: List[int], n: int):
    """
    Dynamic Programming algorithm to find MSS.
    :param a:
    :param n:
    :return:
    """

    maxscore = 0
    l = 1
    r = 0

    # dynamic programming matrix
    S = np.zeros(shape=(n, n), dtype=int)

    for i in range(n):
        for j in range(n):
            if i == j:
                S[i, i] = a[i]
            else:
                S[i, j] = S[i, j - 1] + a[j]
            if S[i, j] > maxscore:
                maxscore = S[i, j]
                l = i
                r = j


def MSS_dc(a: List[int], i: int, j: int):
    """
    Divide and conquer algorithm for finding MSS
    :param a:
    :param n:
    :return:
    """

    if i == j:
        if a[i] > 0:
            return (a[i], i, i)
        else:
            return 0, i, i - 1
    else:
        m = np.floor((i + j - 1) / 2)
        s_1, i_1, j_1 = MSS_dc(a, i, m)
        s_2, i_2, j_2 = MSS_dc(a, m + 1, j)
        i_3 = m
        s = a[i_3]
        simax = s
        for k in range(i_3 - 1, i, -1):
            s = s + a[k]
            if (s > simax):
                simax = s
                i_3 = k
        j_3 = m + 1
        s = a[j_3]
        sjmax = s
        for k in range(j_3 + 1, j):
            s = s + a[k]
            if (s > sjmax):
                sjmax = s
                j_3 = k
        s_3 = simax + sjmax  # s_3 = sigma(i_3, j_3)
        if np.max([s_1, s_2, s_3]) == s_1:
            return s_1, i_1, j_1
        elif np.max([s_1, s_2, s_3]) == s_2:
            return s_2, i_2, j_2
        else:
            return s_3, i_3, j_3


def MSS_clever(a: List[int], n: int):
    """

    :param a:
    :param n:
    :return:
    """

    maxscore = 0
    l = 1
    r = 0
    rmaxscore = 0
    rstart = 1

    for i in range(n):
        if (rmaxscore + a[i] > a[i]):
            rmaxscore = rmaxscore + a[i]
        else:
            rmaxscore = a[i]
            rstart = i
        if (rmaxscore > maxscore):
            maxscore = rmaxscore
            l = rstart
            r = i

    return maxscore, l, r


if __name__ == '__main__':
    # non naive algorithm: https://www.aaai.org/Papers/ISMB/1999/ISMB99-027.pdf
    test_MSS_Naive()
    test_MSS_Short_Naive()
    test_MSS_Short_Inverse_Naive()
