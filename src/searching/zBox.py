#!/usr/bin/python

"""
Z-Box algorithm implementations for Knuth-Morris-Pratt and Boyer-Moore algorithms.
Algorithmic Bioinformatics 1

If you want to go through an algorithm step by step the <algorithm_name>_verbose() function print a lot of details for
 that ;)

@author: mpjw
"""

import numpy as np


def zBoxes(s: str, m: int):
    """
    Z-Boxes algorithm found in script chapter 2.5. "Z-Box-Algorithmen"

    :param s: input word as string
    :param m: length of inout word
    :return: Z-Box array
    """

    l = r = 0  # *no Z - box seen so far
    i = 1  # Note that i will never be decremented
    Z = np.zeros(m)

    for k in range(1, m):

        if k > r:                               # case 1
            i = k
            while i < m and s[i] == s[i-k]:
                i += 1
            Z[k] = i - k
            if Z[k] > 0:
                l = k
                r = i - 1
        else:
            if Z[k-l] < r-k+1:                  # case 2a
                Z[k] = Z[k-l]
            else:                               # case 2b
                i = r + 1
                while i < m and s[i] == s[i-k]:
                    i += 1
                Z[k] = i-k
                if i-1 > r:
                    l = k
                    r = i-1

    return Z


def zBoxes_verbose(s: str, m: int):
    """
    Z-Boxes algorithm found in script chapter 2.5. "Z-Box-Algorithmen"
    Detailed verbose version giving variable values each iteration.

    :param s: input word as string
    :param m: length of inout word
    :return: Z-Box array
    """

    l = r = 0  # *no Z - box seen so far
    i = 1  # Note that i will never be decremented
    Z = np.zeros(m)

    for k in range(1, m):
        print('iteration k=', str(k), ':', 'i=', str(i), 'l=', str(l), 'r=', str(r))

        if k > r:                               # case 1
            print('case 1: k>r', str(k), '>', str(r))
            i = k
            while i < m and s[i] == s[i-k]:
                i += 1
            Z[k] = i - k
            if Z[k] > 0:
                l = k
                r = i - 1
        else:
            if Z[k-l] < r-k+1:                  # case 2a
                print('case 2a: k<=r', str(k), '<=', str(r), 'Z[k-l] < r-k+1:', str(Z[k-l]), '<', str(r-k+1))
                Z[k] = Z[k-l]
            else:                               # case 2b
                print('case 2b: k<=r', str(k), '<=', str(r), 'Z[k-l] >= r-k+1:', str(Z[k-l]), '<', str(r-k+1))
                i = r + 1
                while i < m and s[i] == s[i-k]:
                    i += 1
                Z[k] = i-k
                if i-1 > r:
                    l = k
                    r = i-1

    print('final Z-Boxes:', ', '.join(map(str, Z[1:])))
    return Z


if __name__ == '__main__':
    # run your test cases here
    zBoxes_verbose('aabaababaaa', len('aabaababaaa'))


