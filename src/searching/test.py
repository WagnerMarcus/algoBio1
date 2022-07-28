#!/usr/bin/python

"""
Unit testing for text searching algorithms from Algorithmic Bioinformatics 1.
"""
import os
import sys
import time
import unittest

import numpy as np

from searching.KMP import naive_text_search, knuth_morris_pratt, compute_borders


class TestKMPMethods(unittest.TestCase):
    text_match_cases = ['ababaabb',
                        'abbbabbaababaabbababbabbaabbbabbaababaabbababbabbaabbbabbaababaabbababbabbaabbbabbaababaab',
                        'bbbbbbbababaaabbbabbaababaabbababbabbababaabbababbabbaabbbabbaababaabbababbabbaabbbabbaababba',
                        'aaaababaabbaaab' * 42,
                        'ababaabbabbbabbaababbbabbaababaabbababbabbaabbbabbaababaabbabaabbbabbaababaabbabbbabbaab' * 42,
                        'ababababaabbabbbabbaababaabbaabbbabbaababaabbaabbbabbaababaabbaabbbabbaababaabbaabbbab' * 3,
                        'ababaabbababaabbabbabbbabbaababaabbabaabbbabbaababaabbaabbbabbaababaabbaabbbabbaababaabbaabbb']
    text_no_match_cases = ['',
                           'a',
                           'aaaaaaab',
                           'a' * 1000,
                           'b' * 1000,
                           'aabbaabbaaaabbbb' * 42,
                           'ababaaba' * 42,
                           'ababaabababaab' * 42]
    word_to_find = 'ababaabb'

    def test_naive(self):
        for t in self.text_no_match_cases:
            self.assertFalse(naive_text_search(t, self.word_to_find))

        for t in self.text_match_cases:
            self.assertTrue(naive_text_search(t, self.word_to_find))

    def test_KMP(self, border_table=None):
        for t in self.text_no_match_cases:
            self.assertFalse(knuth_morris_pratt(t, self.word_to_find, border_table=border_table))

        for t in self.text_match_cases:
            self.assertTrue(knuth_morris_pratt(t, self.word_to_find, border_table=border_table))

    def test_runtime(self, k=1000):
        start = time.time()
        for i in range(k):
            self.test_naive()
        print('Runtime for {} naive tests: {} sec'.format(k, (time.time() - start)))

        start = time.time()
        # compute border table once
        border_table = np.zeros(len(self.word_to_find) + 1, dtype=int)
        compute_borders(border_table, self.word_to_find)

        for i in range(k):
            # reuse border table each call
            self.test_KMP(border_table=border_table)
        print('Runtime for {} KMP tests: {} sec'.format(k, (time.time() - start)))


# if __name__ == '__main__':
#     sys.path.insert(0, os.getcwd())
#     unittest.main()



