# -*- coding: utf-8 -*-

import itertools

import numpy as np


class SmithWatermanAlgorithm(object):
    _GAP_BASE = "-"

    def __init__(self, qry: str, ref: str, match_score: int = 3, gap_cost: int = 2):
        self._qry = qry
        self._ref = ref
        self._match_score = match_score
        self._gap_cost = gap_cost
        self._n_match = 0
        self._n_mismatch = 0
        self._n_gap = 0
        self._H = np.zeros((len(self._qry) + 1, len(self._ref) + 1), np.int)
        self._construct()
        self._traversal()

    def get_scoring_matrix(self):
        return self._H

    def get_qry_alignment(self):
        return "".join(self._qry_bases)

    def get_ref_alignment(self):
        return "".join(self._ref_bases)

    def get_qry_begin_index(self):
        return self._begin_i

    def get_ref_begin_index(self):
        return self._begin_j

    def get_qry_final_index(self):
        return self._final_i

    def get_ref_final_index(self):
        return self._final_j

    def _construct(self):
        for i, j in itertools.product(range(1, self._H.shape[0]), range(1, self._H.shape[1])):
            match = self._H[i - 1, j - 1] + (
                self._match_score if self._qry[i - 1] == self._ref[j - 1] else - self._match_score
            )
            insert = self._H[i - 1, j] - self._gap_cost
            delete = self._H[i, j - 1] - self._gap_cost
            self._H[i, j] = max(match, insert, delete, 0)

    def _traversal(self):
        self._qry_bases = list()
        self._ref_bases = list()
        self._begin_i, self._begin_j = 0, 0
        self._final_i, self._final_j = np.unravel_index(self._H.argmax(), self._H.shape)
        self._traceback(self._final_i, self._final_j)

    def _traceback(self, i: int, j: int):
        if self._H[i, j] == 0:
            self._begin_i = i
            self._begin_j = j
            return
        qry_base = self._qry[i - 1]
        ref_base = self._ref[j - 1]
        if self._H[i, j] == self._H[i - 1, j] - self._gap_cost:  # insert
            self._qry_bases.insert(0, qry_base)
            self._ref_bases.insert(0, self._GAP_BASE)
            self._n_gap += 1
            self._traceback(i - 1, j)
        elif self._H[i, j] == self._H[i, j - 1] - self._gap_cost:  # delete
            self._qry_bases.insert(0, self._GAP_BASE)
            self._ref_bases.insert(0, ref_base)
            self._n_gap += 1
            self._traceback(i, j - 1)
        elif self._H[i, j] in (self._H[i - 1, j - 1] + self._match_score, self._H[i - 1, j - 1] - self._match_score):
            self._qry_bases.insert(0, qry_base)
            self._ref_bases.insert(0, ref_base)
            if qry_base == ref_base:
                self._n_match += 1
            else:
                self._n_mismatch += 1
            self._traceback(i - 1, j - 1)
