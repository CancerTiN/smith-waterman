# -*- coding: utf-8 -*-

from dataclasses import FrozenInstanceError
import unittest

import numpy as np

from algorithm import ScoringScheme, smith_waterman


class SmithWatermanTestCase(unittest.TestCase):
    def test_finds_best_local_alignment_for_simple_sequences(self):
        result = smith_waterman("AACCGGTT", "TTCCGGAA")

        self.assertEqual(result.query_alignment, "CCGG")
        self.assertEqual(result.reference_alignment, "CCGG")
        self.assertEqual(result.query_span, (2, 6))
        self.assertEqual(result.reference_span, (2, 6))
        self.assertEqual(result.score, 12)
        self.assertEqual(result.matches, 4)
        self.assertEqual(result.mismatches, 0)
        self.assertEqual(result.gaps, 0)

    def test_supports_gap_in_best_alignment(self):
        result = smith_waterman(
            "ACGT",
            "AGT",
            ScoringScheme(match=3, mismatch=-3, gap=-2),
        )

        self.assertEqual(result.query_alignment, "ACGT")
        self.assertEqual(result.reference_alignment, "A-GT")
        self.assertEqual(result.query_span, (0, 4))
        self.assertEqual(result.reference_span, (0, 3))
        self.assertEqual(result.score, 7)
        self.assertEqual(result.matches, 3)
        self.assertEqual(result.mismatches, 0)
        self.assertEqual(result.gaps, 1)

    def test_returns_empty_alignment_when_no_positive_score_exists(self):
        result = smith_waterman("AAAA", "TTTT")

        self.assertEqual(result.query_alignment, "")
        self.assertEqual(result.reference_alignment, "")
        self.assertEqual(result.query_span, (0, 0))
        self.assertEqual(result.reference_span, (0, 0))
        self.assertEqual(result.score, 0)
        self.assertEqual(result.matches, 0)
        self.assertEqual(result.mismatches, 0)
        self.assertEqual(result.gaps, 0)

    def test_handles_empty_sequences(self):
        result = smith_waterman("", "ACGT")

        self.assertEqual(result.query_alignment, "")
        self.assertEqual(result.reference_alignment, "")
        self.assertEqual(result.query_span, (0, 0))
        self.assertEqual(result.reference_span, (0, 0))
        self.assertEqual(result.score, 0)
        self.assertEqual(result.scoring_matrix.shape, (1, 5))

    def test_exposes_numpy_int64_scoring_matrix(self):
        result = smith_waterman("GATTACA", "GCATGCU")

        self.assertIsInstance(result.scoring_matrix, np.ndarray)
        self.assertEqual(result.scoring_matrix.dtype, np.dtype("int64"))
        self.assertEqual(result.scoring_matrix.shape, (8, 8))
        self.assertEqual(result.scoring_matrix.max(), result.score)

    def test_returns_read_only_result(self):
        result = smith_waterman("ACGT", "ACGT")

        with self.assertRaises(FrozenInstanceError):
            result.score = 0

        with self.assertRaises(ValueError):
            result.scoring_matrix[1, 1] = 99

    def test_rejects_invalid_scoring_scheme(self):
        with self.assertRaisesRegex(ValueError, "match must be positive"):
            ScoringScheme(match=0)

        with self.assertRaisesRegex(ValueError, "mismatch must be negative"):
            ScoringScheme(mismatch=0)

        with self.assertRaisesRegex(ValueError, "gap must be negative"):
            ScoringScheme(gap=0)

    def test_preserves_long_sequence_regression_bounds(self):
        query = "CTTGACGTGTTTATGTATTCTTTTGCCAGTATATATTCTACACACCATATTATCTGCTGCAACCAAAAGACACAATGTTC"
        reference = "CCGCTTTTAAGGGCTATATCCGTCCCTAGACCAATATAATAGTTCGTCTATGTGATCTCTTGAATTACGCATTCTATTGG"

        result = smith_waterman(query, reference)

        self.assertEqual(result.query_span, (5, 65))
        self.assertEqual(result.reference_span, (1, 71))
        self.assertEqual(len(result.query_alignment), len(result.reference_alignment))


if __name__ == "__main__":
    unittest.main()
