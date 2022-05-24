# -*- coding: utf-8 -*-

import unittest

from algorithm import SmithWatermanAlgorithm


class SmithWatermanAlgorithmTestCase(unittest.TestCase):
    def test_alpha(self):
        qry = "CTTGACGTGTTTATGTATTCTTTTGCCAGTATATATTCTACACACCATATTATCTGCTGCAACCAAAAGACACAATGTTC"
        ref = "CCGCTTTTAAGGGCTATATCCGTCCCTAGACCAATATAATAGTTCGTCTATGTGATCTCTTGAATTACGCATTCTATTGG"
        smith_waterman_algorithm = SmithWatermanAlgorithm(qry, ref)
        self.assertEqual(smith_waterman_algorithm.get_qry_begin_index(), 5)
        self.assertEqual(smith_waterman_algorithm.get_qry_final_index(), 65)
        self.assertEqual(smith_waterman_algorithm.get_ref_begin_index(), 1)
        self.assertEqual(smith_waterman_algorithm.get_ref_final_index(), 71)

    def test_beta(self):
        qry = "GTAGCTAGAGGTGAGACCCCCGTAAACACCAGCAATGGCAGGATTAAGAGAAGTAGAAGTAGGGGCCGGAGATCCGTCCT"
        ref = "AACATTCTAGATAGTTACTGCAGCGCCCATTGTTCTGGACTAGTGCCTTGTGTGAGAATTCGGAGGTTCCGGGCCAAATC"
        smith_waterman_algorithm = SmithWatermanAlgorithm(qry, ref)
        self.assertEqual(smith_waterman_algorithm.get_qry_begin_index(), 4)
        self.assertEqual(smith_waterman_algorithm.get_qry_final_index(), 74)
        self.assertEqual(smith_waterman_algorithm.get_ref_begin_index(), 6)
        self.assertEqual(smith_waterman_algorithm.get_ref_final_index(), 80)

    def test_gamma(self):
        qry = "ATAAATGATGGGAACGAGATCCCGGAGGCTCGGATTGGTATGACAAGGTGTATCGTGATCGTCGGTGCGTCAGCTTGGGC"
        ref = "GGTAAGGTATAGCTGCATCCTACTTACGATGTGAAGTTACACACCTCAACTCCAGAGTCCCGTTGGGGGAGTGTATTTTT"
        smith_waterman_algorithm = SmithWatermanAlgorithm(qry, ref)
        qry_begin_index = smith_waterman_algorithm.get_qry_begin_index()
        qry_final_index = smith_waterman_algorithm.get_qry_final_index()
        ref_begin_index = smith_waterman_algorithm.get_ref_begin_index()
        ref_final_index = smith_waterman_algorithm.get_ref_final_index()
        self.assertEqual(smith_waterman_algorithm.get_qry_begin_index(), 10)
        self.assertEqual(smith_waterman_algorithm.get_qry_final_index(), 76)
        self.assertEqual(smith_waterman_algorithm.get_ref_begin_index(), 5)
        self.assertEqual(smith_waterman_algorithm.get_ref_final_index(), 74)

        qry_alignment = smith_waterman_algorithm.get_qry_alignment()
        ref_alignment = smith_waterman_algorithm.get_ref_alignment()
        qry_matching = qry[:qry_begin_index] + " (" + qry_alignment + ") " + qry[qry_final_index:]
        ref_matching = ref[:ref_begin_index] + " (" + ref_alignment + ") " + ref[ref_final_index:]

        print("\nHere are the details of sequence matching:")
        # ATAAATGATG (GG-A-A-C-GAGATCCCGGAGGCT--CGGAT-TG--GTATGACA-AGGTGTA-TCGTGA-TC--GTCGGTGCGTCAGCT-T) GGGC
        #      GGTAA (GGTATAGCTGC-ATCCT--A--CTTACG-ATGTGAAGT-T-ACACACCTCAACTCCAGAGTCCCGTTGGGG-G--AG-TGT) ATTTTT
        if qry_begin_index == ref_begin_index:
            print(qry_matching)
            print(ref_matching)
        elif qry_begin_index < ref_begin_index:
            print(" " * (ref_begin_index - qry_begin_index) + qry_matching)
            print(ref_matching)
        elif qry_begin_index > ref_begin_index:
            print(qry_matching)
            print(" " * (qry_begin_index - ref_begin_index) + ref_matching)

    def test_delta(self):
        qry = "GTGGCAACATCTCACAATTGCCAGTTAACGTCTTCCTTCTCTCTCTGTCATAGGGACTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCA" \
              "AGGAATTAAGAGAAGCAACATCTCCGGAAGCCAACAAGGAAATCCTCGATGTGAG"
        ref = "GTGGCAcCATCTCACAATTGCCAGTTAACGTCTTCCTTCTCTCTCTGTCATAGGGACTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCA" \
              "AGGAATTAAGAGAAGCAACATCTCCGaAAGCCAACAAGGAAATCCTCGATGTGAG"
        smith_waterman_algorithm = SmithWatermanAlgorithm(qry, ref)
        qry_begin_index = smith_waterman_algorithm.get_qry_begin_index()
        qry_final_index = smith_waterman_algorithm.get_qry_final_index()
        ref_begin_index = smith_waterman_algorithm.get_ref_begin_index()
        ref_final_index = smith_waterman_algorithm.get_ref_final_index()

        qry_alignment = smith_waterman_algorithm.get_qry_alignment()
        ref_alignment = smith_waterman_algorithm.get_ref_alignment()
        qry_matching = qry[:qry_begin_index] + " (" + qry_alignment + ") " + qry[qry_final_index:]
        ref_matching = ref[:ref_begin_index] + " (" + ref_alignment + ") " + ref[ref_final_index:]

        print("\nHere are the details of sequence matching:")
        if qry_begin_index == ref_begin_index:
            print(qry_matching)
            print(ref_matching)
        elif qry_begin_index < ref_begin_index:
            print(" " * (ref_begin_index - qry_begin_index) + qry_matching)
            print(ref_matching)
        elif qry_begin_index > ref_begin_index:
            print(qry_matching)
            print(" " * (qry_begin_index - ref_begin_index) + ref_matching)

        self.assertEqual(qry_begin_index, 0)
        self.assertEqual(qry_final_index, 157)
        self.assertEqual(ref_begin_index, 0)
        self.assertEqual(ref_final_index, 157)
