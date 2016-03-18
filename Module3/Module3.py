#! /usr/bin/env python
"""Coursera Algorithms Project 4."""

import math
import random
import urllib2
import matplotlib.pyplot as plt

# Project
def set_matrix_value(matrix, key, letter, score1, score2, score3):
    """Helper function for build_scoring_matrix."""
    if key == "-" or letter == "-":
        matrix[key][letter] = score1
    elif key != letter:
        matrix[key][letter] = score2
    else:
        matrix[key][letter] = score3


def build_scoring_matrix(alphabet, diag_score, 
                              off_diag_score, dash_score):
    """Builds scoring matrix."""
    augmented = alphabet | {"-"}
    matrix = {letter: {} for letter in augmented}
    for key in matrix:
        for letter in augmented:
            set_matrix_value(matrix, key, letter,
                             dash_score, off_diag_score, diag_score)
    return matrix


def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """Computes alignment matrix."""
    matrix = [[0 for __ in seq_y] + [0] for _ in seq_x]
    matrix.append([0 for __ in seq_y] + [0])
    for idx, val in enumerate(seq_x):
        value = matrix[idx][0] + scoring_matrix[val]["-"]
        if global_flag:
            matrix[idx + 1][0] = value
        else:
            matrix[idx + 1][0] = max(0, value)
    for idx, val in enumerate(seq_y):
        value = matrix[0][idx] + scoring_matrix["-"][val]
        if global_flag:
            matrix[0][idx + 1] = value
        else:
            matrix[0][idx + 1] = max(0, value)
    for id1, val1 in enumerate(seq_x):
        for id2, val2 in enumerate(seq_y):
            value = max(matrix[id1][id2] + scoring_matrix[val1][val2],
                         matrix[id1][id2 + 1] + scoring_matrix[val1]["-"],
                         matrix[id1 + 1][id2] + scoring_matrix["-"][val2])
            if global_flag:
                matrix[id1 + 1][id2 + 1] = value
            else:
                matrix[id1 + 1][id2 + 1] = max(value, 0)
    return matrix

        
def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """Computes global alignment matrix."""
    idx, idy = len(seq_x), len(seq_y)
    score = alignment_matrix[idx][idy]
    seq_x_prime, seq_y_prime = "", ""
    while idx and idy:
        if alignment_matrix[idx][idy] == alignment_matrix[idx - 1][idy - 1] + scoring_matrix[seq_x[idx - 1]][seq_y[idy - 1]]:
            seq_x_prime = seq_x[idx - 1] + seq_x_prime
            seq_y_prime = seq_y[idy - 1] + seq_y_prime
            idx -= 1
            idy -= 1
        else:
            if alignment_matrix[idx][idy] == alignment_matrix[idx - 1][idy] + scoring_matrix[seq_x[idx - 1]]["-"]:
                seq_x_prime = seq_x[idx - 1] + seq_x_prime
                seq_y_prime = "-" + seq_y_prime
                idx -= 1
            else:
                seq_x_prime = "-" + seq_x_prime
                seq_y_prime = seq_y[idy - 1] + seq_y_prime
                idy -= 1
    while idx:
        seq_x_prime = seq_x[idx - 1] + seq_x_prime
        seq_y_prime = "-" + seq_y_prime
        idx -= 1
    while idy:
        seq_x_prime = "-" + seq_x_prime
        seq_y_prime = seq_y[idy - 1] + seq_y_prime
        idy -= 1
    return score, seq_x_prime, seq_y_prime


def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """Computes local alignment matrix."""
    max_value = max(sum(alignment_matrix, []))
    for idx, row in enumerate(alignment_matrix):
        if max_value in row:
            start_row = idx
            start_col = row.index(max_value)
            break
    seq_x_prime, seq_y_prime = "", ""
    score = alignment_matrix[start_row][start_col]
    idx, idy = start_row, start_col
    while idx and idy:
        if not alignment_matrix[idx][idy]:
            return score, seq_x_prime, seq_y_prime
        if alignment_matrix[idx][idy] == alignment_matrix[idx - 1][idy - 1] + \
                scoring_matrix[seq_x[idx - 1]][seq_y[idy - 1]]:
            seq_x_prime = seq_x[idx - 1] + seq_x_prime
            seq_y_prime = seq_y[idy - 1] + seq_y_prime
            idx -= 1
            idy -= 1
        else:
            if alignment_matrix[idx][idy] == alignment_matrix[idx - 1][idy] + \
                    scoring_matrix[seq_x[idx - 1]]["-"]:
                seq_x_prime = seq_x[idx - 1] + seq_x_prime
                seq_y_prime = "-" + seq_y_prime
                idx -= 1
            else:
                seq_x_prime = "-" + seq_x_prime
                seq_y_prime = seq_y[idy - 1] + seq_y_prime
                idy -= 1
    while idx:
        if not alignment_matrix[idx][idy]:
            return score, seq_x_prime, seq_y_prime
        seq_x_prime = seq_x[idx - 1] + seq_x_prime
        seq_y_prime = "-" + seq_y_prime
        idx -= 1
    while idy:
        if not alignment_matrix[idx][idy]:
            return score, seq_x_prime, seq_y_prime
        seq_x_prime = "-" + seq_x_prime
        seq_y_prime = seq_y[idy - 1] + seq_y_prime
        idy -= 1
    return score, seq_x_prime, seq_y_prime

# Application
# URLs for data files
PAM50_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_PAM50.txt"
HUMAN_EYELESS_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_HumanEyelessProtein.txt"
FRUITFLY_EYELESS_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_FruitflyEyelessProtein.txt"
CONSENSUS_PAX_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_ConsensusPAXDomain.txt"
WORD_LIST_URL = "http://storage.googleapis.com/codeskulptor-assets/assets_scrabble_words3.txt"


def read_scoring_matrix(filename):
    """
    Read a scoring matrix from the file named filename.

    Argument:
    filename -- name of file containing a scoring matrix

    Returns:
    A dictionary of dictionaries mapping X and Y characters to scores
    """
    scoring_dict = {}
    scoring_file = urllib2.urlopen(filename)
    ykeys = scoring_file.readline()
    ykeychars = ykeys.split()
    for line in scoring_file.readlines():
        vals = line.split()
        xkey = vals.pop(0)
        scoring_dict[xkey] = {}
        for ykey, val in zip(ykeychars, vals):
            scoring_dict[xkey][ykey] = int(val)
    return scoring_dict




def read_protein(filename):
    """
    Read a protein sequence from the file named filename.

    Arguments:
    filename -- name of file containing a protein sequence

    Returns:
    A string representing the protein
    """
    protein_file = urllib2.urlopen(filename)
    protein_seq = protein_file.read()
    protein_seq = protein_seq.rstrip()
    return protein_seq


def read_words(filename):
    """
    Load word list from the file named filename.

    Returns a list of strings.
    """
    # load assets
    word_file = urllib2.urlopen(filename)

    # read in files as string
    words = word_file.read()

    # template lines and solution lines list of line string
    word_list = words.split('\n')
    print "Loaded a dictionary with", len(word_list), "words"
    return word_list


human_protein = read_protein(HUMAN_EYELESS_URL)
fly_protein = read_protein(FRUITFLY_EYELESS_URL)
scoring_matrix = read_scoring_matrix(PAM50_URL)

global_alignment_matrix = compute_global_alignment(seq_x=human_protein,
                                                   seq_y=fly_protein,
                                                   scoring_matrix=scoring_matrix,)


print compute_local_alignment(seq_x=human_protein,
                              seq_y=fly_protein,
                              alignment_matrix=)


if __name__ == "__main__":
   import unittest

   class TestProjectFunctions(unittest.TestCase):
       def test_build_scoring_matrix_simple(self):
           outcome = build_scoring_matrix(alphabet=set(), diag_score=0,
                                               off_diag_score=0, dash_score=0)
           expected = {"-": {"-": 0}}
           self.assertEqual(outcome, expected)

       def test_build_scoring_matrix_one_letter(self):
           outcome = build_scoring_matrix(alphabet={"A"},
                                          diag_score=2,
                                          off_diag_score=1,
                                          dash_score=-1)
           expected = {"-": {"-": -1, "A": -1}, "A": {"-": -1, "A": 2}}
           self.assertEqual(outcome, expected)

       def test_build_scoring_matrix_two_letters(self):
           outcome = build_scoring_matrix(alphabet={"A", "B"},
                                          diag_score=2,
                                          off_diag_score=1,
                                          dash_score=-1)
           expected = {"-": {"-": -1, "A": -1, "B": -1},
                       "A": {"-": -1, "A": 2, "B": 1},
                       "B": {"-": -1, "A": 1, "B": 2}}
           self.assertEqual(outcome, expected)

       def test_compute_alignment_matrix_simple(self):
           outcome = compute_alignment_matrix(seq_x="", seq_y="",
               scoring_matrix={"-": {"-": -1, "A": -1},
                               "A": {"-": -1, "A": 2}}, global_flag=True)
           expected = [[0]]
           self.assertEqual(outcome, expected)

       def test_compute_alignment_matrix_one_letters(self):
           outcome = compute_alignment_matrix(seq_x="A", seq_y="-",
               scoring_matrix={"-": {"-": -1, "A": -1},
                               "A": {"-": -1, "A": 2}}, global_flag=True)
           expected = [[0, -1], [-1,  -1]]
           self.assertEqual(outcome, expected)

       def test_compute_alignment_matrix_two_letters(self):
           outcome = compute_alignment_matrix(seq_x="AA", seq_y="TAAT",
               scoring_matrix={"-": {"-": -6, "A": -6, "T": -6},
                               "A": {"-": -6, "A": 10, "T": 4},
                               "T": {"-": -6, "A": 4, "T": 10}},
               global_flag=True)
           expected = [[0, -6, -12, -18, -24],
                       [-6, 4, 4, -2, -8],
                       [-12, -2, 14, 14, 8]]
           self.assertEqual(outcome, expected)

       def test_compute_alignment_matrix_two_letters_local(self):
           outcome = compute_alignment_matrix(seq_x="AA", seq_y="TAAT",
               scoring_matrix={"-": {"-": -6, "A": -6, "T": -6},
                               "A": {"-": -6, "A": 10, "T": 4},
                               "T": {"-": -6, "A": 4, "T": 10}},
               global_flag=False)
           expected = [[0, 0, 0, 0, 0], [0, 4, 10, 10, 4], [0, 4, 14, 20, 14]]
           self.assertEqual(outcome, expected)

       def test_compute_global_alignment_two_letters(self):
           outcome = compute_global_alignment(seq_x="AA", seq_y="TAAT",
               scoring_matrix={"-": {"-": -6, "A": -6, "T": -6},
                               "A": {"-": -6, "A": 10, "T": 4},
                               "T": {"-": -6, "A": 4, "T": 10}},
               alignment_matrix=[[0, -6, -12, -18, -24],
                                 [-6, 4, 4, -2, -8],
                                 [-12, -2, 14, 14, 8]])
           expected = 8, "-AA-", "TAAT"
           self.assertEqual(outcome, expected)

       def test_compute_local_alignment_two_letters(self):
           outcome = compute_local_alignment(seq_x="AA", seq_y="TAAT",
               scoring_matrix={"-": {"-": -6, "A": -6, "T": -6},
                               "A": {"-": -6, "A": 10, "T": 4},
                               "T": {"-": -6, "A": 4, "T": 10}},
               alignment_matrix=[[0, 0, 0, 0, 0],
                                 [0, 4, 10, 10, 4],
                                 [0, 4, 14, 20, 14]])
           expected = 20, "AA", "AA"
           self.assertEqual(outcome, expected)


   unittest.main()
