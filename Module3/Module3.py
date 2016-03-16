#! /usr/bin/env python
"""Coursera Algorithms Project 4."""

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

    
#if __name__ == "__main__":
#    import unittest

#    class TestProjectFunctions(unittest.TestCase):
#        def test_build_scoring_matrix_simple(self):
#            outcome = build_scoring_matrix(alphabet=set(), diag_score=0,
#                                                off_diag_score=0, dash_score=0)
#            expected = {"-": {"-": 0}}
#            self.assertEqual(outcome, expected)
        
#        def test_build_scoring_matrix_one_letter(self):
#            outcome = build_scoring_matrix(alphabet={"A"}, 
#                                           diag_score=2,
#                                           off_diag_score=1,
#                                           dash_score=-1)
#            expected = {"-": {"-": -1, "A": -1}, "A": {"-": -1, "A": 2}}
#            self.assertEqual(outcome, expected)

#        def test_build_scoring_matrix_two_letters(self):
#            outcome = build_scoring_matrix(alphabet={"A", "B"},
#                                           diag_score=2,
#                                           off_diag_score=1,
#                                           dash_score=-1)
#            expected = {"-": {"-": -1, "A": -1, "B": -1},
#                        "A": {"-": -1, "A": 2, "B": 1},
#                        "B": {"-": -1, "A": 1, "B": 2}}
#            self.assertEqual(outcome, expected)

#        def test_compute_alignment_matrix_simple(self):
#            outcome = compute_alignment_matrix(seq_x="", seq_y="",
#                scoring_matrix={"-": {"-": -1, "A": -1},
#                                "A": {"-": -1, "A": 2}}, global_flag=True)
#            expected = [[0]]
#            self.assertEqual(outcome, expected)

#        def test_compute_alignment_matrix_one_letters(self):
#            outcome = compute_alignment_matrix(seq_x="A", seq_y="-",
#                scoring_matrix={"-": {"-": -1, "A": -1},
#                                "A": {"-": -1, "A": 2}}, global_flag=True)
#            expected = [[0, -1], [-1,  -1]]
#            self.assertEqual(outcome, expected)

#        def test_compute_alignment_matrix_two_letters(self):
#            outcome = compute_alignment_matrix(seq_x="AA", seq_y="TAAT",
#                scoring_matrix={"-": {"-": -6, "A": -6, "T": -6},
#                                "A": {"-": -6, "A": 10, "T": 4},
#                                "T": {"-": -6, "A": 4, "T": 10}},
#                global_flag=True)
#            expected = [[0, -6, -12, -18, -24],
#                        [-6, 4, 4, -2, -8],
#                        [-12, -2, 14, 14, 8]]
#            self.assertEqual(outcome, expected)
            
#        def test_compute_alignment_matrix_two_letters_local(self):
#            outcome = compute_alignment_matrix(seq_x="AA", seq_y="TAAT",
#                scoring_matrix={"-": {"-": -6, "A": -6, "T": -6},
#                                "A": {"-": -6, "A": 10, "T": 4},
#                                "T": {"-": -6, "A": 4, "T": 10}},
#                global_flag=False)
#            expected = [[0, 0, 0, 0, 0], [0, 4, 10, 10, 4], [0, 4, 14, 20, 14]]
#            self.assertEqual(outcome, expected)
            
#        def test_compute_global_alignment_two_letters(self):
#            outcome = compute_global_alignment(seq_x="AA", seq_y="TAAT",
#                scoring_matrix={"-": {"-": -6, "A": -6, "T": -6},
#                                "A": {"-": -6, "A": 10, "T": 4},
#                                "T": {"-": -6, "A": 4, "T": 10}},
#                alignment_matrix=[[0, -6, -12, -18, -24],
#                                  [-6, 4, 4, -2, -8],
#                                  [-12, -2, 14, 14, 8]])
#            expected = 8, "-AA-", "TAAT"
#            self.assertEqual(outcome, expected)

#        def test_compute_local_alignment_two_letters(self):
#            outcome = compute_local_alignment(seq_x="AA", seq_y="TAAT",
#                scoring_matrix={"-": {"-": -6, "A": -6, "T": -6},
#                                "A": {"-": -6, "A": 10, "T": 4},
#                                "T": {"-": -6, "A": 4, "T": 10}},
#                alignment_matrix=[[0, 0, 0, 0, 0],
#                                  [0, 4, 10, 10, 4],
#                                  [0, 4, 14, 20, 14]])
#            expected = 20, "AA", "AA"
#            self.assertEqual(outcome, expected)
                

#    unittest.main()
