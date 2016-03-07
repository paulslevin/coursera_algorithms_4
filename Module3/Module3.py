#! /usr/bin/python


def test_build_scoring_matrix(alphabet, diag_score,
                               off_diag_score, dash_score):
    pass

if __name__ == "__main__":
    #  unit tests
    import unittest

    class TestProjectFunctions(unittest.TestCase):
        def test_build_scoring_matrix_simple(self):
            self.assertEqual(test_build_scoring_matrix(alphabet=set(),
                                                       diag_score=0,
                                                       off_diag_score=0,
                                                       dash_score=0),
                             {"-": {"-": 0}})
    
    unittest.main()