from BMI203_HW3 import alignment
from BMI203_HW3 import io
import os
import csv
import pandas as pd
import numpy as np
import os
from random import uniform
import time
import random

def test_smith_waterman_and_scoring_algorithms_round1():
    filename_a = "data/sequences/prot-0031.fa"
    filename_b = "data/sequences/prot-0034.fa"
    seq_a, seq_b, score_1 =  alignment.smith_waterman_alignment(filename_a, filename_b, 'data/scoring/PAM100', 4, 1)
    score_2 =  alignment.smith_waterman(filename_a, filename_b, 'data/scoring/PAM100', 4, 1)
    score_3 =  alignment.smith_waterman(filename_a, filename_b, 'data/scoring/PAM100', 4, 1)
    scoring_matrix = io.read_scoring_matrix('data/scoring/PAM100')
    score_4 =  alignment.score_alignment(seq_a, seq_b, scoring_matrix,4,1)
    # Test that the score equals what you would get from EMBOSS for these sequences.
    assert score_1 == 74
    assert score_2 == 74
    assert score_3 == 74
    assert score_4 == 74

def test_smith_waterman_and_scoring_algorithms_round2():
    filename_a = "data/sequences/prot-0102.fa"
    filename_b = "data/sequences/prot-0098.fa"
    seq_a, seq_b, score_1 =  alignment.smith_waterman_alignment(filename_a, filename_b, 'data/scoring/PAM100', 4, 1)
    score_2 =  alignment.smith_waterman(filename_a, filename_b, 'data/scoring/PAM100', 4, 1)
    score_3 =  alignment.smith_waterman(filename_a, filename_b, 'data/scoring/PAM100', 4, 1)
    scoring_matrix = io.read_scoring_matrix('data/scoring/PAM100')
    score_4 =  alignment.score_alignment(seq_a, seq_b, scoring_matrix,4,1)
    # Test that the score equals what you would get from EMBOSS for these sequences.
    assert score_1 == 38
    assert score_2 == 38
    assert score_3 == 38
    assert score_4 == 38

def test_smith_waterman_and_scoring_algorithms_empty_seq():
    filename_a = "data/tests/empty_seq.fa"
    filename_b = "data/tests/empty_seq.fa"
    seq_a, seq_b, score_1 =  alignment.smith_waterman_alignment(filename_a, filename_b, 'data/scoring/PAM100', 4, 1)
    score_2 =  alignment.smith_waterman(filename_a, filename_b, 'data/scoring/PAM100', 4, 1)
    score_3 =  alignment.smith_waterman(filename_a, filename_b, 'data/scoring/PAM100', 4, 1)
    scoring_matrix = io.read_scoring_matrix('data/scoring/PAM100')
    score_4 =  alignment.score_alignment(seq_a, seq_b, scoring_matrix,4,1)
    # Test that the score equals what you would get from EMBOSS for these sequences.
    assert score_1 == 0
    assert score_2 == 0
    assert score_3 == 0
    assert score_4 == 0

def test_optimize_scoring_matrix():
    #Test that a matrix produced by the function is symmetric
    output = io.read_scoring_matrix('data/optimization/opt_PAM100')

    for i in range(24):
        for j in range(i,24):
            assert output.iloc[i,j] == output.iloc[j,i]
