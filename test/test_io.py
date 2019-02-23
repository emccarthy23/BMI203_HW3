from BMI203_HW3 import io
import pandas as pd
import numpy as np
import os
from random import uniform
import time
import random

#Not testing write_alignment
def test_read_pairs():
    pairs = io.read_pairs('data/pairs/Pospairs.txt')
    length = 50
    assert len(pairs) == length

def test_read_sequence():
    output = io.read_sequence('data/sequences/prot-0004.fa')

    assert output[1] == 'SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM'


def test_read_scoring_matrix():
    #Check if the matrix is symmetric
    output = io.read_scoring_matrix('data/scoring/BLOSUM50')
    for i in range(len(output)):
        for j in range(i,len(output)):
            assert output.loc[i,j] == output.loc[j,i]
def test_write_optimal_matrix():
    input_matrix = io.read_scoring_matrix('data/scoring/BLOSUM50')
    output = io.write_optimal_matrix('test_write_optimal_matrix_file', 0, 'data/scoring/BLOSUM50', input_matrix, 0, 0,0)
    output_matrix = io.read_scoring_matrix('test_write_optimal_matrix_file')
    assert input_matrix == output_matrix
