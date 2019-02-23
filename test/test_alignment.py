from BMI203_HW3 import alignment
from BMI203_HW3 import io
import os
import csv
smith_waterman_alignment, smith_waterman,smith_waterman_len_adj, score_alignment, score_performance, optimize_scoring_matrix
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
    filename_b = "data/tests/empty_seq.fa
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

def test_score_performance():
    #Make sure that it gives a value between 0 and 4
    with open('tests/pos_seq.csv', 'r') as f:
        reader = csv.reader(f)
        pos_alignments = list(reader)
        f.close()
    with open('tests/neg_seq.csv', 'r') as f:
        reader = csv.reader(f)
        neg_alignments = list(reader)
        f.close()
    scoring_matrix = io.read_scoring_matrix('data/scoring/PAM100')
    score = score_performance(pos_alignments,neg_alignments,scoring_matrix, gap,ext)[0]
    # Test that the score equals what you would get from EMBOSS for these sequences.
    assert score <= 4
    assert score >= 0

def test_optimize_scoring_matrix():
    #Make sure outputed matrix is symmetric
    with open('tests/pos_seq.csv', 'r') as f:
        reader = csv.reader(f)
        pos_alignments = list(reader)
        f.close()
    with open('tests/neg_seq.csv', 'r') as f:
        reader = csv.reader(f)
        neg_alignments = list(reader)
        f.close()
    output = optimize_scoring_matrix(pos_alignments,neg_alignments, 'data/scoring/PAM100', 5, 3, 1)
    for index_1 in range(output[2].shape[0]):
        for index_2 in range(index_1,output[2].shape[0]):
            assert output[2].iloc[index_1,index_2] == output[2].iloc[index_2,index_1]
