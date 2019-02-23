from .io import read_pairs, read_sequence, read_scoring_matrix, write_alignment, write_optimal_matrix
import pandas as pd
import numpy as np
import os
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
from random import uniform
import time
import random


def smith_waterman_alignment(path_seq_a, path_seq_b, path_scoring_matrix, gap_open, gap_extend):
    #Read sequences and scoring matrix
    name_a, seq_a = read_sequence(path_seq_a)
    name_b, seq_b = read_sequence(path_seq_b)
    cost_matrix = read_scoring_matrix(path_scoring_matrix)

    #Store value for -Inf
    MIN = -float("inf")

    """
    We will initialize the values for three matrices:
    A - score for aligning seq_a and seq_b up to i,j with a gap in seq_a at position i
    B - score for aligning seq_a and seq_b up to i,j with a gap in seq_b at position j
    M - score for aligning seq_a and seq_b with a alignment at position i,j
    """
    dim_a = len(seq_a)+1
    dim_b = len(seq_b)+1

    ####Initialize matrices

    A = np.zeros((dim_a, dim_b))
    B = np.zeros((dim_a, dim_b))
    M = np.zeros((dim_a, dim_b))

    A_direction = {}
    B_direction = {}
    M_direction = {}
    for i in range(dim_a):
        A_direction[i,0] = "end"
        B_direction[i,0] = "end"
        M_direction[i,0] = "end"
    for j in range(dim_b):
        A_direction[0,j] = "end"
        B_direction[0,j] = "end"
        M_direction[0,j] = "end"

    #Fill in top row and left row for each gap matrix with -Inf because we will not allow
    #a gap at the start for our local alignment
    for i in range(1,dim_a):
        A[i,0] = MIN
        B[i,0] = MIN
    for j in range(1,dim_b):
        A[0,j] = MIN
        B[0,j] = MIN


    #Make list to keep track of direction
    event_gap_a = ["open_gap_a", "extend_gap_a", "open_gap_a_from_b"]
    event_gap_b = ["open_gap_b", "extend_gap_b", "open_gap_b_from_a"]
    event_match = ["match", "close_gap_a","close_gap_b","end"]
    """
    Now we will fill in the values for these three matrices
    """


    for i in range(1,dim_a):
        for j in range(1,dim_b):
            #For A (putting a gap in seq_a at position i), we have three possibilities for how seq_a up to i
            #and seq_b up to j-1 (since j is aligning to the gap in seq_a) could have been aligned. We can come from:
            # 1) a previous alignment
            # 2) a previous gap in seq_a
            # 3) a previous gap in seq_b (so a new gap in seq_a)
            #dic_A["M":M[i,j-1] + gap_open + gap_extend, "A":gap_extend + A[i,j-1]]
            values = [M[i,j-1] - gap_open - gap_extend ,  A[i,j-1] - gap_extend, B[i,j-1] - gap_open - gap_extend]
            A[i,j] = max(values)
            A_direction[i,j] = event_gap_a[values.index(A[i,j])]
            #For B (putting a gap in seq_b at position j), we have three possibilities for how seq_a up to i-1 (since
            # i is aligning with a gap) and seq_b up to j could have been aligned. We can come from:
            # 1) a previous alignment
            # 2) a previous gap in seq_a (so a new gap in seq_b)
            # 3) a previous gap in seq_b
            values = [M[i-1,j] - gap_open - gap_extend , B[i-1,j] - gap_extend, A[i-1,j] - gap_open - gap_extend, ]
            B[i,j] = max(values)
            B_direction[i,j] = event_gap_b[values.index(B[i,j])]
            #For M alinging position i and j from seq_a and seq_b respectively we can come from:
            # 1) a previous alignment
            # 2) a previous gap in seq_a
            # 3) a previous gap in seq_b
            #Cost for aligning seq_a and seq_b at position i,j (need to account for zero indexing)
            #Let 0 be the minimum score in order to create local rather than global optimum alignments
            values = [cost_matrix.loc[seq_a[i-1],seq_b[j-1]] + M[i-1,j-1], cost_matrix.loc[seq_a[i-1],seq_b[j-1]] + A[i-1,j-1], cost_matrix.loc[seq_a[i-1],seq_b[j-1]] + B[i-1,j-1],0]
            M[i,j] = max(values)
            M_direction[i,j] = event_match[values.index(max(values))]
    if M.max() == 0:
        return "","",0

    #Do traceback to get aligned sequence
    #Initalize strings to contan alignment

    #Find index of max score position in Z
    indices = np.where(M == M.max())

    #If there are multiple alignments take the first position
    index_a = indices[0][0]
    index_b = indices[1][0]
    #Store first traceback value and initial direction
    alignment_seqa, alignment_seqb = seq_a[index_a-1], seq_b[index_b-1]

    direction = M_direction[index_a,index_b]

    #Do traceback and store sequence

    while direction != "end":
        #Move in recorded direction and update alligned sequence
        if index_a == 1 or index_b == 1:
            break
        elif direction == 'close_gap_a':
            index_a, index_b = index_a-1, index_b-1
            alignment_seqa, alignment_seqb = "-"+ alignment_seqa, seq_b[index_b-1]+ alignment_seqb
            direction = A_direction[index_a, index_b]
        elif direction == 'close_gap_b':
            index_a, index_b = index_a-1, index_b-1
            alignment_seqa, alignment_seqb = seq_a[index_a-1]+ alignment_seqa, "-"+ alignment_seqb
            direction = B_direction[index_a, index_b]
        elif direction == 'match':
            index_a, index_b = index_a-1, index_b-1
            alignment_seqa, alignment_seqb = seq_a[index_a-1] + alignment_seqa, seq_b[index_b-1]+ alignment_seqb
            direction = M_direction[index_a, index_b]
        elif direction == 'end':
            break
        elif direction == 'open_gap_a':
            index_a, index_b = index_a, index_b-1
            alignment_seqa, alignment_seqb = seq_a[index_a-1]+ alignment_seqa,seq_b[index_b-1]+ alignment_seqb
            direction = M_direction[index_a, index_b]
        elif direction == 'extend_gap_a':
            index_a, index_b = index_a, index_b-1
            alignment_seqa, alignment_seqb = "-"+ alignment_seqa,seq_b[index_b-1]+ alignment_seqb
            direction = A_direction[index_a, index_b]
        elif direction == 'open_gap_a_from_b':
            index_a, index_b = index_a, index_b-1
            alignment_seqa, alignment_seqb = seq_a[index_a-1]+ alignment_seqa,"-"+ alignment_seqb
            direction = B_direction[index_a, index_b]
        elif direction == 'open_gap_b':
            index_a, index_b = index_a-1, index_b
            alignment_seqa, alignment_seqb = seq_a[index_a-1]+ alignment_seqa,seq_b[index_b-1]+ alignment_seqb
            direction = M_direction[index_a, index_b]
        elif direction == 'extend_gap_b':
            index_a, index_b = index_a-1, index_b
            alignment_seqa, alignment_seqb = seq_a[index_a-1]+ alignment_seqa,"-"+ alignment_seqb
            direction = B_direction[index_a, index_b]
        elif direction == 'open_gap_b_from_a':
            index_a, index_b = index_a-1, index_b
            alignment_seqa, alignment_seqb = "-"+ alignment_seqa,seq_b[index_b-1]+ alignment_seqb
            direction = A_direction[index_a, index_b]
    #Need to remove last match if it actually was making the alignment worse
    if cost_matrix.loc[seq_a[index_a-1],seq_b[index_b-1]] < 0:
        return alignment_seqa[1:], alignment_seqb[1:], M.max()
    else:
        return alignment_seqa, alignment_seqb, M.max()

def smith_waterman(path_seq_a, path_seq_b, path_scoring_matrix, gap_open, gap_extend):
    #Read sequences and scoring matrix
    name_a, seq_a = read_sequence(path_seq_a)
    name_b, seq_b = read_sequence(path_seq_b)
    cost_matrix = read_scoring_matrix(path_scoring_matrix)

    #Store value for -Inf
    MIN = -float("inf")

    """
    We will initialize the values for three matrices:
    A - score for aligning seq_a and seq_b up to i,j with a gap in seq_a at position i
    B - score for aligning seq_a and seq_b up to i,j with a gap in seq_b at position j
    M - score for aligning seq_a and seq_b with a alignment at position i,j
    """
    dim_a = len(seq_a)+1
    dim_b = len(seq_b)+1

    ####Initialize matrices

    A = np.zeros((dim_a, dim_b))
    B = np.zeros((dim_a, dim_b))
    M = np.zeros((dim_a, dim_b))

    #Fill in top row and left row for each gap matrix with -Inf because we will not allow
    #a gap at the start for our local alignment
    for i in range(1,dim_a):
        A[i,0] = MIN
        B[i,0] = MIN
    for j in range(1,dim_b):
        A[0,j] = MIN
        B[0,j] = MIN

    """
    Now we will fill in the values for these three matrices
    """


    for i in range(1,dim_a):
        for j in range(1,dim_b):
            #For A (putting a gap in seq_a at position i), we have three possibilities for how seq_a up to i
            #and seq_b up to j-1 (since j is aligning to the gap in seq_a) could have been aligned. We can come from:
            # 1) a previous alignment
            # 2) a previous gap in seq_a
            # 3) a previous gap in seq_b (so a new gap in seq_a)
            #dic_A["M":M[i,j-1] + gap_open + gap_extend, "A":gap_extend + A[i,j-1]]
            A[i,j] = max(M[i,j-1] - gap_open - gap_extend ,  A[i,j-1] - gap_extend, B[i,j-1] - gap_open - gap_extend)


            #For B (putting a gap in seq_b at position j), we have three possibilities for how seq_a up to i-1 (since
            # i is aligning with a gap) and seq_b up to j could have been aligned. We can come from:
            # 1) a previous alignment
            # 2) a previous gap in seq_a (so a new gap in seq_b)
            # 3) a previous gap in seq_b
            B[i,j] = max(M[i-1,j] - gap_open - gap_extend , B[i-1,j] - gap_extend, A[i-1,j] - gap_open - gap_extend)


            #For M alinging position i and j from seq_a and seq_b respectively we can come from:
            # 1) a previous alignment
            # 2) a previous gap in seq_a
            # 3) a previous gap in seq_b
            #Cost for aligning seq_a and seq_b at position i,j (need to account for zero indexing)
            #Let 0 be the minimum score in order to create local rather than global optimum alignments
            M[i,j] = max(cost_matrix.loc[seq_a[i-1],seq_b[j-1]] + M[i-1,j-1], cost_matrix.loc[seq_a[i-1],seq_b[j-1]] + A[i-1,j-1], cost_matrix.loc[seq_a[i-1],seq_b[j-1]] + B[i-1,j-1],0)

    return np.max(M)


def smith_waterman_len_adj(path_seq_a, path_seq_b, path_scoring_matrix, gap_open, gap_extend):
    #Read sequences and scoring matrix
    name_a, seq_a = read_sequence(path_seq_a)
    name_b, seq_b = read_sequence(path_seq_b)
    cost_matrix = read_scoring_matrix(path_scoring_matrix)

    #Store value for -Inf
    MIN = -float("inf")

    """
    We will initialize the values for three matrices:
    A - score for aligning seq_a and seq_b up to i,j with a gap in seq_a at position i
    B - score for aligning seq_a and seq_b up to i,j with a gap in seq_b at position j
    M - score for aligning seq_a and seq_b with a alignment at position i,j
    """
    dim_a = len(seq_a)+1
    dim_b = len(seq_b)+1

    ####Initialize matrices

    A = np.zeros((dim_a, dim_b))
    B = np.zeros((dim_a, dim_b))
    M = np.zeros((dim_a, dim_b))

    #Fill in top row and left row for each gap matrix with -Inf because we will not allow
    #a gap at the start for our local alignment
    for i in range(1,dim_a):
        A[i,0] = MIN
        B[i,0] = MIN
    for j in range(1,dim_b):
        A[0,j] = MIN
        B[0,j] = MIN

    """
    Now we will fill in the values for these three matrices
    """


    for i in range(1,dim_a):
        for j in range(1,dim_b):
            #For A (putting a gap in seq_a at position i), we have three possibilities for how seq_a up to i
            #and seq_b up to j-1 (since j is aligning to the gap in seq_a) could have been aligned. We can come from:
            # 1) a previous alignment
            # 2) a previous gap in seq_a
            # 3) a previous gap in seq_b (so a new gap in seq_a)
            #dic_A["M":M[i,j-1] + gap_open + gap_extend, "A":gap_extend + A[i,j-1]]
            A[i,j] = max(M[i,j-1] - gap_open - gap_extend ,  A[i,j-1] - gap_extend, B[i,j-1] - gap_open - gap_extend)


            #For B (putting a gap in seq_b at position j), we have three possibilities for how seq_a up to i-1 (since
            # i is aligning with a gap) and seq_b up to j could have been aligned. We can come from:
            # 1) a previous alignment
            # 2) a previous gap in seq_a (so a new gap in seq_b)
            # 3) a previous gap in seq_b
            B[i,j] = max(M[i-1,j] - gap_open - gap_extend , B[i-1,j] - gap_extend, A[i-1,j] - gap_open - gap_extend)


            #For M alinging position i and j from seq_a and seq_b respectively we can come from:
            # 1) a previous alignment
            # 2) a previous gap in seq_a
            # 3) a previous gap in seq_b
            #Cost for aligning seq_a and seq_b at position i,j (need to account for zero indexing)
            #Let 0 be the minimum score in order to create local rather than global optimum alignments
            M[i,j] = max(cost_matrix.loc[seq_a[i-1],seq_b[j-1]] + M[i-1,j-1], cost_matrix.loc[seq_a[i-1],seq_b[j-1]] + A[i-1,j-1], cost_matrix.loc[seq_a[i-1],seq_b[j-1]] + B[i-1,j-1],0)

    return np.max(M)/min(len(seq_a),len(seq_b))

def score_alignment(seq_a, seq_b, scoring_matrix,gap,ext):
    if len(seq_a) != len(seq_b):
        return "Please return valid local alignments"
    else:
        i = 0
        score = 0
        while i in range(len(seq_a)):
            if seq_a[i] != "-" and seq_b[i] != "-":
                score = scoring_matrix.loc[seq_a[i],seq_b[i]] + score
            elif seq_a[i] == "-" and seq_b[i] != "-":
                if seq_a[i-1] != "-":
                    score = -gap -ext + score
                else:
                    score = -ext + score
            elif seq_b[i] == "-": #and seq_a[i] != "-":
                if seq_b[i-1] != "-":
                    score = -gap -ext + score
                else:
                    score = -ext + score
            elif seq_b[i] == "-" and seq_a[i] == "-":
                if seq_b[i-1] == "-" and seq_a[i-1] == "-":
                    score = -ext*2 + score
                elif seq_b[i-1] != "-" and seq_a[i-1] != "-":
                    score = -ext*2 + -gap*2 + score
                else:
                    score = -ext*2 + -gap + score
            #print(score)
            i += 1
    return score


def score_performance(pos_sequences,neq_sequences,scoring_matrix, gap,ext):
    ##Caluclate the performance for the current scoring matrix
    pos_scores = [score_alignment(x[0], x[1], scoring_matrix, gap, ext) for x in pos_sequences]
    neg_scores = [score_alignment(x[0], x[1], scoring_matrix, gap, ext) for x in neq_sequences]
    #Find thresholds for 0, 0.1, 0.2, and 0.3 false positive rate. (there are 50 Pos_pairs)
    thresholds = []
    thresholds.append(sorted(neg_scores)[-1])
    thresholds.append(sorted(neg_scores)[-6])
    thresholds.append(sorted(neg_scores)[-11])
    thresholds.append(sorted(neg_scores)[-16])
    #Calculate true_pos rate at each score
    false_pos = []
    true_pos = []
    for value in thresholds:
        false_pos.append(np.sum(np.array(neg_scores) > value)/len(neg_scores))
        true_pos.append(np.sum(np.array(pos_scores) > value)/len(pos_scores))
    overall_score = np.sum(true_pos)
    return overall_score, false_pos, true_pos

def optimize_scoring_matrix(alignments_pos, alignments_neg, starting_matrix_path, gap_open, gap_ext, num_iterations):
    start = time.time()
    score_mat = read_scoring_matrix(starting_matrix_path)
    #Store original matrix
    original_mat = score_mat
    #I will perform a Monte Carlo inspired optimization
    #Intialize the iteration matrices and iteration score lists
    iteration_mat = [score_mat]*10

    original_score = score_performance(alignments_pos, alignments_neg, score_mat, gap_open, gap_ext)[0]
    iteration_scores = [original_score]*10
    print("Loaded initial matrix and scores")

#     for k in range(99):
#         for i in range(score_mat.shape[0]):
#             for j in range(i,score_mat.shape[0]):
#                 score_mat.iloc[i,j] = score_mat.iloc[i,j] + uniform(-1,1)
#                 score_mat.iloc[j,i] = score_mat.iloc[i,j]
#         iteration_scores.append(score_performance(alignments_pos, alignments_neg, score_mat, gap_open, gap_ext)[0])
#         iteration_mat.append(score_mat)
#         score_mat = original_mat
    print("Created matrix family")
    #Perform optimization with pool of matrices
    #Take top 10 scoring matrices
    #Mutate each one ten times to create pool of 100 matrices
    #Repeat
    iteration_counter = []
    iteration_counter.append(0)
    iteration_score_counter = []
    iteration_score_counter.append([np.mean(iteration_scores),np.std(iteration_scores),np.max(iteration_scores)])
    for i in range(num_iterations):
        #Find top 10 matrices
        best_mat_scores = sorted(iteration_scores)[-10:]
        best_mat_indices = [iteration_scores.index(x) for x in best_mat_scores]
        new_matrices = []
        original_scores = iteration_scores
        iteration_scores = [0]*100
        #Make list of 100 matrices with each of the best ten repeated ten times
        for index in best_mat_indices:
            new_matrices = new_matrices + [iteration_mat[index]]*10
        #Mutate the 9 copies but keep one of each original
        for additive in range(10):
            iteration_scores[9+10*additive]= (original_scores[best_mat_indices[additive]])
            for mat_index in range(9):
                mat = new_matrices[additive+mat_index]
                for b in range(mat.shape[0]):
                    for c in range(b,mat.shape[0]):
                        mat.iloc[b,c] = mat.iloc[b,c] + uniform(-1,1)
                        mat.iloc[c,b] = mat.iloc[b,c]
                new_matrices[additive*10+mat_index] = mat
                iteration_scores[additive*10+mat_index] = score_performance(alignments_pos, alignments_neg, mat, gap_open, gap_ext)[0]
        iteration_mat = new_matrices
        iteration_counter.append(i+1)
        iteration_score_counter.append([np.mean(iteration_scores),np.std(iteration_scores), np.max(iteration_scores)])
        print("Finished iteration", i+1)
    new_matrix_score = max(iteration_scores)
    new_matrix = iteration_mat[iteration_scores.index(new_matrix_score)]
    end = time.time()
    total_time = (end-start)/60
    return iteration_scores, iteration_mat, new_matrix, new_matrix_score, original_score, iteration_counter, iteration_score_counter, total_time
