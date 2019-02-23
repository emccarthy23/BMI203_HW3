import os
import numpy as np
import pandas as pd
import pandas as pd
import numpy as np
import os
from random import uniform
import time
import random



def read_pairs(file):
    #Open the file
    f = open ( file , 'r')
    #Read all the lines
    lines = f.read().splitlines()
    #Create a list of the paired filenames
    files = [line.split()for line in lines]
    return files

def read_sequence(file):
    #Store protein name
    protein_name = file.split('/')[-1][:-3]
    #Open the file
    f = open ( file , 'r')
    #Read all the lines
    lines = f.read().splitlines()
    #Remove the header
    lines = lines[1:]
    #Initialize array to store sequence
    sequence = ''
    for line in lines:
        sequence = sequence + line
    sequence = sequence.upper()
    return protein_name, sequence

def read_scoring_matrix(file):

    #Open the file
    f = open ( file , 'r')
    #Read all the lines
    lines = f.read().splitlines()
    #Find where the header ends
    aa_names_index = 0
    while lines[aa_names_index][0] == '#':
        aa_names_index = aa_names_index +1
    #Remove header
    lines = lines[aa_names_index:]
    #Initialize DataFrame to store scoring matrix
    amino_acids = lines[0].split()
    scoring_df = pd.DataFrame(index=amino_acids, columns=amino_acids)
    #Store scores into dataframe
    for line_index in range(1,len(lines)):
        scoring_df.iloc[line_index-1] = lines[line_index].split()
    #Convert scores to numbers from strings
    scoring_df = scoring_df.apply(pd.to_numeric)
    return scoring_df

def write_alignment(filename, path_seq_a, path_seq_b, path_scoring_matrix, gap_open, gap_extend, alignment_seqa, alignment_seqb, score):
    """
    Write Smith Waterman alignment for two sequences out to a file.

    Input: a filename and Smith Waterman alignment inputs and outputs
    Output: none
    """

    out = open(filename, 'w')
    out.write("\nScoring Matrix %s\n------------\n" % path_scoring_matrix)
    out.write("\nSequence A %s\n------------\n" % path_seq_a)
    out.write("\nSequence B %s\n------------\n" % path_seq_b)
    out.write("\nGap open penalty %d\n------------\n" % gap_open)
    out.write("\nGap extension penalty %d\n------------\n" % gap_extend)
    out.write("\nScore %d\n------------\n" % score)

    out.write("\nLocal alignment for Seq A %s\n------------\n" % alignment_seqa)
    out.write("\nLocal alignment for Seq B %s\n------------\n" % alignment_seqa)

    out.close()

new_matrix, new_matrix_score, original_score, iteration_counter, iteration_score_counter, total_time
def write_optimal_matrix(filename, k, path, opt_matrix, opt_matrix_score, original_score,time):
    """
    Write the optimal scoring matrix out to a file.

    Input: a filename and the output from optimize_scoring_matrix
    Output: none
    """

    out = open(filename, 'w')
    out.write("# Optimal Matrix after %d iterations on %s scoring matrix\n" % (k, path))
    out.write("# Time (min) %d \n" % time)
    out.write("# Original matrix score %d \n" % original_score)
    out.write("# New matrix score %d \n" % opt_matrix_score)
    out.close()
    opt_matrix.to_csv(filename, index=None, sep=' ', mode='a')
