import sys
from .io import read_pairs, read_sequence, read_scoring_matrix, write_alignment, write_optimal_matrix
from .alignment import smith_waterman, smith_waterman_len_adj, smith_waterman_alignment, optimize_scoring_matrix

# Some quick stuff to make sure the program is called correctly
if sys.argv[1] != '-O'
    if len(sys.argv) != 8:
        print("Usage: python -m BMI203_HW3 [-S | -L | -A] <path_seq_a>  <path_seq_b> <path_scoring matrix> <gap_open> <gap_extend> <output_file>")
        sys.exit(0)
else:
    if len(sys.argv) != 9:
        print("Usage: python -m BMI203_HW3 [-O] <alignments_seq_a>  <alignments_seq_b> <path_scoring matrix> <gap_open> <gap_extend> <num_iterations> <output_file>")
        sys.exit(0)

# Choose Smith Waterman algorithm
if sys.argv[1][0:2] == '-S':
    print("Outputting score for Smith Waterman algorithm")
    score = smith_waterman(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    write_alignment(sys.argv[7], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], "Not produced", "Not produced", score)

if sys.argv[1][0:2] == '-L':
    print("Outputting score for length-adjusted Smith Waterman algorithm")
    score = smith_waterman_len_adj(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    write_alignment(sys.argv[7], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], "Not produced", "Not produced", score)

if sys.argv[1][0:2] == '-A':
    print("Outputing alignment and score for Smith Waterman algorithm")
    alignment_a, alignment_b, score = smith_waterman_alignment(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    write_alignment(sys.argv[7], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], alignment_a, alignment_b, score)

# Choose Optimization
if sys.argv[1][0:2] == '-O':
    print("Outputing optimized scoring matrix")
    output = optimize_scoring_matrix(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
    write_optimal_matrix(sys.argv[8], sys.argv[7], sys.argv[4], output[0],output[1],output[2],output[7])
