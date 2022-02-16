from itertools import combinations
from Bio import pairwise2

from input_parser import InputParser

def msa_objective_real(df):
    total_energy = 0
    sequences = InputParser.dataframe_to_sequences(df)
    seq_combinations = combinations(sequences, 2)

    for combo in seq_combinations:
        seq_a, seq_b = combo
        # Identical characters are given 5 points, 4 point is deducted for each non-identical character
        # 3 points are deducted when opening a gap, and 0.1 points are deducted when extending it
        blosum_score = pairwise2.align.globalms(seq_a, seq_b, 5, -4, -3, -0.1, score_only=True)
        total_energy += blosum_score

    return total_energy
