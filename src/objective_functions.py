from itertools import combinations
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62 as blosum
from input_parser import InputParser

blosum.update(((b,a),val) for (a,b),val in list(blosum.items()))

# def msa_objective_real(df, sequences_comparer_name = "global_ms"):

#     sequences_comparer = SequencesComparerFactory.from_name(sequences_comparer_name)

#     score_total = 0
#     sequences = InputParser.dataframe_to_sequences(df)
#     seq_combinations = combinations(sequences, 2)

#     for combo in seq_combinations:
#         seq_a, seq_b = combo

#         combination_score = sequences_comparer.compare(seq_a, seq_b)
#         score_total += combination_score

#     return score_total


class SequencesComparerFactory:
    @staticmethod
    def from_name(name):
        if name == "global_ms":
            return GlobalMs()
        elif name == "blosum":
            return Blosum()
        else:
            raise "Wrong name"

class SequencesComparer:
    def calculate_score(self, alignment_dataframe):
        score_total = 0
        sequences = InputParser.dataframe_to_sequences(alignment_dataframe)
        seq_combinations = combinations(sequences, 2)

        for combo in seq_combinations:
            seq_a, seq_b = combo

            combination_score = self.compare(seq_a, seq_b)
            score_total += combination_score

        return score_total

    def compare(self, seq_a, seq_b):
        pass

class GlobalMs(SequencesComparer):
    # Identical characters are given 5 points, 4 point is deducted for each non-identical character
    # 3 points are deducted when opening a gap, and 0.1 points are deducted when extending it

    def compare(self, seq_a, seq_b):
        return pairwise2.align.globalms(seq_a, seq_b, 5, -4, -3, -0.1, score_only=True)

class Blosum(SequencesComparer):
    # https://stackoverflow.com/questions/5686211/is-there-a-function-that-can-calculate-a-score-for-aligned-sequences-given-the-a

    def __init__(self):
        self.matrix = blosum
        self.first_gap_score = -2
        self.gap_continuation_score = -1
        self.error_score=-10

    def compare(self, seq_a, seq_b):
        return self.__score_pairwise(seq_a, seq_b)

    def __score_pairwise(self, seq1, seq2):
        score = 0
        gap = False

        for i in range(len(seq1)):
            pair = (seq1[i], seq2[i])

            if not gap:
                if '-' in pair:
                    gap = True
                    score += self.first_gap_score
                else:
                    score += self.__score_match(pair)
            else:
                if '-' not in pair:
                    gap = False

                    score += self.__score_match(pair)
                else:
                    score += self.gap_continuation_score
        return score

    def __score_match(self, pair):
        try:
            if pair in self.matrix:
                return self.matrix[pair]
            else:
                return self.matrix[(tuple(reversed(pair)))]
        except:
            return self.error_score
