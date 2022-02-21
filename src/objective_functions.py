from abc import abstractmethod
from itertools import combinations
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62 as blosum
from input_parser import InputParser

blosum.update(((b,a),val) for (a,b),val in list(blosum.items()))


class SequencesComparerFactory:
    @staticmethod
    def from_name(name):
        if name == "global_ms":
            return GlobalMs()
        elif name == "global_ms_min":
            return GlobalMsMin()
        elif name == "blosum":
            return Blosum()
        elif name == "matching":
            return MatchingCount()
        else:
            raise NameError("Wrong name: {0}".format(name))


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

    def np_calculate_score(self, alignment_nparray):
        score_total = 0
        seq_combinations = combinations(alignment_nparray, 2)

        for combo in seq_combinations:
            raw_seq_a, raw_seq_b = combo

            combination_score = self.np_compare(raw_seq_a, raw_seq_b)
            score_total += combination_score

        return score_total


    @abstractmethod
    def compare(self, seq_a, seq_b) -> float:
        pass

    @abstractmethod
    def np_compare(self, seq_a, seq_b) -> float:
        pass

class GlobalMs(SequencesComparer):
    # Identical characters are given 5 points, 4 point is deducted for each non-identical character
    # 3 points are deducted when opening a gap, and 0.1 points are deducted when extending it

    def compare(self, seq_a, seq_b):
        return pairwise2.align.globalms(seq_a, seq_b, 5, -4, -3, -0.1, score_only=True)


    def np_compare(self, seq_a, seq_b) -> float:
        pass

class GlobalMsMin(SequencesComparer):
    # Identical characters are given 5 points, 4 point is deducted for each non-identical character
    # 3 points are deducted when opening a gap, and 0.1 points are deducted when extending it

    def compare(self, seq_a, seq_b):
        return -1 * pairwise2.align.globalms(seq_a, seq_b, 5, -4, -3, -0.1, score_only=True)


    def np_compare(self, seq_a, seq_b) -> float:
        pass

class Blosum(SequencesComparer):
    # https://stackoverflow.com/questions/5686211/is-there-a-function-that-can-calculate-a-score-for-aligned-sequences-given-the-a

    def __init__(self):
        self.matrix = blosum
        self.first_gap_score = -0.5
        self.gap_continuation_score = -0.25
        self.error_score= -1

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

    def np_compare(self, seq_a, seq_b) -> float:
        pass

class MatchingCount(SequencesComparer):
    def __init__(self):
        self.opening_gap_penalty = 2
        self.continuation_gap_penalty = 3
        self.residue_match = -1
        self.error_penalty = 10
        self.mismatch_penalty = 1
        self.half_gap = self.opening_gap_penalty / 2.0


    def compare(self, seq_a, seq_b) -> float:
        score_total = 0
        gap_column = False

        for i in range(len(seq_a)):
            pos_a = seq_a[i]
            pos_b = seq_b[i]

            if pos_a == pos_b and pos_a == '-':
                if gap_column:
                    score_local = self.continuation_gap_penalty
                else:
                    score_local = self.opening_gap_penalty

                gap_column = True
            elif pos_a == pos_b:
                gap_column = False
                score_local = self.residue_match
            elif pos_a == '-' or pos_b == '-':
                gap_column = False
                score_local = self.half_gap
            elif pos_a != pos_b:
                score_local = self.mismatch_penalty
            else:
                gap_column = False
                score_local = self.error_penalty

            score_total += score_local

        return score_total


    def np_compare(self, seq_a, seq_b) -> float:
        score_total = 0
        gap_column = False

        for i in range(len(seq_a)):
            pos_a = seq_a[i]
            pos_b = seq_b[i]

            if pos_a == pos_b and pos_a == b'-':
                if gap_column:
                    score_local = self.continuation_gap_penalty
                else:
                    score_local = self.opening_gap_penalty

                gap_column = True
            elif pos_a == pos_b:
                gap_column = False
                score_local = self.residue_match
            elif pos_a == '-' or pos_b == b'-':
                gap_column = False
                score_local = self.half_gap
            elif pos_a != pos_b:
                score_local = self.mismatch_penalty
            else:
                gap_column = False
                score_local = self.error_penalty

            score_total += score_local

        return score_total
