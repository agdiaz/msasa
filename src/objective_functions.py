from abc import abstractmethod, ABCMeta
from itertools import combinations
# from Bio.SubsMat.MatrixInfo import blosum62 as blosum
# from input_parser import InputParser
# import numpy as np
from numpy import array, apply_along_axis, unique, stack
# blosum.update(((b,a),val) for (a,b),val in list(blosum.items()))

GAP = '-'
B_GAP = b'-'
ALL_GAPS = [B_GAP]

char_blosum62 = {
    (b'W', b'F'): 1, (b'L', b'R'): -2, (b'S', b'P'): -1, (b'V', b'T'): 0,
    (b'Q', b'Q'): 5, (b'N', b'A'): -2, (b'Z', b'Y'): -2, (b'W', b'R'): -3,
    (b'Q', b'A'): -1, (b'S', b'D'): 0, (b'H', b'H'): 8, (b'S', b'H'): -1,
    (b'H', b'D'): -1, (b'L', b'N'): -3, (b'W', b'A'): -3, (b'Y', b'M'): -1,
    (b'G', b'R'): -2, (b'Y', b'I'): -1, (b'Y', b'E'): -2, (b'B', b'Y'): -3,
    (b'Y', b'A'): -2, (b'V', b'D'): -3, (b'B', b'S'): 0, (b'Y', b'Y'): 7,
    (b'G', b'N'): 0, (b'E', b'C'): -4, (b'Y', b'Q'): -1, (b'Z', b'Z'): 4,
    (b'V', b'A'): 0, (b'C', b'C'): 9, (b'M', b'R'): -1, (b'V', b'E'): -2,
    (b'T', b'N'): 0, (b'P', b'P'): 7, (b'V', b'I'): 3, (b'V', b'S'): -2,
    (b'Z', b'P'): -1, (b'V', b'M'): 1, (b'T', b'F'): -2, (b'V', b'Q'): -2,
    (b'K', b'K'): 5, (b'P', b'D'): -1, (b'I', b'H'): -3, (b'I', b'D'): -3,
    (b'T', b'R'): -1, (b'P', b'L'): -3, (b'K', b'G'): -2, (b'M', b'N'): -2,
    (b'P', b'H'): -2, (b'F', b'Q'): -3, (b'Z', b'G'): -2, (b'X', b'L'): -1,
    (b'T', b'M'): -1, (b'Z', b'C'): -3, (b'X', b'H'): -1, (b'D', b'R'): -2,
    (b'B', b'W'): -4, (b'X', b'D'): -1, (b'Z', b'K'): 1, (b'F', b'A'): -2,
    (b'Z', b'W'): -3, (b'F', b'E'): -3, (b'D', b'N'): 1, (b'B', b'K'): 0,
    (b'X', b'X'): -1, (b'F', b'I'): 0, (b'B', b'G'): -1, (b'X', b'T'): 0,
    (b'F', b'M'): 0, (b'B', b'C'): -3, (b'Z', b'I'): -3, (b'Z', b'V'): -2,
    (b'S', b'S'): 4, (b'L', b'Q'): -2, (b'W', b'E'): -3, (b'Q', b'R'): 1,
    (b'N', b'N'): 6, (b'W', b'M'): -1, (b'Q', b'C'): -3, (b'W', b'I'): -3,
    (b'S', b'C'): -1, (b'L', b'A'): -1, (b'S', b'G'): 0, (b'L', b'E'): -3,
    (b'W', b'Q'): -2, (b'H', b'G'): -2, (b'S', b'K'): 0, (b'Q', b'N'): 0,
    (b'N', b'R'): 0, (b'H', b'C'): -3, (b'Y', b'N'): -2, (b'G', b'Q'): -2,
    (b'Y', b'F'): 3, (b'C', b'A'): 0, (b'V', b'L'): 1, (b'G', b'E'): -2,
    (b'G', b'A'): 0, (b'K', b'R'): 2, (b'E', b'D'): 2, (b'Y', b'R'): -2,
    (b'M', b'Q'): 0, (b'T', b'I'): -1, (b'C', b'D'): -3, (b'V', b'F'): -1,
    (b'T', b'A'): 0, (b'T', b'P'): -1, (b'B', b'P'): -2, (b'T', b'E'): -1,
    (b'V', b'N'): -3, (b'P', b'G'): -2, (b'M', b'A'): -1, (b'K', b'H'): -1,
    (b'V', b'R'): -3, (b'P', b'C'): -3, (b'M', b'E'): -2, (b'K', b'L'): -2,
    (b'V', b'V'): 4, (b'M', b'I'): 1, (b'T', b'Q'): -1, (b'I', b'G'): -4,
    (b'P', b'K'): -1, (b'M', b'M'): 5, (b'K', b'D'): -1, (b'I', b'C'): -1,
    (b'Z', b'D'): 1, (b'F', b'R'): -3, (b'X', b'K'): -1, (b'Q', b'D'): 0,
    (b'X', b'G'): -1, (b'Z', b'L'): -3, (b'X', b'C'): -2, (b'Z', b'H'): 0,
    (b'B', b'L'): -4, (b'B', b'H'): 0, (b'F', b'F'): 6, (b'X', b'W'): -2,
    (b'B', b'D'): 4, (b'D', b'A'): -2, (b'S', b'L'): -2, (b'X', b'S'): 0,
    (b'F', b'N'): -3, (b'S', b'R'): -1, (b'W', b'D'): -4, (b'V', b'Y'): -1,
    (b'W', b'L'): -2, (b'H', b'R'): 0, (b'W', b'H'): -2, (b'H', b'N'): 1,
    (b'W', b'T'): -2, (b'T', b'T'): 5, (b'S', b'F'): -2, (b'W', b'P'): -4,
    (b'L', b'D'): -4, (b'B', b'I'): -3, (b'L', b'H'): -3, (b'S', b'N'): 1,
    (b'B', b'T'): -1, (b'L', b'L'): 4, (b'Y', b'K'): -2, (b'E', b'Q'): 2,
    (b'Y', b'G'): -3, (b'Z', b'S'): 0, (b'Y', b'C'): -2, (b'G', b'D'): -1,
    (b'B', b'V'): -3, (b'E', b'A'): -1, (b'Y', b'W'): 2, (b'E', b'E'): 5,
    (b'Y', b'S'): -2, (b'C', b'N'): -3, (b'V', b'C'): -1, (b'T', b'H'): -2,
    (b'P', b'R'): -2, (b'V', b'G'): -3, (b'T', b'L'): -1, (b'V', b'K'): -2,
    (b'K', b'Q'): 1, (b'R', b'A'): -1, (b'I', b'R'): -3, (b'T', b'D'): -1,
    (b'P', b'F'): -4, (b'I', b'N'): -3, (b'K', b'I'): -3, (b'M', b'D'): -3,
    (b'V', b'W'): -3, (b'W', b'W'): 11, (b'M', b'H'): -2, (b'P', b'N'): -2,
    (b'K', b'A'): -1, (b'M', b'L'): 2, (b'K', b'E'): 1, (b'Z', b'E'): 4,
    (b'X', b'N'): -1, (b'Z', b'A'): -1, (b'Z', b'M'): -1, (b'X', b'F'): -1,
    (b'K', b'C'): -3, (b'B', b'Q'): 0, (b'X', b'B'): -1, (b'B', b'M'): -3,
    (b'F', b'C'): -2, (b'Z', b'Q'): 3, (b'X', b'Z'): -1, (b'F', b'G'): -3,
    (b'B', b'E'): 1, (b'X', b'V'): -1, (b'F', b'K'): -3, (b'B', b'A'): -2,
    (b'X', b'R'): -1, (b'D', b'D'): 6, (b'W', b'G'): -2, (b'Z', b'F'): -3,
    (b'S', b'Q'): 0, (b'W', b'C'): -2, (b'W', b'K'): -3, (b'H', b'Q'): 0,
    (b'L', b'C'): -1, (b'W', b'N'): -4, (b'S', b'A'): 1, (b'L', b'G'): -4,
    (b'W', b'S'): -3, (b'S', b'E'): 0, (b'H', b'E'): 0, (b'S', b'I'): -2,
    (b'H', b'A'): -2, (b'S', b'M'): -1, (b'Y', b'L'): -1, (b'Y', b'H'): 2,
    (b'Y', b'D'): -3, (b'E', b'R'): 0, (b'X', b'P'): -2, (b'G', b'G'): 6,
    (b'G', b'C'): -3, (b'E', b'N'): 0, (b'Y', b'T'): -2, (b'Y', b'P'): -3,
    (b'T', b'K'): -1, (b'A', b'A'): 4, (b'P', b'Q'): -1, (b'T', b'C'): -1,
    (b'V', b'H'): -3, (b'T', b'G'): -2, (b'I', b'Q'): -3, (b'Z', b'T'): -1,
    (b'C', b'R'): -3, (b'V', b'P'): -2, (b'P', b'E'): -1, (b'M', b'C'): -1,
    (b'K', b'N'): 0, (b'I', b'I'): 4, (b'P', b'A'): -1, (b'M', b'G'): -3,
    (b'T', b'S'): 1, (b'I', b'E'): -3, (b'P', b'M'): -2, (b'M', b'K'): -1,
    (b'I', b'A'): -1, (b'P', b'I'): -3, (b'R', b'R'): 5, (b'X', b'M'): -1,
    (b'L', b'I'): 2, (b'X', b'I'): -1, (b'Z', b'B'): 1, (b'X', b'E'): -1,
    (b'Z', b'N'): 0, (b'X', b'A'): 0, (b'B', b'R'): -1, (b'B', b'N'): 3,
    (b'F', b'D'): -3, (b'X', b'Y'): -1, (b'Z', b'R'): 0, (b'F', b'H'): -1,
    (b'B', b'F'): -3, (b'F', b'L'): 0, (b'X', b'Q'): -1, (b'B', b'B'): 4
}

class SequencesComparer(metaclass=ABCMeta):
    def calculate_score(self, alignment_dataframe):
        pass

    @abstractmethod
    def np_calculate_score(self, alignment_ndarray):
        pass

    @abstractmethod
    def compare(self, seq_a, seq_b):
        pass

    @abstractmethod
    def np_compare(self, i, combo):
        pass


class SingleMS(SequencesComparer):
    # Groups of identical characters are given 1 points * size of the group,
    # 5 points are deducted for a each gap

    def calculate_score(self, alignment_dataframe):
        pass

    def compare(self, seq_a, seq_b):
        pass

    def np_calculate_score(self, alignment_ndarray):
        columns = alignment_ndarray.shape[1]
        result_objects = apply_along_axis(self.np_compare, axis=0, arr=alignment_ndarray)

        return columns * sum(result_objects)

    def np_compare(self, column):
        a, b = unique(column, return_counts=True)
        c = stack((a, b), axis = 1)

        local_scores = apply_along_axis(lambda elem: int(elem[1]) * (10 if elem[0] == B_GAP else 1), axis=1, arr=c)

        return len(local_scores) * sum(local_scores)


class SingleMatching(SequencesComparer):
    # Groups of identical characters are given 1 points * size of the group,
    # 5 points are deducted for a each gap

    def calculate_score(self, alignment_dataframe):
        pass

    def compare(self, seq_a, seq_b):
        pass

    def np_calculate_score(self, alignment_ndarray):
        columns = alignment_ndarray.shape[1]
        result_objects = apply_along_axis(self.np_compare, axis=0, arr=alignment_ndarray)

        return columns * len(result_objects) * sum(result_objects)

    def np_compare(self, column):
        unique_residues = unique(column)
        unique_residues_count = len(unique_residues)

        # ALL GAPS IN A COLUMN:
        if unique_residues_count == 1 and unique_residues[0] == B_GAP:
            return 100
        else:
            return unique_residues_count


class SingleBlosum(SequencesComparer):
    # Groups of identical characters are given 1 points * size of the group,
    # 5 points are deducted for a each gap

    def calculate_score(self, alignment_dataframe):
        pass

    def compare(self, seq_a, seq_b):
        pass

    def np_calculate_score(self, alignment_ndarray):
        columns = alignment_ndarray.shape[1]
        result_objects = apply_along_axis(self.np_compare, axis=0, arr=alignment_ndarray)

        return columns * sum(result_objects)

    def np_compare(self, column):
        cs = list(combinations(column, 2))
        cs_arr = array(cs)

        c = apply_along_axis(self.score, axis=1, arr=cs_arr)

        return len(list(set(cs))) * sum(c)

    def score(self, residues_tuple):
        a, b = residues_tuple

        if B_GAP in (a, b):
            return 25
        else:
            res = char_blosum62.get((a, b), char_blosum62.get((b, a), -20))
            return res * -1


# class GlobalMs(SequencesComparer):
#     # Identical characters are given 5 points, 4 point is deducted for each non-identical character
#     # 3 points are deducted when opening a gap, and 0.1 points are deducted when extending it

#     def compare(self, seq_a, seq_b):
#         return pairwise2.align.globalms(seq_a, seq_b, 5, -4, -3, -0.1, score_only=True)


#     def np_compare(self, i, combo):
#         pass


# class GlobalMsMin(SequencesComparer):
#     # Identical characters are given 5 points, 4 point is deducted for each non-identical character
#     # 3 points are deducted when opening a gap, and 0.1 points are deducted when extending it

#     def compare(self, seq_a, seq_b):
#         return -1 * pairwise2.align.globalms(seq_a, seq_b, 5, -4, -3, -0.1, score_only=True)


#     def np_compare(self, i, combo):
#         pass


# class Blosum(SequencesComparer):
#     # https://stackoverflow.com/questions/5686211/is-there-a-function-that-can-calculate-a-score-for-aligned-sequences-given-the-a

#     def __init__(self):
#         self.matrix = blosum
#         self.first_gap_score = -0.5
#         self.gap_continuation_score = -0.25
#         self.error_score= -1

#     def compare(self, seq_a, seq_b):
#         return self.__score_pairwise(seq_a, seq_b)

#     def __score_pairwise(self, seq1, seq2):
#         score = 0
#         gap = False

#         for i in range(len(seq1)):
#             pair = (seq1[i], seq2[i])

#             if not gap:
#                 if GAP in pair:
#                     gap = True
#                     score += self.first_gap_score
#                 else:
#                     score += self.__score_match(pair)
#             else:
#                 if GAP not in pair:
#                     gap = False

#                     score += self.__score_match(pair)
#                 else:
#                     score += self.gap_continuation_score
#         return score

#     def __score_match(self, pair):
#         try:
#             if pair in self.matrix:
#                 return self.matrix[pair]
#             else:
#                 return self.matrix[(tuple(reversed(pair)))]
#         except:
#             return self.error_score

#     def np_compare(self, i, combo):
#         pass


# class MatchingCount(SequencesComparer):
#     def __init__(self):
#         self.opening_gap_penalty = 2
#         self.continuation_gap_penalty = 3
#         self.residue_match = -2
#         self.error_penalty = 10
#         self.mismatch_penalty = 0
#         self.half_gap = self.opening_gap_penalty / 2.0


#     def compare(self, seq_a, seq_b):
#         score_total = 0
#         gap_column = False

#         for i in range(len(seq_a)):
#             pos_a = seq_a[i]
#             pos_b = seq_b[i]

#             if pos_a == pos_b and pos_a == GAP:
#                 if gap_column:
#                     score_local = self.continuation_gap_penalty
#                 else:
#                     score_local = self.opening_gap_penalty

#                 gap_column = True
#             elif pos_a == pos_b:
#                 gap_column = False
#                 score_local = self.residue_match
#             elif pos_a == GAP or pos_b == GAP:
#                 gap_column = False
#                 score_local = self.half_gap
#             elif pos_a != pos_b:
#                 score_local = self.mismatch_penalty
#             else:
#                 gap_column = False
#                 score_local = self.error_penalty

#             score_total += score_local

#         return score_total


#     def np_compare(self, i, combo):
#         vec_a, vec_b = combo

#         seq_a = [b for b in vec_a]
#         seq_b = [b for b in vec_b]

#         score_total = 0
#         gap_column = False

#         for i in range(len(seq_a)):
#             pos_a = seq_a[i]
#             pos_b = seq_b[i]

#             if pos_a == pos_b and pos_a == B_GAP:
#                 if gap_column:
#                     score_local = self.continuation_gap_penalty
#                 else:
#                     score_local = self.opening_gap_penalty

#                 gap_column = True
#             elif pos_a == pos_b:
#                 gap_column = False
#                 score_local = self.residue_match
#             elif pos_a == B_GAP or pos_b == B_GAP:
#                 gap_column = False
#                 score_local = self.half_gap
#             elif pos_a != pos_b:
#                 score_local = self.mismatch_penalty
#             else:
#                 gap_column = False
#                 score_local = self.error_penalty

#             score_total += score_local

#         return score_total

class SequencesComparerFactory:
    @staticmethod
    def from_name(name):
        # if name == "global_ms":
        #     return GlobalMs()
        # elif name == "global_ms_min":
        #     return GlobalMsMin()
        # elif name == "blosum":
        #     return Blosum()
        # elif name == "matching":
        #     return MatchingCount()
        if name == "single_ms":
            return SingleMS()
        elif name == "single_blosum":
            return SingleBlosum()
        elif name == "single_matching":
            return SingleMatching()
        else:
            raise NameError("Wrong name: {0}".format(name))