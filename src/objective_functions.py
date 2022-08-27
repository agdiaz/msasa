from abc import abstractmethod, ABCMeta
from itertools import combinations
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62 as blosum
from input_parser import InputParser
import multiprocessing as mp
import numpy as np

blosum.update(((b,a),val) for (a,b),val in list(blosum.items()))

GAP = '-'
B_GAP = b'-'


class SequencesComparer(metaclass=ABCMeta):
    def calculate_score(self, alignment_dataframe):
        score_total = 0
        sequences = InputParser.dataframe_to_sequences(alignment_dataframe)
        seq_combinations = combinations(sequences, 2)

        for combo in seq_combinations:
            seq_a, seq_b = combo

            combination_score = self.compare(seq_a, seq_b)
            score_total += combination_score

        return score_total

    def np_calculate_score(self, alignment_ndarray):
        pool = mp.Pool(mp.cpu_count())

        seq_combinations = combinations(alignment_ndarray, 2)

        result_objects = [pool.apply_async(self.np_compare, args=(i, combo)) for i, combo in enumerate(seq_combinations)]
        results = [r.get() for r in result_objects]

        pool.join()
        pool.close()

        return sum(results)

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

    def np_calculate_score(self, alignment_ndarray):
        # pool = mp.Pool(mp.cpu_count())
        columns = alignment_ndarray.transpose()

        # result_objects = [pool.apply_async(self.np_compare, args=(i, combo)) for i, combo in enumerate(columns)]
        # results = [r.get() for r in result_objects]

        results = []
        for i, combo in enumerate(columns):
            r = self.np_compare(i, combo)
            results.append(r)

        return sum(results)

    def compare(self, seq_a, seq_b):
        pass

    def np_compare(self, i, column):
        unique, counts = np.unique(column, return_counts=True)
        listOfUniqueValues = zip(unique, counts)

        groups = {}
        total_score = 0 # The less groups, the better score

        for elem in listOfUniqueValues:
            multiplier = 10 if elem[0] == B_GAP else 1
            local_score = multiplier * elem[1]
            total_score += local_score
            groups[elem[0]] = { 'count': elem[1], 'score': local_score }

        return len(unique) * total_score

class GlobalMs(SequencesComparer):
    # Identical characters are given 5 points, 4 point is deducted for each non-identical character
    # 3 points are deducted when opening a gap, and 0.1 points are deducted when extending it

    def compare(self, seq_a, seq_b):
        return pairwise2.align.globalms(seq_a, seq_b, 5, -4, -3, -0.1, score_only=True)


    def np_compare(self, i, combo):
        pass


class GlobalMsMin(SequencesComparer):
    # Identical characters are given 5 points, 4 point is deducted for each non-identical character
    # 3 points are deducted when opening a gap, and 0.1 points are deducted when extending it

    def compare(self, seq_a, seq_b):
        return -1 * pairwise2.align.globalms(seq_a, seq_b, 5, -4, -3, -0.1, score_only=True)


    def np_compare(self, i, combo):
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
                if GAP in pair:
                    gap = True
                    score += self.first_gap_score
                else:
                    score += self.__score_match(pair)
            else:
                if GAP not in pair:
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

    def np_compare(self, i, combo):
        pass


class MatchingCount(SequencesComparer):
    def __init__(self):
        self.opening_gap_penalty = 2
        self.continuation_gap_penalty = 3
        self.residue_match = -2
        self.error_penalty = 10
        self.mismatch_penalty = 0
        self.half_gap = self.opening_gap_penalty / 2.0


    def compare(self, seq_a, seq_b):
        score_total = 0
        gap_column = False

        for i in range(len(seq_a)):
            pos_a = seq_a[i]
            pos_b = seq_b[i]

            if pos_a == pos_b and pos_a == GAP:
                if gap_column:
                    score_local = self.continuation_gap_penalty
                else:
                    score_local = self.opening_gap_penalty

                gap_column = True
            elif pos_a == pos_b:
                gap_column = False
                score_local = self.residue_match
            elif pos_a == GAP or pos_b == GAP:
                gap_column = False
                score_local = self.half_gap
            elif pos_a != pos_b:
                score_local = self.mismatch_penalty
            else:
                gap_column = False
                score_local = self.error_penalty

            score_total += score_local

        return score_total


    def np_compare(self, i, combo):
        vec_a, vec_b = combo

        seq_a = [b for b in vec_a]
        seq_b = [b for b in vec_b]

        score_total = 0
        gap_column = False

        for i in range(len(seq_a)):
            pos_a = seq_a[i]
            pos_b = seq_b[i]

            if pos_a == pos_b and pos_a == B_GAP:
                if gap_column:
                    score_local = self.continuation_gap_penalty
                else:
                    score_local = self.opening_gap_penalty

                gap_column = True
            elif pos_a == pos_b:
                gap_column = False
                score_local = self.residue_match
            elif pos_a == B_GAP or pos_b == B_GAP:
                gap_column = False
                score_local = self.half_gap
            elif pos_a != pos_b:
                score_local = self.mismatch_penalty
            else:
                gap_column = False
                score_local = self.error_penalty

            score_total += score_local

        return score_total

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
        elif name == "single_ms":
            return SingleMS()
        else:
            raise NameError("Wrong name: {0}".format(name))