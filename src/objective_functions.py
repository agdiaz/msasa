from abc import abstractmethod, ABCMeta
from itertools import combinations
# from Bio.SubsMat.MatrixInfo import blosum62 as blosum
# from input_parser import InputParser
# import numpy as np
from numpy import array, apply_along_axis, unique, stack
import numexpr as ne

from r_blosum69 import B_GAP, tuple_score

# blosum.update(((b,a),val) for (a,b),val in list(blosum.items()))

# https://github.com/biotite-dev/biotite/blob/master/src/biotite/sequence/align/matrix_data/RBLOSUM69_13p.mat
class SequencesComparer(metaclass=ABCMeta):
    def __init__(self, alignment_shape):
        self.alignment_shape = alignment_shape

    @abstractmethod
    def np_calculate_score(self, alignment_ndarray):
        pass

    @abstractmethod
    def np_compare(self, i, combo):
        pass


class SingleMS(SequencesComparer):

    def np_calculate_score(self, alignment_ndarray):
        original_rows, original_columns = self.alignment_shape
        _rows, columns = alignment_ndarray.shape

        normalized_count_of_columns = 1 + (columns - original_columns) / ((2 * original_columns) - original_columns)

        result_objects = apply_along_axis(self.np_compare, axis=0, arr=alignment_ndarray)
        column_energies = ne.evaluate("sum(result_objects)", local_dict = {'result_objects': result_objects})

        return normalized_count_of_columns * column_energies

    def np_compare(self, column):
        column_length = len(column)
        GAP_WEIGHT = 0.8
        RES_WEIGHT = 0.2

        def local_energy(elem):
            residue_type, count = elem
            proportion = int(count) / column_length
            return proportion * GAP_WEIGHT if residue_type == B_GAP else proportion * RES_WEIGHT

        unique_residues, count_residues = unique(column, return_counts=True)
        grouped_residues = stack((unique_residues, count_residues), axis = 1)

        local_scores = apply_along_axis(local_energy, axis=1, arr=grouped_residues)

        return ne.evaluate("sum(local_scores)", local_dict={'local_scores': local_scores})


class SingleMatching(SequencesComparer):
    # Groups of identical characters are given 1 points * size of the group,
    # 5 points are deducted for a each gap

    def np_calculate_score(self, alignment_ndarray):
        columns = alignment_ndarray.shape[1]
        result_objects = apply_along_axis(self.np_compare, axis=0, arr=alignment_ndarray)

        return columns * len(result_objects) * ne.evaluate("sum(result_objects)", local_dict={"result_objects": result_objects})

    def np_compare(self, column):
        unique_residues = unique(column)
        unique_residues_count = len(unique_residues)

        # ALL GAPS IN A COLUMN:
        if unique_residues_count == 1 and unique_residues[0] == B_GAP:
            return 1
        else:
            return unique_residues_count


class SingleBlosum(SequencesComparer):
    # Groups of identical characters are given 1 points * size of the group,
    # 5 points are deducted for a each gap

    def np_calculate_score(self, alignment_ndarray):
        _original_rows, original_columns = self.alignment_shape
        _rows, columns = alignment_ndarray.shape
        normalized_count_of_columns = 1 + (columns - original_columns) / ((2 * original_columns) - original_columns)

        result_objects = apply_along_axis(self.np_compare, axis=0, arr=alignment_ndarray)
        total_score = ne.evaluate("sum(result_objects)", local_dict={'result_objects': result_objects})

        return normalized_count_of_columns * -1 * total_score

    def np_compare(self, column):
        combination_of_residues = array(list(combinations(column, 2)))
        combinations_count = len(combination_of_residues)

        uniques, counts = unique(combination_of_residues, return_counts=True, axis=0)
        blosum_scores = apply_along_axis(self.score, axis=1, arr=uniques)
        pondered_blosum_sum = ne.evaluate("sum(blosum_scores * counts / combinations_count)", local_dict={'blosum_scores': blosum_scores, 'counts': counts, 'combinations_count': combinations_count})

        return float(pondered_blosum_sum)

    def score(self, residues_tuple):
        residue_a, residue_b = residues_tuple
        return tuple_score(residue_a, residue_b)


class SequencesComparerFactory:
    @staticmethod
    def from_name(name, alignment_shape):
        if name == "single_ms":
            return SingleMS(alignment_shape)
        elif name == "single_blosum":
            return SingleBlosum(alignment_shape)
        elif name == "single_matching":
            return SingleMatching(alignment_shape)
        else:
            raise NameError("Wrong name: {0}".format(name))