from abc import abstractmethod, ABCMeta
from functools import cached_property, lru_cache, partial
from itertools import combinations
# from Bio.SubsMat.MatrixInfo import blosum62 as blosum
# from input_parser import InputParser
# import numpy as np
from numpy import array, apply_along_axis, unique, stack, sort
import numexpr as ne

from r_blosum69 import B_GAP, WEIGHTS, tuple_score

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

    @lru_cache(maxsize=None)
    def normalized_column_size(self, columns):
        normalized_count_of_columns = 1 + (columns - self.original_columns) / ((2 * self.original_columns) - self.original_columns)
        return normalized_count_of_columns

    @cached_property
    def original_columns(self):
        return self.alignment_shape[1]

    @lru_cache(maxsize=None)
    def group_by_residues(self, column_as_tuple):
        column = array(column_as_tuple)
        return unique(column, return_counts=True)


class SingleMS(SequencesComparer):
    def np_calculate_score(self, alignment_ndarray):
        columns = alignment_ndarray.shape[1]

        result_objects = apply_along_axis(self.np_compare, axis=0, arr=alignment_ndarray)
        column_energies = ne.evaluate("sum(result_objects)", local_dict = {'result_objects': result_objects})

        return self.normalized_column_size(columns) * column_energies

    def local_energy(self, column_length, elem):
        residue_type, residue_count = elem

        return self.proportioned_energy(column_length, residue_type, residue_count)

    @lru_cache(maxsize=None)
    def proportioned_energy(self, column_length, residue_type, residue_count):
        proportion = int(residue_count) / column_length

        return proportion * WEIGHTS[residue_type]


    def np_compare(self, column):
        return self.cached_np_compare(tuple(sort(column)))

    @lru_cache(maxsize=None)
    def cached_np_compare(self, column_as_tuple):
        column_length = len(column_as_tuple)
        comparison_function = partial(self.local_energy, column_length)

        unique_residues, count_residues = self.group_by_residues(column_as_tuple)
        grouped_residues = stack((unique_residues, count_residues), axis = 1)

        local_scores = apply_along_axis(comparison_function, axis=1, arr=grouped_residues)

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
    def np_calculate_score(self, alignment_ndarray):
        columns = alignment_ndarray.shape[1]

        result_objects = apply_along_axis(self.np_compare, axis=0, arr=alignment_ndarray)
        total_score = ne.evaluate("sum(result_objects)", local_dict={'result_objects': result_objects})

        return self.normalized_column_size(columns) * -1 * total_score

    def np_compare(self, column):
        return self.cached_blosum_score(tuple(sort(column)))

    @lru_cache(maxsize=None)
    def cached_blosum_score(self, column_as_tuple):
        combinations_count, uniques, counts = self.cached_combinations(column_as_tuple)
        blosum_scores = apply_along_axis(self.score, axis=1, arr=uniques)

        pondered_blosum_sum = ne.evaluate("sum(blosum_scores * (counts / combinations_count))", local_dict={
            'blosum_scores': blosum_scores,
            'counts': counts,
            'combinations_count': combinations_count
        })

        return pondered_blosum_sum

    @lru_cache(maxsize=None)
    def cached_combinations(self, column_as_tuple):
        combination_of_residues = array(list(combinations(column_as_tuple, 2)))
        combinations_count      = len(combination_of_residues)
        uniques, counts         = unique(combination_of_residues, return_counts=True, axis=0)

        return (combinations_count, uniques, counts)

    def score(self, residues_tuple):
        return tuple_score(residues_tuple[0], residues_tuple[1])


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