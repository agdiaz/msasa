# import numpy as np
from numpy import array, all, delete, insert, array_equal, nonzero
from random import choice, randint, random
from r_blosum69 import B_GAP
from functools import lru_cache, partial

class AlignmentGenerator:
    def realign_sequences(self, alignment, changes: int = 1):
        sequences_count = len(alignment)
        range_of_sequences_count = range(sequences_count)

        new_alignment = alignment

        for _change_index in range(changes):
            neighbor_alignment = new_alignment.copy()

            while array_equal(new_alignment, neighbor_alignment):
                random_seq_index = choice(range_of_sequences_count)
                random_seq = neighbor_alignment[random_seq_index]
                has_changed = False

                if random() < 0.5:
                    position_to_edit = randint(0, len(random_seq) - 1)
                    modifier = partial(self.__modify_row_adding_gap, random_seq_index, position_to_edit)
                    neighbor_alignment = self.__add_gap(neighbor_alignment, modifier)
                    has_changed = True
                else:
                    gaps = nonzero(random_seq == B_GAP)[0]
                    if len(gaps) > 0:
                        position_to_edit = choice(gaps)
                        modifier = partial(self.__modify_row_removing_gap, random_seq_index, position_to_edit)
                        neighbor_alignment = self.__remove_gap(neighbor_alignment, modifier)
                        has_changed = True

                if has_changed:
                    neighbor_alignment = self.clean_alignment_gaps(tuple(neighbor_alignment.flatten()), neighbor_alignment.shape) # array([column for column in neighbor_alignment.T if not all(column == B_GAP)]).T

            new_alignment = neighbor_alignment

        return new_alignment

    @lru_cache(maxsize=None)
    def clean_alignment_gaps(self, neighbor_alignment_as_tuple: tuple, original_shape: tuple):
        neighbor_alignment = array(neighbor_alignment_as_tuple).reshape(original_shape)

        return array([column for column in neighbor_alignment.T if not all(column == B_GAP)]).T

    def __add_gap(self, alignment, modifier):
        row_list = [modifier(row_index, row) for row_index, row in enumerate(alignment)]

        return array(row_list)

    def __remove_gap(self, alignment, modifier):
        row_list = [modifier(row_index, row) for row_index, row in enumerate(alignment)]

        return array(row_list)

    def __modify_row_adding_gap(self, index_to_alter, pos, row_index, row):
        insert_position = pos if row_index == index_to_alter else 0 if random() < 0.5 else -1

        return insert(row, insert_position, B_GAP)

    def __modify_row_removing_gap(self, index_to_alter, pos, row_index, row):
        if row_index == index_to_alter:
            insert_position = 0 if random() < 0.5 else -1

            return insert(delete(row, pos), insert_position, B_GAP)
        else:
            return row
