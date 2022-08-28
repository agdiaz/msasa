import numpy as np
from numpy.random import rand
from random import choice, choices, randint, getrandbits
from input_parser import InputParser
import re

B_GAP = b'-'

def msa_neighbor_add_remove(df, changes=1):
    seq_count = len(df.index)
    if changes > seq_count:
        changes = seq_count

    random_seq_indexes = set(choices(list(df.index), k=changes))

    sequences_dictionary = { 'sequences': {}, 'max_length': None, 'min_length': None }

    dft = df.transpose()
    for index, row in dft.items():
        values = dft[index]
        sequence = ''.join([x for x in values])

        if index in random_seq_indexes:
            if rand() > 0.5:
                # Add a GAP
                pos = randint(0, len(sequence) - 1)
                new_sequence = sequence[:pos] + "-" + sequence[pos:]
            else:
                # Remove a GAP
                gap_positions = [i.start() for i in re.finditer("-", sequence)]
                if len(gap_positions) > 0:
                    pos = choice(gap_positions)
                    new_sequence = sequence[0:pos:] + sequence[pos + 1::]
                else:
                    # No GAPS to remove
                    new_sequence = sequence

            sequence = new_sequence

        sequence_length = len(sequence)

        sequences_dictionary['sequences'][index] = { 'sequence': sequence, 'count': sequence_length }
        if sequences_dictionary['max_length'] == None or sequence_length > sequences_dictionary['max_length']:
            sequences_dictionary['max_length'] = sequence_length
        if sequences_dictionary['min_length'] == None or sequence_length < sequences_dictionary['min_length']:
            sequences_dictionary['min_length'] = sequence_length

    neighbor = InputParser.build_dataframe(sequences_dictionary)

    return neighbor


def __add_gap(array, index_to_alter = 0, pos = 0):

    _rows, cols = array.shape
    new_array = np.ndarray((0, cols + 1), dtype='|S1')
    l_trim = pos >= int((cols + 1) / 2)

    for row_index, row in enumerate(array):
        if row_index == index_to_alter:
            new_array = np.vstack([new_array, np.insert(row, pos, b'-')])
        elif l_trim:
            new_array = np.vstack([new_array, np.insert(row, 0, b'-')])
        else:
            new_array = np.vstack([new_array, np.insert(row, cols, b'-')])

    return new_array

def __remove_gap(array, index_to_alter = 0, pos = 0):
    _rows, cols = array.shape
    new_array = np.ndarray((0, cols), dtype='|S1')
    l_trim = pos < int((cols - 1) / 2)

    for row_index, row in enumerate(array):
        if row_index == index_to_alter:
            new_row = np.delete(row, pos)

            if l_trim:
                new_row = np.insert(new_row, 0, b'-')
            else:
                new_row = np.insert(new_row, cols - 1, b'-')

            new_array = np.vstack([new_array, new_row])
        else:
            new_array = np.vstack([new_array, row])

    return new_array

def np_msa_neighbor_add_remove(np_array, iteration: int, changes: int =1):
    sequences_count = len(np_array)
    range_of_sequences_count = range(sequences_count)

    new_np_array = np_array.copy()

    for _change_index in range(changes):
        neighbor = new_np_array.copy()

        while np.array_equal(new_np_array, neighbor):
            # print("Neighbor and current are still the same", iteration)
            attempts_to_safe_remove_gaps = 10
            remove_gap = False

            while attempts_to_safe_remove_gaps > 0:
                random_seq_index = choice(range_of_sequences_count)
                random_seq = np_array[random_seq_index]

                remove_gap = bool(getrandbits(1))
                if remove_gap:
                    gaps = [pos for pos, char in enumerate(random_seq) if char == B_GAP]
                    if len(gaps) > 0:
                        break
                    else:
                        attempts_to_safe_remove_gaps -= 1
                else:
                    break

            if remove_gap:
                position_to_edit = choice(gaps)
                neighbor = __remove_gap(np_array, random_seq_index, position_to_edit)
                # print("Changed applied removing one GAP", iteration, _change_index)
            else:
                position_to_edit = randint(0, len(random_seq) - 1)
                neighbor = __add_gap(np_array, random_seq_index, position_to_edit)
                # print("Changed applied adding one GAP", iteration, _change_index)

            while np.all(neighbor[:, -1] == B_GAP):
                neighbor = np.delete(neighbor, -1, 1)

            while np.all(neighbor[:, 0] == B_GAP):
                neighbor = np.delete(neighbor, 0, 1)

        # print("Neighbor and current are finally distinct!", iteration)
        new_np_array = neighbor

    return new_np_array
