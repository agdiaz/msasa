import numpy as np
from numpy.random import rand
from random import choice, choices, randint, getrandbits
from input_parser import InputParser
import re

B_GAP = 45

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

def __is_last_column_gap(sequences, max_sequence_length):
    return all(s[max_sequence_length - 1] == B_GAP for s in sequences)

def __is_first_column_gap(sequences, max_sequence_length):
    return all(s[0] == B_GAP for s in sequences)

def __ltrim_sequences(sequences):
    for index, elem in enumerate(sequences):
        sequences[index] = elem[:-1]

    # return sequences

def __rtrim_sequences(sequences):
    for index, elem in enumerate(sequences):
        sequences[index] = elem[1:]

    # return sequences

def __add_gap(array, index_to_alter = 0, pos = 0):
    sequence_to_alter = array[index_to_alter]

    new = sequence_to_alter[:pos] + b'-' + sequence_to_alter[pos:]
    new_max_length = len(new)

    if pos >= int(new_max_length / 2):
        adjusted_sequences = np.char.ljust(array, new_max_length, fillchar="-")
    else:
        adjusted_sequences = np.char.rjust(array, new_max_length, fillchar="-")

    adjusted_sequences[index_to_alter] = new

    return adjusted_sequences

def __remove_gap(array, index_to_alter = 0, pos = 0):
    sequence_to_alter = array[index_to_alter]

    new = sequence_to_alter[0 : pos : ] + sequence_to_alter[pos + 1 : :]
    new_max_length = max(len(new), len(max(array, key=len)))

    if pos < int(new_max_length / 2):
        adjusted_sequences = np.char.ljust(array, new_max_length, fillchar="-")
        new = np.char.ljust(new, new_max_length, fillchar="-")
    else:
        adjusted_sequences = np.char.rjust(array, new_max_length, fillchar="-")
        new = np.char.rjust(new, new_max_length, fillchar="-")

    adjusted_sequences[index_to_alter] = new

    return adjusted_sequences

def np_msa_neighbor_add_remove(np_array, changes=1):
    sequences_count = len(np_array)
    range_of_sequences_count = range(sequences_count)

    while True:
        random_seq_index = choice(range_of_sequences_count)
        random_seq = np_array[random_seq_index]

        remove_gap = bool(getrandbits(1))
        if remove_gap:
            gaps = [pos for pos, char in enumerate(random_seq) if char == B_GAP]
            if len(gaps) > 0:
                break
        else:
            break

    if remove_gap:
        pos = choice([pos for pos, char in enumerate(random_seq) if char == B_GAP])
        neighbor = __remove_gap(np_array, random_seq_index, pos)
    else:
        pos = randint(0, len(random_seq) - 1)
        neighbor = __add_gap(np_array, random_seq_index, pos)

    max_sequence_length = len(max(neighbor, key=len))
    while __is_last_column_gap(neighbor, max_sequence_length):
        __ltrim_sequences(neighbor)
        max_sequence_length -= 1

    while __is_first_column_gap(neighbor, max_sequence_length):
        __rtrim_sequences(neighbor)
        max_sequence_length -= 1

    return neighbor
