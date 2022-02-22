import numpy as np
from numpy.random import rand
from random import choice, choices, randint
from input_parser import InputParser
import re


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

def __create_neighbor_removing_gap(sequence):
    gaps = [pos for pos, char in enumerate(sequence) if char == '-']
    pos = choice(gaps)

    return sequence[0:pos:] + sequence[pos + 1::]


def __create_neighbor_adding_gap(sequence, changes):
    pos = randint(0, len(sequence) - 1)

    return sequence[:pos] + "-" + sequence[pos:]

def __is_last_column_gap(sequences, max_sequence_length):
    return all(s[max_sequence_length - 1] == '-' for s in sequences)

def __trim_sequences(sequences):
    for seq_index in range(len(sequences)):
        sequences[seq_index] = sequences[seq_index][:-1]

    return sequences

def np_msa_neighbor_add_remove(np_array, changes=1):
    to_str = lambda vector: "".join(vector.decode("utf-8"))
    sequences = [to_str(s) for s in np_array]
    sequences_count = len(sequences)
    range_of_sequences_count = range(sequences_count)

    while True:
        random_seq_index = choice(range_of_sequences_count)
        random_seq = sequences[random_seq_index]

        remove_gap = randint(0, 1) == 0
        if not remove_gap:
            break
        else:
            if random_seq.count('-') > 0:
                break

    if remove_gap:
        new_random_seq = __create_neighbor_removing_gap(random_seq)
    else:
        new_random_seq = __create_neighbor_adding_gap(random_seq, changes)

    sequences[random_seq_index] = new_random_seq

    max_sequence_length = len(max(sequences, key = len))
    for seq_index in range_of_sequences_count:
        sequences[seq_index] = sequences[seq_index].ljust(max_sequence_length, '-')

    while(__is_last_column_gap(sequences, max_sequence_length)):
        sequences = __trim_sequences(sequences)
        max_sequence_length -= 1

    neighbor = np.chararray([sequences_count, max_sequence_length])
    for seq_index, sequence in enumerate(sequences):
        neighbor[seq_index] = [c for c in sequence]

    return neighbor
