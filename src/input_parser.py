from Bio import SeqIO
import pandas as pd
import numpy as np
from random import random

class InputParser():
    @staticmethod
    def read_fasta_to_bioseq(path_to_files):
        sequences = []

        for file_path in path_to_files:
            for record in SeqIO.parse(file_path, "fasta"):
                sequences.append(record.seq)

        return sequences

    @staticmethod
    def read_fasta_to_dict(path_to_files):
        sequences_dictionary = { 'sequences': {}, 'max_length': None, 'min_length': None }
        max_length = None
        min_length = None

        for file_path in path_to_files:
            for record in SeqIO.parse(file_path, "fasta"):
                sequence_length = record.seq.__len__()
                sequences_dictionary['sequences'][record.id] = { 'sequence': record.seq.__str__(), 'count': sequence_length }

                if max_length == None or sequence_length > max_length:
                    max_length = sequence_length
                if min_length == None or sequence_length < min_length:
                    min_length = sequence_length

        sequences_dictionary['max_length'] = max_length
        sequences_dictionary['min_length'] = min_length

        return sequences_dictionary

    @staticmethod
    def sequences_to_dict(sequences):
        sequences_dictionary = { 'sequences': {}, 'max_length': None, 'min_length': None }
        max_length = None
        min_length = None

        for record in sequences:
            sequence_length = record.seq.__len__()
            sequences_dictionary['sequences'][record.id] = { 'sequence': record.seq.__str__(), 'count': sequence_length }

            if max_length == None or sequence_length > max_length:
                max_length = sequence_length
            if min_length == None or sequence_length < min_length:
                min_length = sequence_length

        sequences_dictionary['max_length'] = max_length
        sequences_dictionary['min_length'] = min_length

        return sequences_dictionary

    @staticmethod
    def build_np_array(sequences_dictionary):
        sequences = [seq_dict['sequence'] for seq_dict in sequences_dictionary['sequences'].values()]
        max_length = sequences_dictionary['max_length']

        if (random() > 0.5):
            adjusted_sequences = np.char.ljust(sequences, max_length, fillchar="-")
        else:
            adjusted_sequences = np.char.rjust(sequences, max_length, fillchar="-")

        x = np.array(adjusted_sequences, dtype=bytes)
        return x.view('S1').reshape((x.size, -1))


    @staticmethod
    def build_dataframe(sequences_dictionary):
        max_length = sequences_dictionary['max_length']
        index = range(max_length)
        columns = {}

        for (entry, value) in sequences_dictionary['sequences'].items():
            residues = value['sequence']
            if len(residues) < max_length:
                residues = residues.ljust(max_length, '-')

            columns[entry] = [c for c in residues]

        df = pd.DataFrame.from_dict(columns, orient='index', columns=index)

        last_col_index = max_length - 1
        while (df[last_col_index] == "-").all():
            df.drop(df.columns[last_col_index], axis=1, inplace=True)
            last_col_index -= 1

        return df

    @staticmethod
    def dataframe_to_sequences(df):
        sequences = []
        for index, row in df.transpose().items():
            values = df.transpose()[index]
            sequence = ''.join([x for x in values])
            sequences.append(sequence)
        return sequences

    @staticmethod
    def np_array_to_sequences(np_array):
        seqs, max_length = np_array.shape

        return np_array.view("S%i" % max_length)

    @staticmethod
    def dataframe_to_msa_file(df, file_name):
        last_col_index = len(df.columns) - 1
        while (df[last_col_index] == "-").all():
            df.drop(df.columns[last_col_index], axis=1, inplace=True)
            last_col_index -= 1

        with open(file_name, "w") as output_file:
            for index, row in df.transpose().items():
                values = df.transpose()[index]
                sequence = ''.join([x for x in values])

                output_file.write(">{0}\n{1}\n".format(index, sequence))

    def np_array_to_msa_file(np_array, sequences_dict, file_name):
        sequences = []
        for seq_ind, row in enumerate(np_array):
            aligned_sequence = "".join([column.decode() for column in row])

            sequences.append(aligned_sequence)

        sequence_names = [*sequences_dict['sequences'].keys()]

        with open(file_name, "w") as output_file:
            for seq_ind, sequence in enumerate(sequences):
                output_file.write(">{0}\n{1}\n".format(sequence_names[seq_ind], sequence))
