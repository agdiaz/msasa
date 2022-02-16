from Bio import SeqIO
# from Bio.SubsMat import MatrixInfo
# from Bio.SubsMat.MatrixInfo import blosum62 as blosum
from Bio import pairwise2
# from numpy import asarray
from numpy import exp
# from numpy.random import randn
from numpy.random import rand
from numpy.random import seed
from math import floor
from random import choices, randint, choice
from itertools import combinations
import re
# from matplotlib import pyplot as plt
import pandas as pd
import sys

class InputParser():
    def read_fasta_to_bioseq(self, path_to_files):
        sequences = []

        for file_path in path_to_files:
            for record in SeqIO.parse(file_path, "fasta"):
                sequences.append(record.seq)

        return sequences

    def read_fasta_to_dict(self, path_to_files):
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

    def build_dataframe(self, sequences_dictionary):
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

    def dataframe_to_sequences(self, df):
        sequences = []
        for index, row in df.transpose().items():
            values = df.transpose()[index]
            sequence = ''.join([x for x in values])
            sequences.append(sequence)
        return sequences

    def dataframe_to_msa_file(self, df, file_name):
        last_col_index = len(df.columns) - 1
        while (df[last_col_index] == "-").all():
            df.drop(df.columns[last_col_index], axis=1, inplace=True)
            last_col_index -= 1

        with open(file_name, "w") as output_file:
            for index, row in df.transpose().items():
                values = df.transpose()[index]
                sequence = ''.join([x for x in values])

                output_file.write(">{0}\n{1}\n".format(index, sequence))

def simulated_annealing_max(objective, neighbor, initial, n_iterations, temp):
	# for the records
	candidates = []
	currents = []
	bests = []
	temperatures = []

	# generate an initial point
	best = initial
	# evaluate the initial point
	best_eval = objective(best)
	# current working solution
	curr, curr_eval = best, best_eval

	# run the algorithm
	for i in range(n_iterations):
		changes = floor(n_iterations/(i + 1))

		# take a step
		candidate = neighbor(best, changes=changes)
		# evaluate candidate point
		candidate_eval = objective(candidate)

		if i % 10 == 0:
			print("{0},{1},{2},{3},{4}".format(i, best_eval, curr_eval, candidate_eval, changes))

		# check for new best solution
		if candidate_eval > best_eval:
			# store new best point
			best, best_eval = candidate, candidate_eval

		# difference between candidate and current point evaluation
		diff = candidate_eval - curr_eval

		# calculate temperature for current epoch
		t = temp / float(i + 1)

		# calculate metropolis acceptance criterion
		metropolis = exp(diff / t)

		# check if we should keep the new point
		if diff > 0 or rand() < metropolis:
			# store the new current point
			curr, curr_eval = candidate, candidate_eval

		bests.append(best_eval)
		temperatures.append(t)
		currents.append(curr_eval)
		candidates.append(candidate_eval)

	return [best, best_eval, bests, currents, candidates, temperatures]

def msa_objective_real(df):
    total_energy = 0
    sequences = InputParser().dataframe_to_sequences(df)
    seq_combinations = combinations(sequences, 2)

    for combo in seq_combinations:
        seq_a, seq_b = combo
        # Identical characters are given 5 points, 4 point is deducted for each non-identical character
        # 3 points are deducted when opening a gap, and 0.1 points are deducted when extending it
        blosum_score = pairwise2.align.globalms(seq_a, seq_b, 5, -4, -3, -0.1, score_only=True)
        total_energy += blosum_score

    return total_energy

def msa_neighbor_add_remove(df, changes=1):
    seq_count = len(df.index)
    if changes > seq_count:
        changes = seq_count

    random_seq_indexes = choices(list(df.index), k=changes)

    with open("/tmp/msa_neighbor_add_remove.fasta", "w") as tmp_file:
        dft = df.transpose()
        for index, row in dft.items():
            values = dft[index]
            sequence = ''.join([x for x in values])

            if index in random_seq_indexes:
                if randint(0, 1) == 0:
                    pos = randint(0, len(sequence) - 1)
                    new_sequence = "".join((sequence[:pos], "-", sequence[pos:]))
                else:
                    gap_positions = [i.start() for i in re.finditer("-", sequence)]
                    if len(gap_positions) > 0:
                        pos = choice(gap_positions)
                        new_sequence = sequence[0:pos:] + sequence[pos + 1::]
                    else:
                        new_sequence = sequence

                sequence = new_sequence

            tmp_file.write(">{0}\n{1}\n".format(index, sequence))

    sequences_dictionary = input_parser.read_fasta_to_dict(["/tmp/msa_neighbor_add_remove.fasta"])
    neighbor = input_parser.build_dataframe(sequences_dictionary)

    return neighbor

if __name__ == "__main__":
	print("WELCOME")
	input_parser = InputParser()
	sequences_dictionary = input_parser.read_fasta_to_dict([sys.argv[1]])
	df = input_parser.build_dataframe(sequences_dictionary)

	# seed the pseudorandom number generator
	seed(2)
	# Init df
	initial_df = df
	initial_energy = msa_objective_real(initial_df)

	# define the total iterations
	n_iterations = 50
	# initial temperature
	temp = 10
	# perform the simulated annealing search
	best, score, bests, currents, candidates, temperatures = simulated_annealing_max(msa_objective_real, msa_neighbor_add_remove, initial_df, n_iterations, temp)
	print('Initial Score = %i; Final Score = %i' % (initial_energy, score))

	input_parser = InputParser()
	input_parser.dataframe_to_msa_file(best, sys.argv[2])

	with open(sys.argv[2], "r") as f:
		content = f.read()

	print(content)

	# fig, ax = plt.subplots(figsize=(25, 6))  # Create a figure containing a single axes.
	# ax.plot(bests, color="green")  # Plot some data on the axes.
	# ax.plot(candidates, color="orange")  # Plot some data on the axes.
	# plt.show()