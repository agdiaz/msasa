import sys
import argparse
import pathlib
from numpy.random import seed

from input_parser import InputParser
from simulated_annealing import SimulatedAnnealing, Results

from objective_functions import SequencesComparerFactory
from generators import msa_neighbor_add_remove
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
plt.ioff()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Process some integers.')
	parser.add_argument('--input', dest='input_file', required=True, type=str, help='path to the input file')
	parser.add_argument('--output', dest='output_file', required=True, type=str, help='path to the output file')
	parser.add_argument('--comparer', dest='sequences_comparer', required=True, type=str, help='method to compare two sequences')
	parser.add_argument('--optimization', dest='optimization_type', required=False, default='min', type=str, choices=['min', 'max'], help='method to compare two sequences')
	parser.add_argument('--n-iterations', dest='n_iterations', type=int, required=False, default=50, help='number of iterations')
	parser.add_argument('--output-plot', dest='output_plot', required=False, help='path to the results plot file', type=str)
	args = parser.parse_args()

	input_file = args.input_file
	output_file = args.output_file
	sequences_comparer_name = args.sequences_comparer
	n_iterations = args.n_iterations
	score_function = SequencesComparerFactory.from_name(sequences_comparer_name)

	try:
		output_plot = args.output_plot
	except IndexError:
		output_plot = None

	sequences_dictionary = InputParser.read_fasta_to_dict([input_file])
	alignment_dataframe = InputParser.build_dataframe(sequences_dictionary)

	# seed the pseudorandom number generator
	seed(2)

	# Init df
	initial_alignment_dataframe = alignment_dataframe
	initial_energy = score_function.calculate_score(alignment_dataframe)

	# initial temperature
	temp = 10

	# perform the simulated annealing search
	sa = SimulatedAnnealing(initial_alignment_dataframe, n_iterations, temp)

	if args.optimization_type == 'max':
		results = sa.maximize(score_function, msa_neighbor_add_remove)
	else:
		results = sa.minimize(score_function, msa_neighbor_add_remove)

	best, score = results.best()

	print('%s;%i;%i' % (output_file, initial_energy, score))

	InputParser.dataframe_to_msa_file(best, output_file)

	if output_plot != None:
		fig, ax = plt.subplots(figsize=(25, 10))
		plt.title("Simulated Annealing - MSA prediction ({0})".format(sequences_comparer_name))
		ax.plot(results.records("bests"), color="green")
		ax.plot(results.records("candidates"), color="orange")

		axTemp = ax.twinx()
		axTemp.plot(results.records("temperatures"), color='red')

		ax.set_xlabel('Iterations')
		ax.set_ylabel('Score', color='g')
		axTemp.set_ylabel('Temperature', color='b')

		fig.savefig(output_plot)
