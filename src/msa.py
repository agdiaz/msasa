import sys

from numpy.random import seed

from input_parser import InputParser
from simulated_annealing import SimulatedAnnealing

from objective_functions import SequencesComparerFactory
from generators import msa_neighbor_add_remove
from matplotlib import pyplot as plt
plt.ioff()


if __name__ == "__main__":
	input_file = sys.argv[1]
	output_file = sys.argv[2]
	sequences_comparer_name = sys.argv[3]
	n_iterations = int(sys.argv[4])
	score_function = SequencesComparerFactory.from_name(sequences_comparer_name)

	try:
		output_plot = sys.argv[5]
		# from matplotlib import pyplot as plt
	except IndexError:
		output_plot = None
		print("NO PLOT")

	sequences_dictionary = InputParser.read_fasta_to_dict([sys.argv[1]])
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
	best, score, bests, currents, candidates, temperatures = sa.maximize(score_function, msa_neighbor_add_remove)

	print('%s;%i;%i' % (sys.argv[2], initial_energy, score))

	InputParser.dataframe_to_msa_file(best, sys.argv[2])

	with open(sys.argv[2], "r") as f:
		content = f.read()

	if output_plot != None:
		print(bests)
		print(candidates)
		fig, ax = plt.subplots(figsize=(25, 6))  # Create a figure containing a single axes.
		plt.title("Simulated Annealing - MSA prediction ({0})".format(sequences_comparer_name))
		ax.plot(range(n_iterations), bests, color="green")  # Plot some data on the axes.
		ax.plot(range(n_iterations), candidates, color="orange")  # Plot some data on the axes.

		axTemp = ax.twinx()
		axTemp.plot(range(n_iterations), temperatures, color='red')

		ax.set_xlabel('Iterations')
		ax.set_ylabel('Score', color='g')
		axTemp.set_ylabel('Temperature', color='b')

		fig.savefig(sys.argv[3])
		# plt.close()
