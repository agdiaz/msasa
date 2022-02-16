import sys

from numpy.random import seed

from input_parser import InputParser
from simulated_annealing import SimulatedAnnealing

from objective_functions import msa_objective_real
from generators import msa_neighbor_add_remove

if __name__ == "__main__":
	input_file = sys.argv[1]
	output_file = sys.argv[2]
	try:
		output_plot = sys.argv[3]
		from matplotlib import pyplot as plt
	except IndexError:
		output_plot = None

	sequences_dictionary = InputParser.read_fasta_to_dict([sys.argv[1]])
	df = InputParser.build_dataframe(sequences_dictionary)

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
	sa = SimulatedAnnealing(initial_df, n_iterations, temp)
	best, score, bests, currents, candidates, temperatures = sa.maximize(msa_objective_real, msa_neighbor_add_remove)

	print('%s;%i;%i' % (sys.argv[2], initial_energy, score))

	InputParser.dataframe_to_msa_file(best, sys.argv[2])

	with open(sys.argv[2], "r") as f:
		content = f.read()

	if output_plot != None:
		fig, ax = plt.subplots(figsize=(25, 6))  # Create a figure containing a single axes.

		ax.plot(bests, color="green")  # Plot some data on the axes.
		ax.plot(candidates, color="orange")  # Plot some data on the axes.

		axTemp = ax.twinx()
		axTemp.plot(temperatures, color='red')

		ax.set_xlabel('Iteration')
		ax.set_ylabel('Score', color='g')
		axTemp.set_ylabel('Temperature', color='b')

		fig.savefig(sys.argv[3])
		plt.close()
