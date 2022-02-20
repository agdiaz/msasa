import argparse

from runner import Runner
from simulated_annealing import SimulatedAnnealing

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Process some integers.')
	parser.add_argument('--input', dest='input_file', required=True, type=str, help='path to the input file')
	parser.add_argument('--output', dest='output_file', required=True, type=str, help='path to the output file')
	parser.add_argument('--comparer', dest='sequences_comparer', required=True, type=str, help='method to compare two sequences')
	parser.add_argument('--optimization', dest='optimization', required=False, default='min', type=str, choices=['min', 'max'], help='method to compare two sequences')
	parser.add_argument('--n-iterations', dest='n_iterations', type=int, required=False, default=50, help='number of iterations')
	parser.add_argument('--output-best-plot', dest='output_best_plot', required=False, default=None, help='path to the best results plot file', type=str)
	parser.add_argument('--output-temp-plot', dest='output_temp_plot', required=False, default=None, help='path to the temperatures plot file', type=str)
	parser.add_argument('--temperature', dest='initial_temp', required=False, default=10, type=int, help='initial temperature')
	args = parser.parse_args()

	Runner(args, SimulatedAnnealing).start()
