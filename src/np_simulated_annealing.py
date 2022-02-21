from math import floor
from numpy.random import rand
from numpy import ndarray
import pandas as pd

from input_parser import InputParser

from optimization_type import Maximization, Minimization
from results import Results

# https://machinelearningmastery.com/simulated-annealing-from-scratch-in-python/
class NpSimulatedAnnealing():


	def __init__(self, input_file, optimization) -> None:
		self.sequences_dictionary = InputParser.read_fasta_to_dict([input_file])

		self.initial: ndarray = InputParser.build_np_array(self.sequences_dictionary)

		if optimization == "max":
			self.optimization = Maximization()
		elif optimization == "min":
			self.optimization = Minimization()
		else:
			raise BaseException("Invalid optimization name")

	def execute(self, n_iterations: int, initial_temp: int, score_function, generate_neighbor) -> Results:
		results: Results = Results()
		iterations_range: range = range(n_iterations)

		# generate an initial point
		best = generate_neighbor(self.initial, changes=1)

		# evaluate the initial point
		best_eval = score_function.np_calculate_score(best)

		# current working solution
		curr = best
		curr_eval = best_eval

		for i in iterations_range:
			changes: int = floor(n_iterations/(i + 1))

			# take a step
			candidate: pd.DataFrame = generate_neighbor(curr, changes=changes)
			for j in range(10):
				candidate: pd.DataFrame = generate_neighbor(candidate, changes=changes)

			# evaluate candidate point
			candidate_eval: float = score_function.np_calculate_score(candidate)

			# difference between candidate and current point evaluation
			diff = candidate_eval - best_eval
			optimization_condition = self.optimization.is_better_than_best(diff)

			# store new best point
			if optimization_condition:
				best = candidate
				best_eval = candidate_eval

			# calculate temperature for current epoch
			current_temp = float(initial_temp / float(i + 1))
			metropolis_condition = self.optimization.metropolis(diff, current_temp)

			# check if we should keep the new point
			if optimization_condition or rand() < metropolis_condition:
				# store the new current point
				curr = candidate
				curr_eval = candidate_eval

			results.register_iteration(i, candidate_eval, curr_eval, best_eval, metropolis_condition, current_temp, diff)

		results.set_best(best, best_eval)
		return results
