from numpy.random import rand
from numpy import ndarray
from math import floor
from input_parser import InputParser

from optimization_type import Maximization, Minimization
from results import Results

# https://machinelearningmastery.com/simulated-annealing-from-scratch-in-python/
class NpSimulatedAnnealing():


	def __init__(self, input_file, optimization):
		self.sequences_dictionary = InputParser.read_fasta_to_dict([input_file])

		self.initial: ndarray = InputParser.build_np_array(self.sequences_dictionary)

		if optimization == "max":
			self.optimization = Maximization()
		elif optimization == "min":
			self.optimization = Minimization()
		else:
			raise BaseException("Invalid optimization name")


	def execute(self, n_iterations: int, initial_temp: float, score_function, generate_neighbor, is_debugging=False):
		results: Results = Results()
		iterations_range: range = range(n_iterations)
		curr_temp: float = initial_temp

		if is_debugging:
			print(
				"i\tchanges\tcurrent_temp\tcurr_eval\tcandidate_eval\tdiff\tcriteria"
			)

		# sequences:
		rows, cols = self.initial.shape

		# generate an initial point
		curr = self.initial.copy()

		# evaluate the initial point
		curr_eval: float = score_function.np_calculate_score(curr)

		for i in iterations_range:
			curr_temp *= 1 - 0.003

			# take a step
			changes = max([1, floor(rows / (0.0025 * i + 1))])
			candidate = generate_neighbor(curr.copy(), i, changes=changes)

			# evaluate candidate point
			candidate_eval: float = score_function.np_calculate_score(candidate)

			# calculate temperature for current epoch
			metropolis_condition: bool = None

			# check if we should keep the new point:
			criteria = "\033[91mCANDIDATE IS NOT BETTER THAN CURRENT"

			# difference between candidate and current point evaluation
			original_curr_eval = curr_eval
			diff = candidate_eval - curr_eval

			optimization_condition = self.optimization.is_better_than_best(diff)
			if optimization_condition:
				criteria = "\033[93mCANDIDATE IS EQUAL TO CURRENT" if diff == 0 else "\033[92mCANDIDATE IS BETTER THAN CURRENT"

				# store the new current point
				curr = candidate
				curr_eval = candidate_eval
			else:
				metropolis_condition = self.optimization.metropolis(diff, curr_temp)
				random_value = rand()
				if random_value < metropolis_condition:
					criteria = "\033[95mMETROPOLIS CONDITION"
					# store the new current point
					# even if it is not better than the current one
					curr = candidate
					curr_eval = candidate_eval

			if is_debugging:
				print(
					"{color}{i:0>5d}\t{changes:0>5d}\t{current_temp:0>8f}\t{curr_eval:0>8f}\t{candidate_eval:0>8f}\t{diff:0>+8f}\t{criteria}\033[0m".format(
						color='\033[92m' if diff < 0 else '\033[93m' if diff == 0 else '\033[91m',
						i=i,
						changes=changes,
						current_temp=curr_temp,
						curr_eval=original_curr_eval,
						candidate_eval=candidate_eval,
						diff=diff,
						criteria=criteria
					)
				)

			results.register_iteration(i, changes, candidate_eval, curr_eval, metropolis_condition, curr_temp, diff)

		results.set_best(curr, curr_eval)
		return results
