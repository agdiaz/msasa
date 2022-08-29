import numpy as np
from input_parser import InputParser
from random import randint

from optimization_type import Maximization, Minimization, metro
from results import Results

# https://machinelearningmastery.com/simulated-annealing-from-scratch-in-python/
class NpSimulatedAnnealing():
	CRITERIA_BETTER = "\033[91mCANDIDATE IS BETTER THAN CURRENT"
	CRITERIA_WORSE = "\033[91mCANDIDATE IS NOT BETTER THAN CURRENT"
	CRITERIA_EQUAL = "\033[93mCANDIDATE IS EQUAL TO CURRENT"
	CRITERIA_METROPOLIS = "\033[95mMETROPOLIS CONDITION"
	CRITERIA_OMITTED = ""

	def __init__(self, input_file, optimization):
		self.sequences_dictionary = InputParser.read_fasta_to_dict([input_file])

		self.initial = InputParser.build_np_array(self.sequences_dictionary)

		if optimization == "max":
			self.optimization = Maximization()
		elif optimization == "min":
			self.optimization = Minimization()
		else:
			raise BaseException("Invalid optimization name")


	def execute(self, n_iterations: int, initial_temp: float, score_function, generate_neighbor, is_debugging=False):
		results: Results = Results()

		if is_debugging:
			print(
				"i\tavail_changes\tchanges\tcurrent_temp\tcurr_eval\tcandidate_eval\tdiff\tcriteria"
			)

		# sequences:
		sequences_count = self.initial.shape[0]

		iterations_array = np.arange(n_iterations, dtype=int)
		changes_array = np.floor(sequences_count / ((0.0025 * iterations_array) + 1))
		temperatures_array = initial_temp / ((0.5 * iterations_array) + 1)
		randoms = np.random.rand(n_iterations)
		iterations = np.stack((iterations_array, changes_array, temperatures_array, randoms), axis = 1)

		# generate an initial point and evaluate the initial point
		curr = self.initial
		curr_eval: float = score_function.np_calculate_score(curr)

		results.set_initial(curr, curr_eval)

		for iteration_index, available_changes, curr_temp, iteration_random in iterations:
			available_changes_int = int(available_changes)
			changes = 1 if available_changes_int < 2 else randint(1, available_changes_int)

			# evaluate candidate point
			candidate      = generate_neighbor(curr, changes)
			candidate_eval = score_function.np_calculate_score(candidate)

			# check if we should keep the new point:
			criteria = NpSimulatedAnnealing.CRITERIA_WORSE

			# difference between candidate and current point evaluation
			original_curr_eval = curr_eval
			diff               = candidate_eval - curr_eval

			optimization_condition = self.optimization.is_better_than_best(diff)
			if optimization_condition:
				criteria = NpSimulatedAnnealing.CRITERIA_EQUAL if diff == 0 else NpSimulatedAnnealing.CRITERIA_BETTER

				# store the new current point
				if diff != 0:
					curr, curr_eval = candidate, candidate_eval
			elif metro(curr_eval, candidate_eval, curr_temp, iteration_random): # iteration_random < self.optimization.metropolis(diff, curr_temp):
				# store the new current point even if it is not better than the current one
				criteria = NpSimulatedAnnealing.CRITERIA_METROPOLIS
				curr, curr_eval  = candidate, candidate_eval

			if is_debugging:
				print(
					"{color}{i:0>5d}\t{available_changes:0>5d}\t{changes:0>5d}\t{current_temp:0>8f}\t{curr_eval:0>8f}\t{candidate_eval:0>8f}\t{diff:0>+8f}\t{criteria}\033[0m".format(
						color='\033[92m' if diff < 0 else '\033[93m' if diff == 0 else '\033[91m',
						i=int(iteration_index),
						available_changes=available_changes_int,
						changes=changes,
						current_temp=curr_temp,
						curr_eval=original_curr_eval,
						candidate_eval=candidate_eval,
						diff=diff,
						criteria=criteria
					)
				)

				# results.register_iteration(i, changes, candidate_eval, curr_eval, metropolis_condition, curr_temp, diff)

		results.set_best(curr, curr_eval)
		return results
