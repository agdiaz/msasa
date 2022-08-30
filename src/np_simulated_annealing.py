# import numpy as np
from numpy import arange, floor, int32, float64
from numpy.random import rand
from numpy.lib.recfunctions import append_fields
from random import randint

from input_parser import InputParser
from optimization_type import metro
from results import Results

# https://machinelearningmastery.com/simulated-annealing-from-scratch-in-python/
class NpSimulatedAnnealing():
	CRITERIA_BETTER = "\033[92mBETTER"
	CRITERIA_WORSE = "\033[91mWORSE"
	CRITERIA_EQUAL = "\033[93mEQUAL"
	CRITERIA_METROPOLIS = "\033[95mMETROPOLIS"

	def __init__(self, input_file, optimization):
		self.sequences_dictionary = InputParser.read_fasta_to_dict([input_file])
		self.initial = InputParser.build_np_array(self.sequences_dictionary)


	def execute(self, n_iterations: int, initial_temp: float, score_function, generate_neighbor, is_debugging=False):
		results: Results = Results()

		if is_debugging:
			print("i\tavail_changes\tchanges\tcurrent_temp\tcurr_eval\tcandidate_eval\tdiff\tcriteria")

		sequences_count    = self.initial.shape[0]
		iterations_array   = arange(n_iterations, dtype=int32)
		changes_array      = floor(sequences_count / ((0.0025 * iterations_array) + 1))
		temperatures_array = initial_temp / ((0.5 * iterations_array) + 1) # -iterations_array * (initial_temp/n_iterations) + initial_temp : TEMP_LINEAR
		randoms            = rand(n_iterations)

		# iterations = np.stack((iterations_array, changes_array, temperatures_array, randoms), axis = 1)
		iterations = append_fields(
			iterations_array,
			['changes', 'temperatures', 'randoms'],
			[changes_array, temperatures_array, randoms],
			usemask=False,
			dtypes=[int32, float64, float64]
		)

		# generate an initial point and evaluate the initial point
		curr, curr_eval = self.initial, score_function.np_calculate_score(self.initial)
		results.set_initial(curr, curr_eval)

		for iteration_index, available_changes, curr_temp, iteration_random in iterations:
			# evaluate candidate point
			changes        = 1 if available_changes < 2 else randint(1, available_changes)
			candidate      = generate_neighbor(curr, changes)
			candidate_eval = score_function.np_calculate_score(candidate)

			criteria = NpSimulatedAnnealing.CRITERIA_WORSE

			# difference between candidate and current point evaluation
			original_curr_eval = curr_eval
			diff               = candidate_eval - curr_eval

			if diff < 0:# self.optimization.is_better_than_best(diff):
				criteria = NpSimulatedAnnealing.CRITERIA_EQUAL

				if diff != 0:
					# store the new current point
					curr, curr_eval, criteria = candidate, candidate_eval, NpSimulatedAnnealing.CRITERIA_BETTER

			elif metro(curr_eval, candidate_eval, curr_temp, iteration_random):
				# store the new current point even if it is not better than the current one
				curr, curr_eval, criteria  = candidate, candidate_eval, NpSimulatedAnnealing.CRITERIA_METROPOLIS

			if is_debugging:
				print(
					"{i:0>5d}\t{available_changes:0>4d}\t{changes:0>4d}\t{current_temp:0>8f}\t{curr_eval:0>8f}\t{candidate_eval:0>8f}\t{diff:0>+8f}\t{criteria}\033[0m".format(
						i=iteration_index,
						available_changes=available_changes,
						changes=changes,
						current_temp=curr_temp,
						curr_eval=original_curr_eval,
						candidate_eval=candidate_eval,
						diff=diff,
						criteria=criteria
					)
				)

		results.set_best(curr, curr_eval)
		return results
