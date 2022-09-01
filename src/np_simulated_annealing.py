# import numpy as np
from sys import stderr
from numpy import arange, int32, float64, geomspace, linspace
from numpy.random import rand
from numpy.lib.recfunctions import append_fields
from random import randint
import numexpr as ne

from input_parser import InputParser
from objective_functions import SequencesComparerFactory
from optimization_type import metro
from results import Results

# https://machinelearningmastery.com/simulated-annealing-from-scratch-in-python/
class NpSimulatedAnnealing():
	CRITERIA_BETTER = "\033[92mBETTER"
	CRITERIA_WORSE = "\033[91mWORSE"
	CRITERIA_EQUAL = "\033[93mEQUAL"
	CRITERIA_METROPOLIS = "\033[95mMETRO"
	BEST_UPDATED = "\033[96mBEST"

	def __init__(self, input_file, sequences_comparer, inner_loop = False):
		self.sequences_dictionary = InputParser.read_fasta_to_dict([input_file])
		self.initial = InputParser.build_np_array(self.sequences_dictionary)
		self.score_function  = SequencesComparerFactory.from_name(sequences_comparer, self.initial.shape)
		self.inner_loop = inner_loop


	def execute(self, n_iterations: int, initial_temp: float, generate_neighbor, is_debugging=False):
		if is_debugging:
			stderr.write("iteration,changes_done,temperature,best_score,candidate_score,score_difference,change_criteria,iterations_without_changes\n")

		sequences_count    = self.initial.shape[0]
		iterations_array   = arange(n_iterations, dtype=int32)
		no_changes_limit   = int(0.20 * n_iterations)
		no_changes_count   = 0

		randoms = rand(n_iterations)
		temperatures = linspace(initial_temp, 0.001, n_iterations)
		available_changes = geomspace(sequences_count, 1, n_iterations)
		changes = ne.evaluate("(randoms * available_changes)", local_dict={'randoms': randoms, 'available_changes': available_changes})

		iterations = append_fields(
			iterations_array,
			['changes', 'temperatures', 'randoms'],
			[changes, temperatures, randoms],
			usemask=False,
			dtypes=[int32, float64, float64]
		)

		# generate an initial point and evaluate the initial point
		curr = generate_neighbor(self.initial)
		curr_eval = self.score_function.np_calculate_score(curr)

		results = Results()
		results.set_initial(curr, curr_eval)

		for iteration_index, changes, curr_temp, iteration_random in iterations:
			# evaluate candidate point

			iteration_candidate, iteration_candidate_score = curr, curr_eval
			if self.inner_loop:
				for _inner_loop_index in range(10):
					candidate      = generate_neighbor(curr, changes)
					candidate_eval = self.score_function.np_calculate_score(candidate)

					if candidate_eval <= iteration_candidate_score:
						iteration_candidate, iteration_candidate_score = candidate, candidate_eval

				candidate, candidate_eval = iteration_candidate, iteration_candidate_score
			else:
				candidate      = generate_neighbor(curr, changes)
				candidate_eval = self.score_function.np_calculate_score(candidate)

			criteria = NpSimulatedAnnealing.CRITERIA_WORSE

			original_curr_eval = curr_eval
			diff               = candidate_eval - curr_eval

			if diff == 0:# self.optimization.is_better_than_best(diff):
				criteria = NpSimulatedAnnealing.CRITERIA_EQUAL
				no_changes_count += 1
			elif diff < 0:
				# store the new current point
				curr, curr_eval, criteria = candidate, candidate_eval, NpSimulatedAnnealing.CRITERIA_BETTER
				no_changes_count = 0
			elif metro(diff, curr_temp, iteration_random):
				# store the new current point even if it is not better than the current one
				curr, curr_eval, criteria = candidate, candidate_eval, NpSimulatedAnnealing.CRITERIA_METROPOLIS
				no_changes_count = 0
			else:
				no_changes_count += 1

			if is_debugging:
				stderr.write(
					"{iteration_index:0>5d},{changes:0>4d},{current_temp:0>6f},{curr_eval:0>6f},{candidate_eval:0>6f},{diff:0>+6f},{criteria}\033[0m,{no_changes_count:0>5d}\n".format(
						iteration_index=iteration_index,
						no_changes_count=no_changes_count,
						changes=changes,
						current_temp=curr_temp,
						curr_eval=original_curr_eval,
						candidate_eval=candidate_eval,
						diff=diff,
						criteria=criteria
					)
				)

			if no_changes_count == no_changes_limit:
				break

		results.set_best(curr, curr_eval)
		return results
