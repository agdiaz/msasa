from numpy.random import rand
from numpy import ndarray

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


	def execute(self, n_iterations: int, initial_temp: int, score_function, generate_neighbor, is_debugging=False):
		results: Results = Results()
		iterations_range: range = range(n_iterations)

		if is_debugging:
			print(
				"i\tcurrent_temp\tbest_eval\tcurr_eval\tcandidate_eval\tdiff"
			)

		# generate an initial point
		best = generate_neighbor(self.initial, changes=1)

		# evaluate the initial point
		best_eval: float = score_function.np_calculate_score(best)

		# current working solution
		curr = best.copy()
		curr_eval: float = best_eval

		for i in iterations_range:
			# take a step
			candidate = generate_neighbor(curr.copy())
			# breakpoint()

			# evaluate candidate point
			candidate_eval: float = score_function.np_calculate_score(candidate)

			# # difference between candidate and current point evaluation
			# diff = candidate_eval - curr_eval
			# optimization_condition = self.optimization.is_better_than_best(diff)

			# store new best point
			if candidate_eval < best_eval: # optimization_condition:
				best = candidate
				best_eval = candidate_eval

			# calculate temperature for current epoch
			metropolis_condition = None
			random_value = None

			# check if we should keep the new point
			current_temp = initial_temp / float(i + 1)

			# # difference between candidate and current point evaluation
			diff = candidate_eval - curr_eval
			optimization_condition = self.optimization.is_better_than_best(diff)

			if optimization_condition:
				# store the new current point
				curr = candidate
				curr_eval = candidate_eval
			else:
				metropolis_condition = self.optimization.metropolis(diff, current_temp)
				random_value = rand()
				if random_value < metropolis_condition:
					# store the new current point
					curr = candidate
					curr_eval = candidate_eval

			if is_debugging:
				print(
					"{color}{i:0>5d}\t{current_temp:0>8f}\t{best_eval:0>8f}\t{curr_eval:0>8f}\t{candidate_eval:0>8f}\t{diff:0>+8f}\033[0m".format(
						color='\033[92m' if diff < 0 else '\033[93m' if diff == 0 else '\033[91m',
						i=i,
						candidate_eval=candidate_eval,
						curr_eval=curr_eval,
						best_eval=best_eval,
						random_value=random_value,
						metropolis_condition=metropolis_condition,
						current_temp=current_temp,
						diff=diff,
					)
				)

			results.register_iteration(i, candidate_eval, curr_eval, best_eval, metropolis_condition, current_temp, diff)

		results.set_best(best, best_eval)
		return results
