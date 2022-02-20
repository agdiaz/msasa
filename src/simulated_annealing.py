from math import floor
from numpy import exp
from numpy.random import rand
import pandas as pd

class Results:
	def __init__(self):
		self._records = pd.DataFrame(columns=['i', 'candidates', 'currents', 'bests', 'temperatures'], index=['i'])
		self._best = None
		self._best_eval = None

	def register_iteration(self, i, candidate_eval, curr_eval, best_eval, t):
		self._records = pd.concat([self._records, pd.DataFrame([[i ,candidate_eval, curr_eval, best_eval, t]], columns=self._records.columns)], ignore_index=True)

	def set_best(self, new_best, best_eval):
		self._best = new_best
		self._best_eval = best_eval

	def best(self):
		return [self._best, self._best_eval]

	def records(self, key):
		return self._records[key]

# https://machinelearningmastery.com/simulated-annealing-from-scratch-in-python/
class SimulatedAnnealing():
	def __init__(self, initial, n_iterations, initial_temp):
		self.initial = initial
		self.n_iterations = n_iterations
		self.initial_temp = initial_temp

	def maximize(self, score_function, generate_neighbor):
		results = Results()

		# generate an initial point
		best = self.initial

		# evaluate the initial point
		best_eval = score_function.calculate_score(best)

		# current working solution
		curr, curr_eval = best, best_eval

		# run the algorithm
		for i in range(self.n_iterations):
			changes = floor(self.n_iterations/(i + 1))

			# take a step
			candidate = generate_neighbor(curr, changes=changes)

			# evaluate candidate point
			candidate_eval = score_function.calculate_score(candidate)

			# check for new best solution
			if candidate_eval > best_eval:
				# store new best point
				best, best_eval = candidate, candidate_eval

			# difference between candidate and current point evaluation
			diff = candidate_eval - curr_eval

			# calculate temperature for current epoch
			current_temp = self.initial_temp / float(i + 1)

			# check if we should keep the new point
			try:
				# calculate metropolis acceptance criterion
				metropolis = exp(diff / current_temp)

				if diff > 0 or rand() < metropolis:
					# store the new current point
					curr, curr_eval = candidate, candidate_eval
			except OverflowError:
				pass
				# curr, curr_eval = candidate, candidate_eval

			results.register_iteration(i, candidate_eval, curr_eval, best_eval, current_temp)

		results.set_best(best, best_eval)
		return results


	def minimize(self, score_function, generate_neighbor):
		results = Results()

		# generate an initial point
		best = self.initial

		# evaluate the initial point
		best_eval = score_function.calculate_score(best)

		# current working solution
		curr, curr_eval = best, best_eval

		# run the algorithm
		for i in range(self.n_iterations):
			changes = floor(self.n_iterations/(i + 1))

			# take a step
			candidate = generate_neighbor(curr, changes=changes)
			# evaluate candidate point
			candidate_eval = score_function.calculate_score(candidate)

			# difference between candidate and current point evaluation
			diff = candidate_eval - curr_eval

			# check for new best solution
			if diff < 0:
				# store new best point
				best, best_eval = candidate, candidate_eval

			# calculate temperature for current epoch
			current_temp = self.initial_temp / float(i + 1)



			# check if we should keep the new point
			try:
				# calculate metropolis acceptance criterion
				metropolis = exp(-diff / current_temp)

				if diff < 0 or rand() < metropolis:
					# store the new current point
					curr, curr_eval = candidate, candidate_eval
			except OverflowError:
				pass
				# curr, curr_eval = candidate, candidate_eval

			results.register_iteration(i, candidate_eval, curr_eval, best_eval, current_temp)

		results.set_best(best, best_eval)

		return results