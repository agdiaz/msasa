import pandas as pd

class Results:


	def __init__(self):
		self._records = pd.DataFrame(columns=['i', 'candidates', 'currents', 'bests', 'metropolis', 'temperatures'], index=['i'])
		self._best = None
		self._best_eval = None


	def register_iteration(self, i, candidate_eval, curr_eval, best_eval, metropolis_condition, temp):
		self._records = pd.concat([self._records, pd.DataFrame([[i ,candidate_eval, curr_eval, best_eval, metropolis_condition, temp]], columns=self._records.columns)], ignore_index=True)


	def set_best(self, new_best, best_eval):
		self._best = new_best
		self._best_eval = best_eval


	def best(self):
		return [self._best, self._best_eval]


	def records(self, key):
		return self._records[key]
