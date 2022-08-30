# import pandas as pd

class Results:


	def __init__(self):
		# self._records = pd.DataFrame(columns=['iterations', 'changes', 'candidates', 'currents', 'metropolis', 'temperatures', 'diff'], index=['i'])
		self._best = None
		self._best_eval = None


	# def register_iteration(self, i, changes, candidate_eval, curr_eval, metropolis_condition, temp, diff):
    # 	self._records = pd.concat([self._records, pd.DataFrame([[i , changes, candidate_eval, curr_eval, metropolis_condition, temp, diff]], columns=self._records.columns)], ignore_index=True)


	def set_initial(self, initial, initial_eval):
		self._initial = initial
		self._initial_eval =  initial_eval


	def set_best(self, new_best, best_eval):
		self._best = new_best
		self._best_eval = best_eval


	def initial(self):
		return self._initial, self._initial_eval


	def best(self):
		return self._best, self._best_eval


	# def records(self, key):
	# 	return self._records[key]
