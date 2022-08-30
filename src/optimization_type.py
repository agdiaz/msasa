from abc import abstractmethod
from numpy import exp
import numexpr as ne

class OptimizationType:
	@abstractmethod
	def is_better_than_best(self, diff):
		pass

class Maximization(OptimizationType):
	def is_better_than_best(self, diff):
		return diff >= 0

class Minimization(OptimizationType):
	def is_better_than_best(self, diff):
		return diff <= 0

def metro(current, candidate, temperature, current_random):
	diff = current - candidate
	power = diff / temperature

	return ne.evaluate("exp(power)") > current_random