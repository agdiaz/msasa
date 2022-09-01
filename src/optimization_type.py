from abc import abstractmethod
from math import exp

# from numpy import exp
# import numexpr as ne

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

def metro(diff, temperature, current_random):
	power = - diff / temperature
	return exp(power) > current_random
	# return ne.evaluate("exp(power)", local_dict={"power": power}) > current_random