from abc import abstractmethod
from math import exp

class OptimizationType:
	@abstractmethod
	def is_better_than_best(self, diff):
		pass


	@abstractmethod
	def metropolis(self, diff, current_temp):
		pass

class Maximization(OptimizationType):
	def is_better_than_best(self, diff):
		return diff >= 0


	def metropolis(self, diff, current_temp):
		"""
		exp( -(func(NEW) - func(OLD)) / T )
		exp( -(6 - 9)/T)
		exp( (-3)/T)
		exp(-3/T)
		"""
		try:
			return exp(diff / current_temp)
		except OverflowError:
			print("OVERFLOW WARN: diff=%f temp=%f" % (diff, current_temp))
			return 1

class Minimization(OptimizationType):
	"""
	diff = fx(NEW) - fx(OLD)
	diff < 0 => NEW < OLD
		eg: fx(NEW) - fx(OLD)
	    	9       - 12
			-3 < 0
			True

	diff > 0 => NEW > OLD
		eg: fx(NEW) - fx(OLD)
			12      - 9
			3 > 0
			False
	"""
	def is_better_than_best(self, diff):
		return diff <= 0


	def metropolis(self, diff, current_temp):
		"""
		exp( -(func(NEW) - func(OLD)) / T )
		exp( -(9 - 6) / T )
		exp( -(3) / T)
		exp(-3/T)
		"""

		try:
			if diff <= 0:
				raise OverflowError()

			return exp(-diff / current_temp)
		except OverflowError:
			print("OVERFLOW WARN: diff=%f temp=%f" % (diff, current_temp))
			return 1