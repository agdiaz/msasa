from argparse import Namespace
from input_parser import InputParser
# from objective_functions import SequencesComparer, SequencesComparerFactory

from timeit import default_timer as timer
from datetime import timedelta
from results import Results
# from simulated_annealing import SimulatedAnnealing

# import matplotlib
# matplotlib.use('Agg')
# from matplotlib import pyplot as plt
# plt.ioff()

class Runner:


	def __init__(self, args: Namespace, msa_algorithm):
		self.execution_id: int = args.execution_id
		self.input_file: str = args.input_file
		self.output_file: str = args.output_file
		self.output_best_plot: str = args.output_best_plot
		self.output_temp_plot: str = args.output_temp_plot
		self.sequences_comparer_name: str = args.sequences_comparer
		self.n_iterations: int = args.n_iterations
		self.inner_loop: str = args.inner_loop
		self.initial_temp: int = args.initial_temp
		self.is_debugging: str = args.is_debugging or False
		self.msa = msa_algorithm(self.input_file, args.sequences_comparer, self.inner_loop)


	def start(self):
		self.__start = timer()
		results = self.msa.execute(
			self.n_iterations,
			self.initial_temp,
			self.is_debugging
		)
		self.__end = timer()

		self._save_results_to_file(results)
		self._print_results_csv_row(results)

		# if self.output_best_plot:
		# 	self._plot_best_results(results)
		# 	self._plot_best_temp_results(results)
		# 	self._plot_diff_results(results)


		# if self.output_temp_plot:
		# 	self._plot_temp_results(results)


	def _print_results_csv_row(self, results: Results):
		initial, initial_score = results.initial()
		best, best_score = results.best()

		diff = best_score - initial_score
		percent_decrease = 100 * (diff / abs(initial_score))

		if self.is_debugging:
			print("output_file;initial_score;best_score;execution_id;time;difference;percent decrease")

		print('%s;%f;%f;%d;%s;%f;%f%%' % (
			self.output_file,
			initial_score,
			best_score,
			self.execution_id,
			timedelta(seconds=self.__end - self.__start),
			diff,
			percent_decrease
		))


	def _save_results_to_file(self, results: Results):
		best, __score = results.best()
		# if self.engine == "pandas":
		# 	InputParser.dataframe_to_msa_file(best, self.output_file)
		# else:
		InputParser.np_array_to_msa_file(best, self.msa.sequences_dictionary, self.output_file)

		if self.is_debugging:
			print("Best MSA:")

			with open(self.output_file, "rU") as input_handle:
				print(input_handle.read())

	# def _plot_best_results(self, results: Results):
	# 	fig, ax = plt.subplots(figsize=(25, 10))

	# 	plt.title("Simulated Annealing - MSA prediction ({0}) - Best results over the time".format(self.sequences_comparer_name))
	# 	ax.plot(results.records("bests"), color="green", label="Best")
	# 	# ax.plot(results.records("candidates"), color="orange", label="Candidate")

	# 	ax.set_xlabel('Iterations')
	# 	ax.set_ylabel('Score', color='green')

	# 	ax.legend()
	# 	fig.savefig(self.output_best_plot)

	# def _plot_diff_results(self, results: Results):
	# 	fig, ax = plt.subplots(figsize=(25, 10))

	# 	plt.title("Simulated Annealing - MSA prediction ({0}) - Diff over the time".format(self.sequences_comparer_name))

	# 	ax.plot(results.records("iterations"), results.records("diff"), color="red", label="Diff")
	# 	ax.axhline(y=0, color='black', linestyle='-')

	# 	axBests = ax.twinx()
	# 	axBests.plot(results.records("bests"), color='green', label="Best")

	# 	ax.set_xlabel('Iterations')
	# 	ax.set_ylabel('Diff (new - curr)', color='red')
	# 	axBests.set_ylabel('Best score', color='green')

	# 	ax.legend()
	# 	axBests.legend()

	# 	fig.savefig(self.output_best_plot.replace(".png", ".diff.png"))

	# def _plot_best_temp_results(self, results: Results):
	# 	fig, ax = plt.subplots(figsize=(25, 10))

	# 	plt.title("Simulated Annealing - MSA prediction ({0}) - Best results over the temperature".format(self.sequences_comparer_name))
	# 	ax.plot(results.records("temperatures"), results.records("bests"), color="green", label="Best")

	# 	ax.set_xlabel('Temperature')
	# 	ax.invert_xaxis()
	# 	ax.set_ylabel('Score', color='green')
	# 	ax.legend()
	# 	fig.savefig(self.output_best_plot.replace(".png", ".temp.png"))


	# def _plot_temp_results(self, results: Results):
	# 	fig, ax = plt.subplots(figsize=(25, 10))
	# 	plt.title("Simulated Annealing - MSA prediction ({0}) - Metropolis condition and temperature over the time".format(self.sequences_comparer_name))
	# 	ax.plot(results.records("temperatures"), results.records("metropolis"), color="green", label="Metropolis")

	# 	axBests = ax.twinx()
	# 	axBests.plot(results.records("temperatures"), results.records("bests"), color='red', label="Best")

	# 	ax.set_xlabel('Temperature')
	# 	ax.invert_xaxis()
	# 	ax.set_ylabel('Metropolis', color='green')
	# 	axBests.set_ylabel('Best score', color='red')

	# 	fig.savefig(self.output_temp_plot)