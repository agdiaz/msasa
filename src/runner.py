from argparse import Namespace
from generators import msa_neighbor_add_remove
from input_parser import InputParser
from objective_functions import SequencesComparer, SequencesComparerFactory
import matplotlib

from results import Results
from simulated_annealing import SimulatedAnnealing
matplotlib.use('Agg')
from matplotlib import pyplot as plt
plt.ioff()

class Runner:


	def __init__(self, args: Namespace, msa_type: SimulatedAnnealing):
		self.input_file: str = args.input_file
		self.output_file: str = args.output_file
		self.output_best_plot: str = args.output_best_plot
		self.output_temp_plot: str = args.output_temp_plot
		self.sequences_comparer_name: str = args.sequences_comparer
		self.n_iterations: int = args.n_iterations
		self.optimization: str = args.optimization
		self.initial_temp: int = args.initial_temp
		self.score_function: SequencesComparer = SequencesComparerFactory.from_name(args.sequences_comparer)
		self.msa: SimulatedAnnealing = msa_type(self.input_file, self.optimization)


	def start(self) -> None:
		results = self._execute_msa()

		self._save_results_to_file(results)
		self._print_results_csv_row(results)

		if self.output_best_plot:
			self._plot_best_results(results)
			self._plot_best_temp_results(results)


		if self.output_temp_plot:
			self._plot_temp_results(results)


	def _execute_msa(self) -> Results:
		results: Results = self.msa.execute(self.n_iterations, self.initial_temp, self.score_function, msa_neighbor_add_remove)

		return results


	def _print_results_csv_row(self, results: Results) -> None:
		initial_score = results.records("bests")[1]
		best, best_score = results.best()
		print('%s;%f;%f' % (self.output_file, initial_score, best_score))


	def _save_results_to_file(self, results: Results) -> None:
		best, __score = results.best()
		InputParser.dataframe_to_msa_file(best, self.output_file)


	def _plot_best_results(self, results: Results) -> None:
		fig, ax = plt.subplots(figsize=(25, 10))

		plt.title("Simulated Annealing - MSA prediction ({0}) - Best results over the time".format(self.sequences_comparer_name))
		ax.plot(results.records("bests"), color="green")
		ax.plot(results.records("candidates"), color="orange")

		ax.set_xlabel('Iterations')
		ax.set_ylabel('Score', color='green')

		fig.savefig(self.output_best_plot)

	def _plot_best_temp_results(self, results: Results) -> None:
		fig, ax = plt.subplots(figsize=(25, 10))

		plt.title("Simulated Annealing - MSA prediction ({0}) - Best results over the temperature".format(self.sequences_comparer_name))
		ax.plot(results.records("temperatures"), results.records("bests"), color="green")

		ax.set_xlabel('Temperature')
		ax.invert_xaxis()
		ax.set_ylabel('Score', color='green')

		fig.savefig(self.output_best_plot.replace(".png", ".temp.png"))


	def _plot_temp_results(self, results: Results) -> None:
		fig, ax = plt.subplots(figsize=(25, 10))
		plt.title("Simulated Annealing - MSA prediction ({0}) - Metropolis condition and temperature over the time".format(self.sequences_comparer_name))
		ax.plot(results.records("temperatures"), results.records("metropolis"), color="green")

		axBests = ax.twinx()
		axBests.plot(results.records("temperatures"), results.records("bests"), color='red')

		ax.set_xlabel('Temperature')
		ax.invert_xaxis()
		ax.set_ylabel('Metropolis', color='green')
		axBests.set_ylabel('Best score', color='red')

		fig.savefig(self.output_temp_plot)