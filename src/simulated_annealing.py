from math import floor
from numpy.random import rand
import pandas as pd

from input_parser import InputParser

from optimization_type import Maximization, Minimization
from results import Results

# https://machinelearningmastery.com/simulated-annealing-from-scratch-in-python/
class SimulatedAnnealing:
    def __init__(self, input_file, optimization) -> None:
        self.sequences_dictionary = InputParser.read_fasta_to_dict([input_file])

        self.initial: pd.DataFrame = InputParser.build_dataframe(
            self.sequences_dictionary
        )

        if optimization == "max":
            self.optimization = Maximization()
        elif optimization == "min":
            self.optimization = Minimization()
        else:
            raise BaseException("Invalid optimization name")

    def execute(
        self,
        n_iterations: int,
        initial_temp: int,
        score_function,
        generate_neighbor,
        is_debugging=False,
    ) -> Results:
        results: Results = Results()
        iterations_range: range = range(n_iterations)

        if is_debugging:
            print(
                "i\tcandidate_eval\tcurr_eval\tbest_eval\tmetropolis_condition\tcurrent_temp\tdiff"
            )

        # generate an initial point
        best = generate_neighbor(self.initial, changes=1)

        # evaluate the initial point
        best_eval = score_function.calculate_score(best)

        # current working solution
        curr = best
        curr_eval = best_eval

        for i in iterations_range:
            changes: int = floor(n_iterations / (i + 1))

            # take a step
            # candidate: pd.DataFrame = generate_neighbor(curr, changes=changes)
            candidate: pd.DataFrame = curr
            for j in range(5):
                candidate: pd.DataFrame = generate_neighbor(candidate, changes=changes)

            # evaluate candidate point
            candidate_eval: float = score_function.calculate_score(candidate)

            # difference between candidate and current point evaluation
            diff = candidate_eval - curr_eval
            optimization_condition = self.optimization.is_better_than_best(diff)

            # store new best point
            if optimization_condition:
                best = candidate
                best_eval = candidate_eval

            # calculate temperature for current epoch
            current_temp = initial_temp / float(i + 1)

            # check if we should keep the new point
            if optimization_condition:
                # store the new current point
                curr = candidate
                curr_eval = candidate_eval
                metropolis_condition = None
            else:
                metropolis_condition = self.optimization.metropolis(diff, current_temp)
                if rand() < metropolis_condition:
                    # store the new current point
                    curr = candidate
                    curr_eval = candidate_eval
            if is_debugging:
                print(
                    "{i}\t{candidate_eval}\t{curr_eval}\t{best_eval}\t{metropolis_condition}\t{current_temp}\t{diff}".format(
                        i=i,
                        candidate_eval=candidate_eval,
                        curr_eval=curr_eval,
                        best_eval=best_eval,
                        metropolis_condition=metropolis_condition,
                        current_temp=current_temp,
                        diff=diff,
                    )
                )

            results.register_iteration(
                i,
                candidate_eval,
                curr_eval,
                best_eval,
                metropolis_condition,
                current_temp,
                diff,
            )

        results.set_best(best, best_eval)
        return results
