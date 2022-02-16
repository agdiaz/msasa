from math import floor
from numpy import exp
from numpy.random import rand

class SimulatedAnnealing():
    def __init__(self, initial, n_iterations, temp):
        self.initial = initial
        self.n_iterations = n_iterations
        self.temp = temp

    def maximize(self, score_function, generate_neighbor):
        # for the records
        candidates = []
        currents = []
        bests = []
        temperatures = []

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
            candidate = generate_neighbor(best, changes=changes)

            # evaluate candidate point
            candidate_eval = score_function.calculate_score(candidate)

            # check for new best solution
            if candidate_eval > best_eval:
                # store new best point
                best, best_eval = candidate, candidate_eval

            # difference between candidate and current point evaluation
            diff = candidate_eval - curr_eval

            # calculate temperature for current epoch
            t = self.temp / float(i + 1)

            # calculate metropolis acceptance criterion
            metropolis = exp(diff / t)

            # check if we should keep the new point
            if diff > 0 or rand() < metropolis:
                # store the new current point
                curr, curr_eval = candidate, candidate_eval

            bests.append(best_eval)
            temperatures.append(t)
            currents.append(curr_eval)
            candidates.append(candidate_eval)

        return [best, best_eval, bests, currents, candidates, temperatures]