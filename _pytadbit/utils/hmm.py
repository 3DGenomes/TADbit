from __future__ import print_function
from numpy import log, pi as pi_num, exp
import sys

def best_path(probs, pi, T):
    """
    Viterbi algorithm with backpointers
    """
    n = len(T)
    m = len(probs[0])
    log_V     = [[0. for _ in range(m)] for _ in range(n)]
    backpt    = [[0  for _ in range(m)] for _ in range(n)]
    states    =  [0  for _ in range(m)]
    log_pi    =  [float('-inf') if pi[i] < 0. else log(pi[i]) for i in range(len(pi))]
    log_T     =  [float('-inf') if T [i] < 0. else log(T[i]) for i in range(len(T))]
    log_probs = [[float('-inf') if probs[i][j] < 0. else log(probs[i][j])
                  for j in range(m)]
                 for i in range(n)]
    for i in range(n):
        log_V[i][0] = log_probs[i][0] + log_pi[i]
    for k in range(1, m):
        for stat in range(n):
            log_V[stat][k] = float('-inf')
            for prev in range(n):
                # original state prob times transition prob
                prob = log_V[prev][k-1] + log_T[prev][stat]
                if prob > log_V[stat][k]:
                    log_V[stat][k] = prob
                    backpt[stat][k-1] = prev
            log_V[stat][k] += log_probs[stat][k]
    # get the likelihood of the most probable path
    prob = log_V[0][-1]
    for i in range(1, n):
        if log_V[i][-1] > prob:
            prob = log_V[i][-1]
            states[-1] = i
    # Follow the backtrack: get the path which maximize the path prob.
    for i in range(m - 2, -1, -1):
        states[i] = backpt[states[i + 1]][i]
    return states, prob

def baum_welch_optimization(xh, T, E, new_pi, new_T, corrector,
                            new_E, etas, gammas):
    """
    implementation of the baum-welch algorithm
    """
    n = len(T)
    m = len(gammas[0])
    for i in range(n):
        for j in range(n):
            new_pi[i] += etas[i][j][0]
    for i in range(n):
        for j in range(n):
            new_T[i][j] += sum(etas[i][j][k] for k in range(m - 1))

    for k in range(m):
        for i in range(n):
            gik = gammas[i][k]
            corrector[i] += gik
            new_E[i][0] += gik * xh[k]
            new_E[i][1] += gik * (xh[k] - E[i][0])**2

def update_parameters(corrector, pi, new_pi, T, new_T, E, new_E):
    """
    final round of the baum-welch
    """
    ### update initial probabilities
    n = len(T)
    total = 0.
    delta = 0.
    for i in range(n):
        total += new_pi[i]
    for i in range(n):
        new_pi[i] /= total
        delta = max(delta, abs(new_pi[i] - pi[i]))
        pi[i] = new_pi[i]
    ### update transitions
    for i in range(n):
        total = 0.
        for j in range(n):
            total += new_T[i][j]
        for j in range(n):
            new_T[i][j] /= total
            delta = max(delta, abs(new_T[i][j] - T[i][j]))
            T[i][j] = new_T[i][j]
    ### update emissions
    for i in range(n):
        # update the means
        if corrector[i] > 0.:
            new_E[i][0] /= corrector[i]
            new_E[i][1] /= corrector[i]
            # corrector[i] = 1.
            delta = max(delta, abs(new_E[i][0] - E[i][0]))
            E[i][0] = new_E[i][0]
            # update the stdevs
            delta = max(delta, abs(new_E[i][1] - E[i][1]))
            E[i][1] = new_E[i][1]
    return delta

def train(pi, T, E, observations, verbose=False, threshold=1e-6, n_iter=1000):
    delta = float('inf')
    for it in range(n_iter):
        # reset for new iteration
        new_pi = [0. for _ in pi]
        new_T  = [[0. for _ in i] for i in T]
        new_E  = [[0. for _ in i] for i in E]
        corrector = [0. for _ in range(len(T))]
        for h in range(len(observations)):
            probs  = gaussian_prob(observations[h], E)
            alphas, scalars = get_alpha(probs, pi, T)
            betas  = get_beta(probs, T, scalars)
            etas   = get_eta(probs, T, alphas, betas)
            gammas = get_gamma(T, alphas, betas)
            baum_welch_optimization(observations[h], T, E, new_pi, new_T, corrector,
                                    new_E, etas, gammas)
        delta = update_parameters(corrector, pi, new_pi, T, new_T, E, new_E)
        if verbose:
            print("\rTraining: %03i/%04i (diff: %.8f)" % (it, n_iter, delta), end=' ')
            sys.stdout.flush()
        if delta <= threshold:
            break
    if verbose:
        print("\n")

def get_eta(probs, T, alphas, betas):
    """
    for Baum-Welch: probability of being in states i and j at times t and t+1
    """
    n = len(T)
    m = len(probs[0])
    etas = [[[0.  for _ in range(m - 1)] for _ in range(n)] for _ in range(n)]
    for k in range(m - 1):
        tot = 0.
        for i in range(n):
            for j in range(n):
                etas[i][j][k] += alphas[i][k] * T[i][j] * probs[j][k+1] * betas[j][k+1]
                tot += etas[i][j][k]
        for i in range(n):
            for j in range(n):
                etas[i][j][k] /= tot
    return etas

def get_gamma(T, alphas, betas):
    """
    for Baum-Welch: probability of being in state i at time t
    """
    n = len(T)
    m = len(alphas[0])
    return [[alphas[i][k] * betas[i][k]  for k in range(m)] for i in range(n)]

def gaussian_prob(x, E):
    """
    of x to follow the gaussian with given E
    https://en.wikipedia.org/wiki/Normal_distribution
    """
    probs = []
    for d, (mu, sd) in enumerate(E):
        pi2sd = (2. * pi_num * sd)**-0.5
        inv2sd = 1. / (2. * sd)
        probs.append([])
        for obs in x:
            p = pi2sd * exp(-(obs - mu)**2 * inv2sd)
            probs[d].append(p)
    return probs

def get_alpha(probs, pi, T):
    """
    computes alphas using forward algorithm
    """
    n = len(T)
    m = len(probs[0])
    alphas = [[0. for _ in range(m)] for _ in range(n)]
    scalars = [0. for _ in range(m)]

    # initialize alpha for each state
    for i in range(n):
        alphas[i][0] = pi[i] * probs[i][0]
        scalars[0] += alphas[i][0]
    for i in range(n):
        alphas[i][0] /= scalars[0]
    for k in range(1, m):
        for i in range(n):
            for j in range(n):
                # all transition probabilities to become "i" times previous alpha
                alphas[i][k] += alphas[j][k-1] * T[j][i]
            # times probablity to belong to this states
            alphas[i][k] *= probs[i][k]
            scalars[k] += alphas[i][k]
        for i in range(n):
            alphas[i][k] /= scalars[k]
    return alphas, scalars

def get_beta(probs, T, scalars):
    """
    computes betas using backward algorithm
    """
    n = len(T)
    m = len(probs[0])
    # intialize beta at 1.0
    betas = [[0. for _ in range(m - 1)] + [1.0] for _ in range(n)]
    for k in range(-2, -m - 1, -1):
        for i in range(n):
            for j in range(n):
                betas[i][k] += betas[j][k + 1] * probs[j][k + 1] * T[i][j]
            betas[i][k] /= scalars[k + 1]
    return betas

