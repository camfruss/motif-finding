
#!/usr/bin/env python3

''' Performs Gibbs sampling on the specified set of sequences to identify motifs

Arguments:
    -f: set of sequences to identify a motif in
    -k: the length of the desired motif to find
Outputs:
    - a list of start positions for each sequence, along with the corresponding
      motif identified
    - the consensus motif (based on majority base at each position)
    - the position probability matrix

Example Usage:
    python3 1a.py -f motif1.fasta -k 10
'''

import argparse
from copy import deepcopy
import random
import numpy as np
from math import exp

np.random.seed(2950)

def sumLogProbs(a, b):
    if a > b:
        return a + np.log1p(exp(b-a))
    else:
        return b + np.log1p(exp(a-b))


''' Reads in the sequences from the motif files.

Arguments:
    filename: which filename to read sequences from
Returns:
    output: list of sequences in motif file
'''


def read_fasta(filename):
    with open(filename, "r") as f:
        output = []
        s = ""
        for l in f.readlines():
            if l.strip()[0] == ">":
                # skip the line that begins with ">"
                if s == "":
                    continue
                output.append(s)
                s = ""
            # keep appending to the current sequence when the line doesn't begin
            # with ">"
            else:
                s += l.strip()
        output.append(s)
        return output


''' Returns the majority vote (consensus) motif sequence.

Arguments:
    motif_sequences: list of all motif instances (list of strings)
Returns:
    majority sequence: consensus motif sequence (string)
'''


def majority(motif_sequences):
    base_ordering = {c: i for i, c in enumerate("ACGT")}
    freqs = np.zeros((len(motif_sequences[0]), 4))
    for m_s in motif_sequences:
        for i, b in enumerate(m_s):
            freqs[i, base_ordering[b]] += 1
    return ''.join(['ACGT'[np.argmax(row)] for row in freqs])


ACGT = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

'''Adds the characters in a motif to an overall distribution array.    
Arguments:
    pi: the distribution array to increment
    sequence: the string containing the motif
    start: the start point
    k: the motif length
    add_subtr: whether to add or subtract 1 from pi
Returns:
    none
'''


def add_motif(pi, sequence, start, k, add_subtr):
    for i in range(0, k):
        current_char = sequence[start + i]
        pi[i][ACGT[current_char]] += add_subtr


'''Calculate the probability of a motif given a model.
Arguments: 
    model: the motif probability matrix
    str: the string containing the motif
    start: the start point of the motif
    k: the length of the motif
Returns:
    pr: the probability of the motif given the model
'''


def motif_prob(model, str, start, k, denom):
    pr = 0
    # probability of a motif = product of probability of each character
    for offset in range(0, k):
        current_char = str[start + offset]
        pr += np.log(model[offset][ACGT[current_char]])
    pr -= (k * denom)
    return pr


'''Calculate the likelihood of a set of motifs given their model.
Arguments: 
    model: the motif probability matrix
    starts: the start points of the motifs
    k: the length of the motif
    sequences: the list of sequences for the motifs
Returns:
'''


def likelihood(model, starts, k, sequences, denom):
    pr = 0
    # no need for denom if it's going to be constant between comparisons anyways
    for i in range(0, len(sequences)):
        pr += motif_prob(model, sequences[i],
                         starts[i], k, denom)
    return pr


''' Performs Gibbs sampling to estimate the start positions of the motif in
    each sequence and the motif model.

Arguments:
    sequences: list of sequences to find motifs in
    k: length of motif
    epsilon: pseudocounts for each base
Returns:
    starts: list of start positions for each sequence 
    pi: corresponding motif model (incorporating pseudocounts).
'''


def gibbs_sampling(sequences, k, epsilon):
    # initialize random motif starts (1 for each sequence, uniform from 0 to seqlen - k)
    num_seqs = len(sequences)
    small_denom = np.log(num_seqs - 1 + 4 * epsilon)
    big_denom = np.log(num_seqs + 4 * epsilon)
    starts = np.zeros(num_seqs, dtype=int)
    # for i in range(0, len(starts)):
    #     rand = np.random.randint(0, len(sequences[i]) - k + 1)
    #     starts[i] = rand
    # initialize ACGT distributions for each position in the motif
    pi = np.full((k, 4), epsilon)

    for j in range(0, num_seqs):
        add_motif(pi, sequences[j], starts[j], k, 1)

    best_starts = np.copy(starts)
    best_model = np.copy(pi)
    best_likelihood = likelihood(
        best_model, best_starts, k, sequences, big_denom)
    time_since_change = 0

    motif_distr = []
    for i in range(0, num_seqs):
        motif_distr.append(np.zeros(len(sequences[i]) - k + 1))

    # while condition holds:
    # for each motif n, update based on relative frequency in other sequences:
    # subtract current motif characters from distribution add_motif(n, -1), draw new start from remainder
    # update starts with new start:
    # add new motifs back to distr (add_motif(n, 1))
    first = True
    while(time_since_change < 100):
        for i in range(0, num_seqs):
            curr_str = sequences[i]
            add_motif(pi, sequences[i], starts[i], k, -1)
            if(first):
                print(pi)
            current_distr = motif_distr[i]
            # draw new start from remainder
            for curr_start in range(0, len(current_distr)):
                current_distr[curr_start] = motif_prob(
                    pi, curr_str, curr_start, k, small_denom)
            if(first):
                print(current_distr)
            # normalize before drawing: x - logsum(x1, x2 ... xn)
            logsum = current_distr[0]
            for curr_start in range(1, len(current_distr)):
                logsum = sumLogProbs(logsum, current_distr[curr_start])
            for j in range(0, len(current_distr)):
                current_distr[j] = exp(current_distr[j] - logsum)
            if(first):
                print(current_distr)
                first = False
            new_start = np.random.choice(len(current_distr), p=current_distr)
            starts[i] = new_start
            add_motif(pi, sequences[i], starts[i], k, 1)
        new_likelihood = likelihood(
            pi, starts, k, sequences, big_denom)
        if (new_likelihood > best_likelihood):
            best_likelihood = new_likelihood
            best_starts = np.copy(starts)
            best_model = deepcopy(pi)
            time_since_change = 0
        time_since_change += 1
    return best_starts, best_model / (num_seqs + 4 * epsilon)


def main():
    parser = argparse.ArgumentParser(
        description='Estimate the start positions and motif model via Gibbs sampling.')
    parser.add_argument('-f', action="store", dest="f",
                        type=str, default='sec9.fasta')
    parser.add_argument('-k', action="store", dest="k", type=int, default=4)
    parser.add_argument('-epsilon', action="store",
                        dest="epsilon", type=float, default=0.5)
    args = parser.parse_args()
    sequences = read_fasta(args.f)
    k = args.k
    epsilon = args.epsilon

    starts, pi = gibbs_sampling(sequences, k, epsilon)
    motif_sequences = [s[starts[i]:starts[i]+k]
                       for i, s in enumerate(sequences)]
    print('\n'.join([("Sequence %d: " % i) +
          m for i, m in enumerate(motif_sequences)]))
    print("Consensus motif: %s" % majority(motif_sequences))
    print(pi)


if __name__ == '__main__':
    main()