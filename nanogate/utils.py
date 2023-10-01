import screed
import mmh3
import itertools
import scipy
import scipy.stats
import random
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def build_kmers(sequence, ksize):
    """
    Method to obtain a list of k-mers given a sequence and a k-mer length
    """
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)

    return kmers


def read_kmers_from_file(filename, ksize):
    """
    Method to obtain all k-mers from a fasta or fastaq file
    """
    all_kmers = []
    for record in screed.open(filename):
        sequence = record.sequence
        kmers = build_kmers(sequence, ksize)
        all_kmers += kmers

    return all_kmers


def hash_kmer(kmer):
    """
    method to turn a k-mer into a hash code
    """
    # calculate the reverse complement
    rc_kmer = screed.rc(kmer)

    # determine whether original k-mer or reverse complement is lesser
    if kmer < rc_kmer:
        canonical_kmer = kmer
    else:
        canonical_kmer = rc_kmer

    # calculate murmurhash using a hash seed of 42
    hash = mmh3.hash64(canonical_kmer, 42)[0]
    if hash < 0: hash += 2 ** 64

    # done
    return hash


# def subsample_modulo(kmers, scaled):
#     """
#     Method to filter kmers given a scale value
#     """
#     max_hash =  2**64
#     keep_below = max_hash / scaled
#     keep = []
#     for kmer in kmers:
#         if hash_kmer(kmer) < keep_below:
#             keep.append(kmer)
#         # otherwise, discard
#     return keep


def jaccard_containment(a, b):
    """
    Method to get jaccard contatinet values
    given two set a and b
    """
    a = set(a)
    b = set(b)

    intersection = len(a.intersection(b))
    if len(a) > 0:
        jc = intersection / len(a)
    else:
        jc = 0
    return jc


def poisson_distribution(parts_list, construct_length):
    """
    Method to build a poisson distribution of
    reads length by reading parts length
    """
    part_lenghts = []
    sample_constructs = []
    for part in parts_list:
        part_lenghts.append(part.length)

    first_sample = generate_random_construct(part_lenghts, construct_length)
    sample_constructs.append(first_sample)
    epsilon = np.mean(sample_constructs)
    new_epsilon = 0
    epsilon_diff = epsilon - new_epsilon
    trials = 0
    while epsilon_diff > 0.001 and trials <20000:
        epsilon = np.mean(sample_constructs)
        sample_construct = generate_random_construct(part_lenghts, construct_length)
        sample_constructs.append(sample_construct)
        new_epsilon = np.mean(sample_constructs)
        epsilon_diff = abs(epsilon - new_epsilon)
        trials += 1

    return sample_constructs


def generate_random_construct(part_lengths, construct_length):
    """
    method to generate a possible construct given a parts library
    :return:
    """
    sample_part_lengths = random.sample(part_lengths, construct_length)
    sample_construct = sum(sample_part_lengths)
    return sample_construct




def evalute_likelihood(read, construct_poisson, reads_list):
    """
    Method to evaluate likelihood between reads length and constructs
    lenght reads too short to belong to the distribution are discarded
    """
    filtered_out = 0
    lambda_poisson = np.mean(construct_poisson)
    likelihood = scipy.stats.poisson.cdf(read.length, lambda_poisson)
    if likelihood > 0.05:
        reads_list.append(read)
    else:
         filtered_out = 1
    return filtered_out




def read_theshold(threshold):
    """
    convert string "02" into float 0.2
    :param threshold:
    :return:
    """
    threshold_float = 0
    threshold_float += float(threshold[0])
    threshold_float += (float(threshold[1])/10)
    return threshold_float


def build_parts_length_dictionary(parts):
    """
    build a dictionary {part:length}
    :return:
    """
    parts_length = {}
    for part in parts:
        parts_length[part.name] = part.length
    return parts_length
