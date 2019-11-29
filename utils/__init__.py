import os     # Handle paths
import re     # Generate cartesian products (seeds)
import math   # Log probabilities
import random # Generate random sequence q
import itertools
import numpy as np

from tqdm import tqdm
from Bio import SeqIO
# from importlib.util import spec_from_file_location, module_from_spec


def gen_seeds(S, w):
    '''
    Generate iterable of seeds from alphabet S and word size w
    '''
    return itertools.product(S, repeat=w)


def find(substring, string):
    '''
    Return starting indices of all occurences
    of substring in string with overlap
    '''
    return tuple(
        match.start()
        for match in re.finditer(f'(?={substring})', string)
    )


def prob_find(seed, d, w, S, hit_thres):
    '''
    Parameters:
        d: database probabilistic sequence as numpy array
    '''
    if w > d.shape[1]:
        raise ValueError('Word size w cannot exceed length of database sequence d')

    thres = w * math.log(hit_thres)
    # indices = []

    indices = tuple(
        i for i in range(d.shape[1] - w + 1)
        if sum( log( d[ S.index( seed[j] ), i+j ] ) for j in range(w) ) >= thres
    )

    # for i in range(d.shape[1] - w + 1):
    #     # Calculate sum of log probabilities
    #     score = sum(
    #         math.log( d[ S.index( seed[j] ), i+j ] )
    #         for j in range(w)
    #     )

    #     if score >= thres:
    #         indices.append(i)
    
    return indices


def log(x):
    return math.log(x) if x > 0 else float('-inf')


def random_str(charset, length):
    '''
    Generate random string from charset of fixed length
    '''
    return ''.join(random.choice(charset) for _ in range(length))


def chunks(iterable, size=10):
    '''
    Split generator into generators of size <= 10
    '''
    iterator = iter(iterable)
    for first in iterator:
        yield itertools.chain([first], itertools.islice(iterator, size - 1))


def get_prob_seq(fasta, fasta_conf, S):
    '''
    Return as numpy array
    '''
    seq = seq_from_fasta(fasta)
    conf = get_conf(fasta_conf)

    assert len(seq) == len(conf)

    print('Generating probabilistic sequence...')

    prob_seq = np.zeros((len(S), len(seq)))

    for col in tqdm(range(len(seq))):
        for row in range(len(S)):
            if seq[col] == S[row]:
                prob_seq[row, col] = conf[col]
            else:
                prob_seq[row, col] = (1 - conf[col]) / (len(S) - 1)

    print('Done.')
    return prob_seq


def seq_from_fasta(fasta):
    '''
    Return sequence as string
    '''
    record = list(SeqIO.parse(open(fasta), 'fasta'))[0]
    return str(record.seq)


def get_conf(fasta_conf):
    '''
    Return confidence values as list of floats
    '''
    return tuple(float(conf) for conf in open(fasta_conf, 'r').read().strip().split())


def e_value(hsp_score):
    pass