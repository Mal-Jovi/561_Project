import os     # Handle paths
import re     # Generate cartesian products (seeds)
import random # Generate random sequence q
import itertools
import numpy as np

from tqdm import tqdm
from Bio import SeqIO
from importlib.util import spec_from_file_location, module_from_spec


def find(substring, string):
    '''
    Return starting indices of all occurences
    of substring in string with overlap
    '''
    return tuple(
        match.start()
        for match in re.finditer(f'(?={substring})', string)
    )


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
            if S[row] == seq[col]:
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
    return tuple(float(conf) for conf in open(fasta_conf, 'r').read().split())


def import_from_file(path):
    '''
    Returns a module object from a filepath
    '''
    if not os.path.exists(path):
        raise FileNotFoundError(f'{path} doesn\'t exist')
    
    name = path.split('/')[-1].split('.')[0]
    spec = spec_from_file_location(name, path)
    if spec is None:
        raise ValueError('could not load module from "%s"' % path)
    
    m = module_from_spec(spec)
    spec.loader.exec_module(m)
    return m