import os            # CHeck if database index exists
import argparse
import itertools
import pickle as pkl # Export database index

from concurrent.futures import ProcessPoolExecutor
from functools import partial
from Bio import SeqIO
from time import time

import utils


def main():
    parser = argparse.ArgumentParser(description='BLAST')
    parser.add_argument('--config', '-c', type=str, default='config.py', help='[string] Path of config file')
    args = parser.parse_args()
    params = utils.import_from_file(args.config)

    q = utils.seq_from_fasta(params.q)
    d = utils.seq_from_fasta(params.d)
    # q = utils.random_str(charset=params.S, length=round(len(d)/10000)) # TBD
    # open('query.fa', 'w').write(q)
    
    blast(q, d, params.w, params.eps, params.S)


def blast(q, d, w, eps, S):
    '''
    Paramters:
        q:   query
        d:   database
        w:   word size
        eps: e-value to determine high scoring pairs (HSPs)
        S:   alphabet
    '''
    
    index = index_table(d, w, S)

    # extend(q, d, w, eps, index)
    # extend_gap()


def index_table(d, w, S):
    '''
    Index database if doesn't already exist at path
    Return index table as dict
    '''
    start = time()
    path = f'data/processed/index_table_w_{w}.pkl'

    if os.path.exists(path):
        return pkl.load(open(path, 'rb'))
     
    print(f'Indexing database with word size {w}...')

    seeds = itertools.product(S, repeat=w)

    with ProcessPoolExecutor() as executor:
        index = dict(filter(
            lambda x: x,
            executor.map(partial(_index_table, d), seeds)
        ))
    pkl.dump(index, open(path, 'wb'))
    
    print(f'Time elapsed: {time() - start:.2f}s')
    return index


def _index_table(d, seed):
    seed = ''.join(seed)
    indices = utils.find(seed, d)
    # return { wmer: indices } if indices else {}
    return (seed, indices) if indices else ()
    # return (wmer, indices)


def index_trie(d, w, S):
    pass


def extend(q, d, w, eps, index):
    '''
    Ungapped extension phase
    '''
    hsps = []
    
    for i in range(len(q) - w + 1):
        seed = q[i : i + w]
        
        for match in index[seed]:
            # Ungapped extension here
            pass



def left():
    pass


def right():
    pass


def extend_gap():
    '''
    Gapped extension phase
    '''
    pass


if __name__ == '__main__':
    main()