import os            # CHeck if database index exists
import itertools
import pickle as pkl # Export database index

from concurrent.futures import ProcessPoolExecutor
from functools import partial, reduce
from Bio import SeqIO
from tqdm import tqdm
from time import time
from pprint import pprint

import utils


def main():
    fasta = 'chr22.maf.ancestors.42000000.complete.boreo.fa'
    record = list(SeqIO.parse(open(f'data/raw/{fasta}'), 'fasta'))[0]

    S = ['A', 'T', 'G', 'C']

    q = utils.random_str(charset=S, length=round(len(record)/10000)) # TBD
    d = str(record.seq)
    w = 4
    eps = None # TBD
    
    blast(q, d, w, eps, S)


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

    wmers = itertools.product(S, repeat=w)

    with ProcessPoolExecutor() as executor:
        index = dict(filter(
            lambda x: x,
            executor.map(partial(_index_table, d), wmers)
        ))

    pkl.dump(index, open(path, 'wb'))
    print(f'Time elapsed: {time() - start:.2f}s')
    return index


def _index_table(d, wmer):
    wmer = ''.join(wmer)
    indices = utils.find(wmer, d)
    # return { wmer: indices } if indices else {}
    return (wmer, indices) if indices else ()


def index_trie(d, w, S):
    pass


def extend(q, d, w, eps, index):
    '''
    Ungapped extension phase
    '''
    hsps = []
    
    for i in range(len(q) - w + 1):
        wmer = q[i : i + w]
        
        for match in index[wmer]:
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