import os # CHeck if database index exists
import itertools
import pickle as pkl # Export database index
from Bio import SeqIO
from tqdm import tqdm
from pprint import pprint

import utils


def main():
    fasta = 'chr22.maf.ancestors.42000000.complete.boreo.fa'
    record = list(SeqIO.parse(open(f'data/raw/{fasta}'), 'fasta'))[0]

    S = ['A', 'T', 'G', 'C']

    q = utils.random_str(charset=S, length=round(len(record)/10000)) # TBD
    d = str(record.seq)
    w = 4
    c = None # TBD
    
    blast(q, d, w, c, S)


def blast(q, d, w, c, S):
    '''
    :param q: query
    :param d: database
    :param w: word size
    :param c: cutoff score to determine high scoring pairs (HSPs)
    :param S: alphabet
    '''
    
    index = index_table(d, w, S)

    ungapped_extension(q, d, w, index)
    gapped_extension()


def index_table(d, w, S):
    '''
    Index database if doesn't already exist at path
    Return index table as dict
    '''
    path = f'data/processed/index_table_w_{w}.pkl'

    if not os.path.exists(path):
        index = {}
        
        print('Indexing database...')

        for wmer in tqdm(itertools.product(S, repeat=w), total=len(S)**w):
            wmer = ''.join(wmer)
            index[wmer] = utils.find(wmer, d)

        pkl.dump(index, open(path, 'wb'))
    
    return pkl.load(open(path, 'rb'))


def index_trie(d, w, S):
    pass


def ungapped_extension(q, d, w, index):
    for i in range(len(q)-w+1):
        wmer = q[i:i+w]
        
        for match in index[wmer]:
            # Ungapped extension here
            pass

def gapped_extension():
    pass


if __name__ == '__main__':
    main()