import argparse
import numpy as np

from utils import import_from_file, seq_from_fasta, get_conf

def main():
    parser = argparse.ArgumentParser(description='Probabilistic Blast')
    parser.add_argument('--config', '-c', type=str, default='config.py', help='[string] Path of config file')
    args = parser.parse_args()

    params = import_from_file(args.config)
    d = get_prob_seq(params.d, params.d_conf, params.S)

    prob_blast()

def prob_blast():
    pass

def get_prob_seq(fasta, fasta_conf, S):
    '''
    Return as numpy array
    '''
    seq = seq_from_fasta(fasta)
    conf = get_conf(fasta_conf)

    assert len(seq) == len(conf)

    prob_seq = np.zeros((len(S), len(seq)))

    for col in range(len(seq)):
        for row in range(len(S)):
            if S[row] == seq[col]:
                prob_seq[row, col] = conf[col]
            else:
                prob_seq[row, col] = (1 - conf[col]) / (len(S) - 1)
    
    print(prob_seq[:,:5])

    return prob_seq

if __name__ == '__main__':
    main()