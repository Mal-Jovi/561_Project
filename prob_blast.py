import argparse
import utils

from Bio.Blast import NCBIWWW, NCBIXML
from pprint import pprint

def main():
    parser = argparse.ArgumentParser(description='Probabilistic BLAST')
    parser.add_argument('--config', '-c', type=str, default='config.py', help='[string] Path of config file')
    args = parser.parse_args()
    params = utils.import_from_file(args.config)

    d = utils.get_prob_seq(params.d, params.d_conf, params.S)
    q = utils.seq_from_fasta(params.q)

    prob_blast(q, d, params.w, params.eps, params.S)

    print(d[:,:5])

    # result_handle = NCBIWWW.qblast('blastn', 'nt', open('query.py', 'r').read())
    # result_handle = NCBIWWW.qblast('blastn', 'nt', q)
    # with open('res.xml', 'w') as fout:
        # fout.write(result_handle.read())

    # blast_recs = NCBIXML.parse(open('res.xml', 'r'))
    # for rec in blast_recs:
    #     pprint(rec.__dict__)


def prob_blast(q, d, w, eps, S):
    pass


def prob_index():
    pass


def prob_extend():
    pass


def prob_extend_gap():
    pass


if __name__ == '__main__':
    main()