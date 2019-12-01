import os
import time
import argparse
import pickle as pkl

from concurrent.futures import ProcessPoolExecutor
from functools import partial
# from Bio.Blast import NCBIWWW, NCBIXML
from pprint import pprint

import utils
from utils.params import params


def main():
    parser = argparse.ArgumentParser(description='Probabilistic BLAST')
    parser.add_argument('--config', '-c', type=str, default='config.py', help='[string] Path of config file')
    args = parser.parse_args()
    params.configure(args.config)

    d = utils.get_prob_seq(params.d, params.d_conf, params.S)
    q = utils.seq_from_fasta(params.q)

    prob_blast(
        q, d,
        params.w, params.S,
        params.hit_thres, params.delta,
        params.hsp_thres, params.e_thres,
    )

    # result_handle = NCBIWWW.qblast('blastn', 'nt', open('query.py', 'r').read())
    # result_handle = NCBIWWW.qblast('blastn', 'nt', q)
    # with open('res.xml', 'w') as fout:
    # fout.write(result_handle.read())

    # blast_recs = NCBIXML.parse(open('res.xml', 'r'))
    # for rec in blast_recs:
    #     pprint(rec.__dict__)


def prob_blast(q, d, w, S, hit_thres, delta, hsp_thres, e_thres):
    '''
    Return local alignments
    '''
    alignments = []
    index = prob_index_table(d, w, S, hit_thres)

    hsps = prob_extend(q, d, w, index, hit_thres, delta, hsp_thres, e_thres)
    prob_extend_gap(hsps, q, d, w, hit_thres, hsp_thres, delta)

    return alignments


def prob_index_table(d, w, S, hit_thres):
    '''
    Index database if doesn't already exist at path
    Return index table as dict
    '''
    start = time.time()
    path = f'data/processed/prob_index_table.w{w}.hit_thres{hit_thres}.pkl'

    if os.path.exists(path):
        return pkl.load(open(path, 'rb'))

    print(f'Indexing database with word size {w}...')

    seeds = utils.gen_seeds(S, w)

    with ProcessPoolExecutor() as executor:
        index = dict(filter(
            lambda x: x,
            executor.map(
                partial(_prob_index_table, d, w, S, hit_thres),
                seeds
            )
        ))
    pkl.dump(index, open(path, 'wb'))

    print(f'Time elapsed: {time.time() - start:.2f}s')
    return index


def _prob_index_table(d, w, S, hit_thres, seed):
    seed = ''.join(seed)
    indices = utils.prob_find(seed, d, w, S, hit_thres)
    return (seed, indices) if indices else ()


def prob_extend(q, d, w, S, index, hit_thres, delta, hsp_thres, e_thres):
    '''
    Ungapped extension
    '''
    hsps = []

    for q_idx in range(len(q) - w + 1):
        seed = q[q_idx:q_idx + w]

        try:  # In case seed not in index table
            for d_idx in index[seed]:
                hsp_left, score_left = prob_left(q_idx, d_idx, q, d, w, S, hit_thres, delta)
                hsp_right, score_right = prob_right(q_idx, d_idx, q, d, w, S, hit_thres, delta)

                if score_left == score_right == float('-inf'):
                    continue

                elif score_left == float('-inf'):
                    hsp_score = score_right

                elif score_right == float('-inf'):
                    hsp_score = score_left

                else:
                    hsp_score = score_right + score_left

                # Need to figure out e-value here
                # Is this right?
                # if utils.e_value(hsp_score) < e_thres:
                #     hsps.append((hsp_left, hsp_right))

                # Or this: hsp_score and ignore e_value
                if hsp_score > hsp_thres:
                    hsps.append((hsp_left, hsp_right))

        except KeyError:
            continue

    return hsps


def prob_left(q_idx, d_idx, q, d, w, S, hit_thres, delta):
    return _prob_extend(q_idx, d_idx, q, d, w, S, hit_thres, delta, step=-1)


def prob_right(q_idx, d_idx, q, d, w, S, hit_thres, delta):
    return _prob_extend(q_idx, d_idx, q, d, w, S, hit_thres, delta, step=1)


def _prob_extend(q_idx, d_idx, q, d, w, S, hit_thres, delta, step):
    max_score = float('-inf')
    max_q_idx = q_idx
    max_d_idx = d_idx
    score = 0

    while 0 < q_idx < len(q) - 1 and 0 < d_idx < d.shape[1] - 1:
        q_idx += step
        d_idx += step

        score += d[S.index(q[q_idx]), d_idx] - hit_thres
        if score > max_score:
            max_score = score
            max_q_idx = q_idx
            max_d_idx = d_idx

        elif max_score - score > 10:
            break

    return (max_q_idx, max_d_idx), max_score


def prob_extend_gap(highest_scoring_pairs, q, d, w, hit_thres, hsp_thres, delta):
    '''
    Gapped extension
    '''

    b = 0.7  # base gap cost
    e = 0.2  # gap extension cost

    hsa = []  # Highest Scoring Alignments

    for hsp in highest_scoring_pairs:
        l_index = hsp[0]
        r_index = hsp[1]

        q_index = q[l_index, r_index]
        # Perform Left gapped extension
        if r_index != 0:
            score_right = prob_right_g(l_index, r_index, q_index, d, w, hit_thres, hsp_thres, delta, b, e)

        if r_index != d.length - 1:
            score_left = prob_left_g(l_index, r_index, q_index, d, w, hit_thres, hsp_thres, delta, b, e)

        hsa_score = score_right + score_left

        # hsa.append((leftIndex, rightIndex, score))

    pass


def prob_right_g(l_index, r_index, q_idx, d_idx, q, d, w, S, hit_thres, delta, b, e):
    m = r_index - l_index + 1
    n = d.shape[1] - 1

    nw_table = [m][n]

    nw_table[0][0] = 0

    i = 1
    j = 1

    while i < m:
        nw_table[i][0] = nw_table[i - 1][0] - b
        i += 1
        

    while j < n:
        nw_table[0][j] = nw_table[0][j - 1] - b
        j+= 1

    for i in range(m):
        for j in range(n):
            '''
            1) perform column row comparison by comparing the character in the query row with the associated probability for the same character in the probablisitc table for the current column position
            2) the hsp right index + column number + 1 is the index for the database
            3) hsp right index + row number + 1 is the index for the query
            4) if the comparison meets the threshold, use the same normal NW scoring scheme?
            5) complex backtrace to generate score or just to note the path? And then perhaps forward trace through the path, generating the score and keeping track of whether or not the current score passes the delta threshold
            6) when delta threshold met or when end of row or column range, then stop
            7) return HSA alignment? How to return exact alignment like BLAST website does?
            '''
            return

    return 0


def prob_left_g(l_index, r_index, q_idx, d_idx, q, d, w, S, hit_thres, delta):
    return 0


if __name__ == '__main__':
    main()
