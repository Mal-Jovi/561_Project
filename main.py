from Bio import SeqIO


def main():
    fasta = 'chr22.maf.ancestors.42000000.complete.boreo.fa'
    record = list(SeqIO.parse(open(f'data/{fasta}'), 'fasta'))[0]


def blast(q, d, w, c):
    '''
    :param q: query
    :param d: database
    :param w: word size
    :param c: cutoff score to determine high scoring pairs (HSPs)
    '''
    index()
    ungapped_extension()
    gapped_extension()


def index():
    pass


def ungapped_extension():
    pass


def gapped_extension():
    pass


if __name__ == '__main__':
    main()