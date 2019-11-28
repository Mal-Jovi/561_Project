# Fasta file for database sequence d
d      = 'data/raw/chr22.maf.ancestors.42000000.complete.boreo.fa'
d_conf = 'data/raw/chr22.maf.ancestors.42000000.complete.boreo.conf'

# Fasta file for query
q = 'query.fa'

# Size of seed
w = 4

# e-value to determine high scoring pairs (HSPs)
eps = None # TBD

# Alphabet
S = ['A', 'T', 'G', 'C']

# Scoring scheme
score = [
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
]

# Threshold for seed hit: w * log(hit_thres)
hit_thres = 0.7

# Delta threshold for ungapped extenion: stop ungapped extension if cur score/max score <= delta
delta = 0.25