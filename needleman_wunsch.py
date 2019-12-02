def needleman_wunsch(seq1, seq2,
                     gap_open_penalty, gap_extend_penalty,
                     match_score, mismatch_score, substitution_matrix):
    
    N = 