package gapped_extension

import (
	"fmt"
	"github.com/Mal-Jovi/561_Project/utils"
)

func Left(q *string, d *[][]float64, hsp *[][]int, hit_thres float64, substitution_matrix *[][]int) {
	fmt.Println(*hsp)
}

func Right(q *string, d *[][]float64, hsp *[][]int, hit_thres float64, substitution_matrix *[][]int) {}

func extend() {}

func NeedlemanWunsch(seq1, seq2 *string,
					 S *[]string,
					 gap_open_penaty, gap_extend_penalty int,
					 substitution_matrix *[][]int) {

	N := init_nw(seq1, seq2, gap_open_penaty, gap_extend_penalty)
	S_idx := S_idx(S)

	for i := 1; i < len(* seq1); i++ {
		for j := 1; j < len(* seq2); j++ {
			(*N)[i][j] = utils.Max(
				(*N)[i-1][j-1] + sub_mat(i, j, seq1, seq2, S_idx, substitution_matrix),
				(*N)[i][j-1] + float64(gap_extend_penalty),
				(*N)[i-1][j] + float64(gap_extend_penalty),
			)
		}
	}
}

func backtrace(N *[][]float64) {

}

func sub_mat(i, j int, seq1, seq2 *string, S_idx *map[string]int, substitution_matrix *[][]int) float64 {
	return float64(
		(*substitution_matrix)[(*S_idx)[ string((*seq1)[i]) ]][(*S_idx)[ string((*seq2)[j]) ]],
	)
}

func SubstitutionMatrix(S *[]string) *[][]int {
	m := utils.MatInt(len(* S), len(* S))
	
	match_score := 1
	mismatch_score := -1

	for i := 0; i < len(* S); i++ {
		for j := 0; j < len(* S); j++ {
			if i == j {
				(* m)[i][j] = match_score
			
			} else {
				(* m)[i][j] = mismatch_score
			}
		}
	}
	return m
}

func S_idx(S *[]string) *map[string]int {
	S_idx := make(map[string]int)

	for i, word := range *S {
		S_idx[word] = i
	}
	return &S_idx
}

func init_nw(seq1, seq2 *string, gap_open_penaty, gap_extend_penalty int) *[][]float64 {
	N := utils.Mat(len(* seq1), len(* seq2))

	// fill first row and first column here
	for i := 0; i < len(*seq1); i++ {
		(*N)[i][0] = float64(i * gap_extend_penalty)
	}
	for j := 1; j < len(*seq2); j++ {
		(*N)[0][j] = float64(j * gap_extend_penalty)
	}
	return N
}