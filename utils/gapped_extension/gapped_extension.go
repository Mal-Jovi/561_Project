package gapped_extension

import (
	"fmt"
	// "math"
	"github.com/kr/pretty"
	"github.com/Mal-Jovi/561_Project/utils"
)

// Gapped extension by incrementally expanding Needleman-Wunsch
// dynamic programming table
//
// Parameters:
// -----------
//
// 	q_idx: Index of query where Needleman-Wunsch begins
// 	d_idx: Index pf databas sequence where Needleman-Wunsch begins
// 	S    : Alphabet e.g. ["A", "T", "G", "C"]
// 	S_idx: Mapping of each character to their index in alphabet `S`
// 	hit_thres: Probability threshold
// 	substitution_matrix: Return match and mismatch scores for each character pairing
// 	direction: Direction of extension; 1 for left-right, -1 for right-left
func Extend(q_idx, d_idx int,
			q *string,
			d *[][]float64,
			S *[]string,
			S_idx *map[string]int,
			hit_thres float64,
			// substitution_matrix *[][]int,
			direction int) {

	if !preconditions(q_idx, d_idx, q, d, direction) {
		fmt.Println("SAD")
		return
	}
	
	gap_penalty := 1.
	
	N := utils.Mat(1, 1)
	backptrs := utils.MatInt(1, 1)
	
	q_len, d_len := get_q_len_d_len(q_idx, d_idx, q, d, direction)

	if direction == 1 {
		fmt.Println("q seg:", (*q)[q_idx:q_idx + q_len])
		fmt.Println("d seg:", *utils.PrettyProbSeg(d, d_idx, d_idx + d_len, S, hit_thres))
	} else {
		fmt.Println("extending left")
		fmt.Println("q seg:", (*q)[:q_idx + 1])
		fmt.Println("d seq:", *utils.PrettyProbSeg(d, 0, d_idx + 1, S, hit_thres))
	}

	fill(N, backptrs, q_idx, d_idx, q_len, d_len, q, d, S_idx, gap_penalty, hit_thres, direction)
	
	fmt.Println("NW table:")
	pretty.Println(*N)
	fmt.Println("backpointers:")
	pretty.Println(*backptrs)
}

func fill(N *[][]float64,
		  backptrs *[][]int,
		  q_idx, d_idx, q_len, d_len int,
		  q *string,
		  d *[][]float64,
		  S_idx *map[string]int,
		  gap_penalty, 
		  hit_thres float64,
		  direction int) {

	// Incrementally expand dynamic programming and backpointers tables by layer
	for layer := 0; layer < max(q_len, d_len); layer++ {
		// Expand, then initialize top row and/or left column
		if layer < q_len {
			fill_bottom(N, backptrs, layer, q_idx, d_idx, q, d, S_idx, gap_penalty, hit_thres, direction)
		}
		if layer < d_len {
			fill_right(N, backptrs, layer, q_idx, d_idx, q, d, S_idx, gap_penalty, hit_thres, direction)
		}
	}
}

func fill_bottom(N *[][]float64,
				 backptrs *[][]int,
				 layer, q_idx, d_idx int,
				 q *string,
				 d *[][]float64,
				 S_idx *map[string]int,
				 gap_penalty,
				 hit_thres float64,
				 direction int) {

	// Expand bottom
	utils.ExpandMatBottom(N, 1)
	utils.ExpandMatIntBottom(backptrs, 1)

	// Initialize cell on leftmost column
	(*N)[layer + 1][0] = -float64(layer + 1) * gap_penalty

	// Fill new bottom row
	for col := 1; col < len((*N)[0]); col++ {
		nw_recurrence(N, backptrs, layer+1, col, q_idx + layer*direction, d_idx + (col-1)*direction, q, d, S_idx, gap_penalty, hit_thres)
	}
}

func fill_right(N *[][]float64,
				backptrs *[][]int,
				layer, q_idx, d_idx int,
				q *string,
				d *[][]float64,
				S_idx *map[string]int,
				gap_penalty,
				hit_thres float64,
				direction int) {
	
	// Expand right
	utils.ExpandMatRight(N, 1)
	utils.ExpandMatIntRight(backptrs, 1)

	// Initialize cell on top column
	(*N)[0][layer + 1] = -float64(layer + 1) * gap_penalty

	// Fill new rightmost column
	for row := 1; row < len(*N); row++ {
		nw_recurrence(N, backptrs, row, layer+1, q_idx + (row-1)*direction, d_idx + layer*direction, q, d, S_idx, gap_penalty, hit_thres)
	}
}

// Needleman-Wunsch recurrence
//
// Fill dynamic programming table and backpointers table at cell (i, j)
func nw_recurrence(N *[][]float64,
	                backptrs *[][]int,
					i, j, q_idx, d_idx int,
					q *string,
					d *[][]float64,
					S_idx *map[string]int,
					gap_penalty,
					hit_thres float64) {

	backptr := -1
	(*N)[i][j], backptr = utils.Max(
		(*N)[i-1][j-1] + sub_mat(q_idx, d_idx, q, d, S_idx, hit_thres),
		(*N)[i][j-1] - gap_penalty,
		(*N)[i-1][j] - gap_penalty,
	)
	(*backptrs)[i][j] = backptr
}

// Substitution matrix
//
// Return match score if probabiliy exceed threshold, else mismatch score
func sub_mat(q_idx, d_idx int,
			 q *string,
			 d *[][]float64,
			 S_idx *map[string]int,
			 hit_thres float64) float64 {
				 
	match_score := 1.
	mismatch_score := -1.

	char := string((*q)[q_idx])
	char_idx := (*S_idx)[char]
	
	if (*d)[char_idx][d_idx] >= hit_thres {
		return match_score
	
	} else {
		return mismatch_score
	}
}

func get_q_len_d_len(q_idx, d_idx int, q *string, d *[][]float64, direction int) (int, int) {
	var q_len int
	var d_len int

	if direction == 1 {
		q_len = len(*q) - q_idx
		d_len = len((*d)[0]) - d_idx
	} else {
		// direction == -1
		q_len = q_idx + 1
		d_len = d_idx + 1
	}
	return q_len, d_len
}

func preconditions(q_idx, d_idx int, q *string, d *[][]float64, direction int) bool {
	if direction != 1 && direction != -1 {
		// step has to be 1 or -1
		fmt.Println("wrong direction")
		return false
	}
	if q_idx < 0 || q_idx >= len(*q) || q_idx < 0 || q_idx >= len((*d)[0]) {
		// out of bounds
		fmt.Println("out of bounds")
		return false
	}
	return true
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func Left(q *string, d *[][]float64, hsp *[][]int, hit_thres, delta float64, substitution_matrix *[][]int) {
	q_idx := (*hsp)[0][0] - 1
	d_idx := (*hsp)[0][1] - 1
	extend(q_idx, d_idx, q, d, hit_thres, delta, substitution_matrix, -1)

	expand_by := 1

	// N := make([][]int, 0)

	N := utils.Mat(1, 1)
	for i := 0; i < 3; i++ {
		fmt.Println("i: ", i)
		utils.ExpandMat(N, expand_by)
		// fmt.Println(N)
		pretty.Print(N)
		// N = append(N, make([]int, 0))
		// fmt.Println(N)
	}
}

func Right(q *string, d *[][]float64, hsp *[][]int, hit_thres, delta float64, substitution_matrix *[][]int) {
	// q_idx := (*hsp)[1][0] - 1
	// d_idx := (*hsp)[1][1] - 1
	// extend(q_idx, d_idx, q, d, hit_thres, delta, substitution_matrix, 1)
}

func extend(q_idx, d_idx int, q *string, d *[][]float64, hit_thres, delta float64, substitution_matrix *[][]int, step int) {
	if step == 0 {
		// step has to be non-zero int
		return
	}
	if q_idx < 0 || q_idx >= len(*q) || q_idx < 0 || q_idx >= len(*d) {
		// out of bounds
		return
	}
}

func backtrace(N *[][]float64) {
	// const diag = 0
	// const left = 1
	// const top = 2
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

func SIdx(S *[]string) *map[string]int {
	S_idx := make(map[string]int)

	for i, char := range *S {
		S_idx[char] = i
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