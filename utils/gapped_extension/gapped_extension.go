package gapped_extension

import (
	"fmt"
	"math"
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
			hit_thres,
			delta float64,
			// substitution_matrix *[][]int,
			direction int) {

	if !preconditions(q_idx, d_idx, q, d, direction) {
		fmt.Println("SAD")
		return
	}
	
	gap_penalty := 1.
	
	N := utils.Mat(1, 1)
	backptrs := utils.MatInt(1, 1)
	
	q_len, d_len := q_len_d_len(q_idx, d_idx, q, d, direction)

	if direction == 1 {
		fmt.Println("q seg:", (*q)[q_idx:q_idx + q_len])
		fmt.Println("d seg:", *utils.PrettyProbSeg(d, d_idx, d_idx + d_len, S, hit_thres))
	} else {
		fmt.Println("extending left")
		fmt.Println("q seg:", (*q)[:q_idx + 1])
		fmt.Println("d seq:", *utils.PrettyProbSeg(d, 0, d_idx + 1, S, hit_thres))
	}

	i_end, j_end := fill(N, backptrs, q_idx, d_idx, q_len, d_len, q, d, S_idx, gap_penalty, hit_thres, delta, direction)
	// alignment := backtrace(backptrs, i_end, j_end, q_idx, d_idx, q, d, direction)
	backtrace(backptrs, i_end, j_end, q_idx, d_idx, q, d, direction)
	
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
		  hit_thres,
		  delta float64,
		  direction int) (int, int) {
	
	has_exceeded_delta := false
	i_max, j_max := -1, -1

	// Incrementally expand dynamic programming and backpointers tables by layer
	for layer := 0; layer < max(q_len, d_len); layer++ {
		bottom_i_max, bottom_j_max := -1, -1
		right_i_max, right_j_max := -1, -1

		// Expand, then initialize top row and/or left column
		if layer < q_len {
			bottom_i_max, bottom_j_max = fill_bottom(N, backptrs, layer, q_idx, d_idx, q, d, S_idx, gap_penalty, hit_thres, direction)
		}
		if layer < d_len {
			right_i_max, right_j_max = fill_right(N, backptrs, layer, q_idx, d_idx, q, d, S_idx, gap_penalty, hit_thres, direction)
		}
		layer_i_max, layer_j_max := layer_cell_max(N, bottom_i_max, bottom_j_max, right_i_max, right_j_max)

		// i_max, j_max = cell_max(N, i_max, j_max, layer_i_max, layer_j_max)
		if i_max < 0 || j_max < 0 || (*N)[layer_i_max][layer_j_max] >= (*N)[i_max][j_max] {
			i_max = layer_i_max
			j_max = layer_j_max
			
		} else if math.Abs((*N)[i_max][j_max] - (*N)[layer_i_max][layer_j_max]) > delta {
			fmt.Println("exceeded delta")
			has_exceeded_delta = true
			break
		}
	}
	fmt.Println("i_max:", i_max)
	fmt.Println("j_max:", j_max)
	if has_exceeded_delta {
		return i_max, j_max
	}
	fmt.Println("last cell:", (*N)[q_len][d_len])
	fmt.Println("q_len:", q_len)
	fmt.Println("d_len:", d_len)
	return q_len, d_len
}

// func cell_max(N *[][]float64, layer_i_max, layer_j_max int) (int, int) {
// 	if i_max == -1 && j_max == -1 {
// 		return layer_i_max, layer_j_max
// 	}
// }

func layer_cell_max(N *[][]float64, bottom_i_max, bottom_j_max, right_i_max, right_j_max int) (int, int) {
	if (bottom_i_max < 0 || bottom_j_max < 0) && right_i_max >= 0 && right_j_max >= 0 {
		return right_i_max, right_j_max
	}
	if (right_i_max < 0 || right_j_max < 0) && bottom_i_max >= 0 && bottom_j_max >= 0 {
		return bottom_i_max, bottom_j_max
	}
	if (right_i_max < 0 || right_j_max < 0) && (bottom_i_max < 0 || bottom_j_max < 0) {
		panic("[layer_cell_max] Index out of range")
	}

	if (*N)[bottom_i_max][bottom_j_max] > (*N)[right_i_max][right_j_max] {
		return bottom_i_max, bottom_j_max
	}
	return right_i_max, right_j_max
}

func fill_bottom(N *[][]float64,
				 backptrs *[][]int,
				 layer, q_idx, d_idx int,
				 q *string,
				 d *[][]float64,
				 S_idx *map[string]int,
				 gap_penalty,
				 hit_thres float64,
				 direction int) (int, int) {

	// Expand bottom
	utils.ExpandMatBottom(N, 1)
	utils.ExpandMatIntBottom(backptrs, 1)

	// Initialize cell on leftmost column
	(*N)[layer + 1][0] = -float64(layer + 1) * gap_penalty

	i_max := layer + 1
	j_max := -1
	max_val := math.Inf(-1)

	// Fill new bottom row
	for col := 1; col < len((*N)[0]); col++ {
		cell_val := nw_recurrence(N, backptrs, layer+1, col, q_idx + layer*direction, d_idx + (col-1)*direction, q, d, S_idx, gap_penalty, hit_thres)

		if cell_val > max_val {
			max_val = cell_val
			j_max = col
		}
	}
	return i_max, j_max
}

func fill_right(N *[][]float64,
				backptrs *[][]int,
				layer,
				q_idx, d_idx int,
				q *string,
				d *[][]float64,
				S_idx *map[string]int,
				gap_penalty,
				hit_thres float64,
				direction int) (int, int) {
	
	// Expand right
	utils.ExpandMatRight(N, 1)
	utils.ExpandMatIntRight(backptrs, 1)

	// Initialize cell on top column
	(*N)[0][layer + 1] = -float64(layer + 1) * gap_penalty

	i_max := -1
	j_max := layer + 1
	max_val := math.Inf(-1)

	// Fill new rightmost column
	for row := 1; row < len(*N); row++ {
		cell_val := nw_recurrence(N, backptrs, row, layer+1, q_idx + (row-1)*direction, d_idx + layer*direction, q, d, S_idx, gap_penalty, hit_thres)

		if cell_val > max_val {
			max_val = cell_val
			i_max = row
		}
	}
	return i_max, j_max
}

// Needleman-Wunsch recurrence
//
// Fill dynamic programming table and backpointers table at cell (i, j)
func nw_recurrence(N *[][]float64,
				   backptrs *[][]int,
				   i, j,
				   q_idx, d_idx int,
				   q *string,
				   d *[][]float64,
				   S_idx *map[string]int,
				   gap_penalty,
				   hit_thres float64) float64 {

	// backptr := -1
	
	cell_val, backptr := utils.Max(
		(*N)[i-1][j-1] + sub_mat(q_idx, d_idx, q, d, S_idx, hit_thres),
		(*N)[i][j-1] - gap_penalty,
		(*N)[i-1][j] - gap_penalty,
	)
	(*N)[i][j] = cell_val
	(*backptrs)[i][j] = backptr

	return cell_val
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
		// return (*d)[char_idx][d_idx]
	
	} else {
		return mismatch_score
		// return -(hit_thres - (*d)[char_idx][d_idx])
	}
}

func q_len_d_len(q_idx, d_idx int, q *string, d *[][]float64, direction int) (int, int) {
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

// func Left(q *string, d *[][]float64, hsp *[][]int, hit_thres, delta float64, substitution_matrix *[][]int) {
func Left(hsp *[][]int, q *string, d *[][]float64, S *[]string, S_idx *map[string]int, hit_thres, delta float64) {
	q_idx := (*hsp)[0][0] - 1
	d_idx := (*hsp)[0][1] - 1
	// extend(q_idx, d_idx, q, d, hit_thres, delta, substitution_matrix, -1)
	Extend(q_idx, d_idx, q, d, S, S_idx, hit_thres, delta, -1)
}

func Right(hsp *[][]int, q *string, d *[][]float64, S *[]string, S_idx *map[string]int, hit_thres, delta float64) {
	q_idx := (*hsp)[1][0] + 1
	d_idx := (*hsp)[1][1] + 1
	Extend(q_idx, d_idx, q, d, S, S_idx, hit_thres, delta, 1)
}

// backtrace
//
// Recover Needleman-Wunsch alignment
func backtrace(backptrs *[][]int,
			   i_end, j_end,
			   q_idx, d_idx int,
			   q *string,
			   d *[][]float64,
			   direcion int) *[]string {

	// const diag = 0
	// const left = 1
	// const top = 2
	q_prime, d_prime := "", ""
	return &[]string{q_prime, d_prime}
}

// func SubstitutionMatrix(S *[]string) *[][]int {
// 	m := utils.MatInt(len(* S), len(* S))
	
// 	match_score := 1
// 	mismatch_score := -1

// 	for i := 0; i < len(* S); i++ {
// 		for j := 0; j < len(* S); j++ {
// 			if i == j {
// 				(* m)[i][j] = match_score
			
// 			} else {
// 				(* m)[i][j] = mismatch_score
// 			}
// 		}
// 	}
// 	return m
// }

func SIdx(S *[]string) *map[string]int {
	S_idx := make(map[string]int)

	for i, char := range *S {
		S_idx[char] = i
	}
	return &S_idx
}

// func init_nw(seq1, seq2 *string, gap_open_penaty, gap_extend_penalty int) *[][]float64 {
// 	N := utils.Mat(len(* seq1), len(* seq2))

// 	// fill first row and first column here
// 	for i := 0; i < len(*seq1); i++ {
// 		(*N)[i][0] = float64(i * gap_extend_penalty)
// 	}
// 	for j := 1; j < len(*seq2); j++ {
// 		(*N)[0][j] = float64(j * gap_extend_penalty)
// 	}
// 	return N
// }