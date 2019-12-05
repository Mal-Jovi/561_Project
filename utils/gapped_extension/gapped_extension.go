package gapped_extension

import (
	"fmt"
	"math"
	// "github.com/kr/pretty"
	"github.com/Mal-Jovi/561_Project/utils"
	. "github.com/Mal-Jovi/561_Project/utils/structs"
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
func extend(q_idx, d_idx int, q *string, d *[][]float64, S_idx *map[string]int, params *Params, direction int) (*string, *string, float64) {

	var q_aligned *string
	var d_aligned *string
	
	if !preconditions(q_idx, d_idx, q, d) {
		// fmt.Println("SAD")
		return q_aligned, d_aligned, 0.
	}

	if direction >= 0 {
		direction = 1
	} else {
		direction = -1
	}
	
	gap_penalty := 1.
	
	N := utils.Mat(1, 1)
	backptrs := utils.MatInt(1, 1)
	
	q_len, d_len := q_len_d_len(q_idx, d_idx, q, d, direction)

	if direction >= 0 {
		fmt.Println("extending right")
		// fmt.Println("q seg:", (*q)[q_idx:q_idx + q_len])
		// fmt.Println("d seg:", *utils.PrettyProbSeg(d, d_idx, d_idx + d_len, S, hit_thres))
	} else {
		fmt.Println("extending left")
		// fmt.Println("q seg:", (*q)[:q_idx + 1])
		// fmt.Println("d seq:", *utils.PrettyProbSeg(d, 0, d_idx + 1, &params.S, params.HitThres))
	}

	i_end, j_end := fill(N, backptrs, q_idx, d_idx, q_len, d_len, q, d, S_idx, params, gap_penalty, direction)
	fmt.Println("exited fill()")
	q_aligned, d_aligned = traceback(backptrs, i_end, j_end, q_idx, d_idx, q, d, &params.S, params.HitThres, direction)
	// traceback(backptrs, i_end, j_end, q_idx, d_idx, q, d, S, hit_thres, direction)
	
	// fmt.Println("NW table:")
	// pretty.Println(*N)
	// fmt.Println("backpointers:")
	// pretty.Println(*backptrs)
	// fmt.Println("alignment:")
	// fmt.Println(*q_aligned)
	// fmt.Println(*d_aligned)

	return q_aligned, d_aligned, (*N)[i_end][j_end]
}

// func Extend(hsp *[][]int, q *string, d *[][]float64, S *[]string, S_idx *map[string]int, hit_thres, delta float64) Alignment {
func Extend(hsp *Hsp, q *string, d *[][]float64, S_idx *map[string]int, params *Params) Alignment {
	var alignment Alignment

	left_q_aligned, left_d_aligned, left_score := left(hsp, q, d, S_idx, params)
	right_q_aligned, right_d_aligned, right_score := right(hsp, q, d, S_idx, params)
	q_aligned, d_aligned := middle(hsp, q, d, &params.S)

	if left_q_aligned == nil || left_d_aligned == nil {
		fmt.Println("didn't extend left")
	} else {
		q_aligned = *left_q_aligned + q_aligned
		*d_aligned = *left_d_aligned + *d_aligned
	}

	if right_q_aligned == nil || right_d_aligned == nil {
		fmt.Println("didn't extend right")
	} else {
		q_aligned += *right_q_aligned
		*d_aligned += *right_d_aligned
	}

	alignment.QAligned = q_aligned
	alignment.DAligned = *d_aligned
	alignment.Hsp = *hsp
	alignment.Score = left_score + right_score + hsp.Score
	// return &q_aligned, d_aligned
	return alignment
}

// func left(hsp *[][]int, q *string, d *[][]float64, S *[]string, S_idx *map[string]int, hit_thres, delta float64) (*string, *string) {
func left(hsp *Hsp, q *string, d *[][]float64, S_idx *map[string]int, params *Params) (*string, *string, float64) {
	// q_idx := (*hsp)[0][0] - 1
	// d_idx := (*hsp)[0][1] - 1
	q_idx := hsp.QIdxLeft - 1
	d_idx := hsp.DIdxLeft - 1
	// extend(q_idx, d_idx, q, d, hit_thres, delta, substitution_matrix, -1)
	return extend(q_idx, d_idx, q, d, S_idx, params, -1)
}

func right(hsp *Hsp, q *string, d *[][]float64, S_idx *map[string]int, params *Params) (*string, *string, float64) {
	// q_idx := (*hsp)[1][0] + 1
	// d_idx := (*hsp)[1][1] + 1
	q_idx := hsp.QIdxRight + 1
	d_idx := hsp.DIdxRight + 1
	return extend(q_idx, d_idx, q, d, S_idx, params, 1)
}

// func middle(hsp *[][]int, q *string, d *[][]float64, S *[]string) (string, *string) {
func middle(hsp *Hsp, q *string, d *[][]float64, S *[]string) (string, *string) {
	// return (*q)[(*hsp)[0][0] : (*hsp)[1][0] + 1], utils.PrettyProbSeg(d, (*hsp)[0][1], (*hsp)[1][1] + 1, S, 0.)
	return (*q)[hsp.QIdxLeft : hsp.QIdxRight + 1], utils.PrettyProbSeg(d, hsp.DIdxLeft, hsp.DIdxRight + 1, S, 0.)
}

func fill(N *[][]float64,
		  backptrs *[][]int,
		  q_idx, d_idx, q_len, d_len int,
		  q *string,
		  d *[][]float64,
		  S_idx *map[string]int,
		  params *Params,
		  gap_penalty float64,
		  direction int) (int, int) {
	
	var i_end int
	var j_end int

	has_exceeded_delta := false
	i_max, j_max := -1, -1

	// Incrementally expand dynamic programming and backpointers tables by layer
	for layer := 0; layer < max(q_len, d_len); layer++ {
		bottom_i_max, bottom_j_max := -1, -1
		right_i_max, right_j_max := -1, -1

		// Expand, initialize top row and/or left column, then fill cells on layer
		if layer < q_len {
			bottom_i_max, bottom_j_max = fill_bottom(N, backptrs, layer, q_idx, d_idx, q, d, S_idx, gap_penalty, params.HitThres, direction)
		}
		if layer < d_len {
			right_i_max, right_j_max = fill_right(N, backptrs, layer, q_idx, d_idx, q, d, S_idx, gap_penalty, params.HitThres, direction)
		}
		// pretty.Println(*N)

		layer_i_max, layer_j_max := layer_cell_max(N, bottom_i_max, bottom_j_max, right_i_max, right_j_max)

		// i_max, j_max = cell_max(N, i_max, j_max, layer_i_max, layer_j_max)
		if i_max < 0 || j_max < 0 || (*N)[layer_i_max][layer_j_max] >= (*N)[i_max][j_max] {
			i_max = layer_i_max
			j_max = layer_j_max
			
		} else if math.Abs((*N)[i_max][j_max] - (*N)[layer_i_max][layer_j_max]) > params.Delta {
			fmt.Println("exceeded delta")
			has_exceeded_delta = true
			break
		}
	}
	if has_exceeded_delta {
		i_end = i_max
		j_end = j_max
		// fmt.Println("i_max:", i_max)
		// fmt.Println("j_max:", j_max)
	} else {
		i_end = q_len
		j_end = d_len
		// fmt.Println("last cell:", (*N)[q_len][d_len])
		// fmt.Println("q_len:", q_len)
		// fmt.Println("d_len:", d_len)
	}
	return i_end, j_end
}

// func cell_max(N *[][]float64, layer_i_max, layer_j_max int) (int, int) {
// 	if i_max == -1 && j_max == -1 {
// 		return layer_i_max, layer_j_max
// 	}
// }

func layer_cell_max(N *[][]float64, bottom_i_max, bottom_j_max, right_i_max, right_j_max int) (int, int) {
	// fmt.Println("bottom_i_max:", bottom_i_max)
	// fmt.Println("bottom_j_max:", bottom_j_max)
	// fmt.Println("right_i_max:", right_i_max)
	// fmt.Println("right_j_max:", right_j_max)

	if (bottom_i_max < 0 || bottom_j_max < 0) && right_i_max >= 0 && right_j_max >= 0 {
		return right_i_max, right_j_max
	}
	if (right_i_max < 0 || right_j_max < 0) && bottom_i_max >= 0 && bottom_j_max >= 0 {
		return bottom_i_max, bottom_j_max
	}
	if (right_i_max < 0 || right_j_max < 0) && (bottom_i_max < 0 || bottom_j_max < 0) {
		panic("[layer_cell_max] Index out of range: -1 present in both right and bottom indices")
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
	j_max := 0
	max_val := math.Inf(-1)

	// Fill new bottom row
	for col := 1; col < len((*N)[0]); col++ {
		cell_val := nw_recurrence(
			N, backptrs,
			layer+1, col,
			q_idx + layer*direction,
			d_idx + (col-1)*direction,
			q, d,
			S_idx,
			gap_penalty,
			hit_thres,
		)
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

	i_max := 0
	j_max := layer + 1
	max_val := math.Inf(-1)

	// Fill new rightmost column
	for row := 1; row < len(*N); row++ {
		cell_val := nw_recurrence(
			N, backptrs,
			row, layer+1,
			q_idx + (row-1)*direction,
			d_idx + layer*direction,
			q, d,
			S_idx,
			gap_penalty,
			hit_thres,
		)
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
				 
	// match_score := 1.
	// mismatch_score := -1.

	char := string((*q)[q_idx])
	char_idx := (*S_idx)[char]
	
	if (*d)[char_idx][d_idx] >= hit_thres {
		// return match_score
		return (*d)[char_idx][d_idx]
	
	} else {
		// return mismatch_score
		return -(hit_thres - (*d)[char_idx][d_idx])
	}
}

func q_len_d_len(q_idx, d_idx int, q *string, d *[][]float64, direction int) (int, int) {
	var q_len int
	var d_len int

	if direction >= 0 {
		q_len = len(*q) - q_idx
		d_len = len((*d)[0]) - d_idx
	} else {
		// direction < 0
		q_len = q_idx + 1
		d_len = d_idx + 1
	}
	return q_len, d_len
}

func preconditions(q_idx, d_idx int, q *string, d *[][]float64) bool {
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

// traceback
//
// Recover Needleman-Wunsch alignment
func traceback(backptrs *[][]int,
			   i_end, j_end,
			   q_idx, d_idx int,
			   q *string,
			   d *[][]float64,
			   S *[]string,
			   hit_thres float64,
			   direction int) (*string, *string) {

	const diag = 0
	const left = 1
	const top = 2

	q_aligned, d_aligned := "", ""

	var d_seg *string

	if direction > 0 {
		d_seg = utils.PrettyProbSeg(d, d_idx, d_idx + j_end, S, 0.)
	} else {
		d_seg = utils.PrettyProbSeg(d, d_idx - j_end + 1, d_idx + 1, S, 0.)
	}

	// fmt.Println("d_seg:", *d_seg)

	// q_idx + i_end*direction
	// d_idx + j_end*direction

	if direction >= 0 {
		q_idx += i_end - 1
		d_idx = j_end - 1
	} else {
		q_idx -= i_end - 1
		d_idx = 0
	}

	for i_end > 0 && j_end > 0 {
		switch (*backptrs)[i_end][j_end] {
		case diag:
			traceback_diag(q_idx, d_idx, &q_aligned, &d_aligned, q, d_seg, direction)
			i_end--; j_end--
			q_idx -= direction
			d_idx -= direction
		case left:
			traceback_left(q_idx, d_idx, &q_aligned, &d_aligned, q, d_seg, direction)
			j_end--
			d_idx -= direction
		case top:
			traceback_top(q_idx, d_idx, &q_aligned, &d_aligned, q, d_seg, direction)
			i_end--
			q_idx -= direction
		}
	}
	for ; j_end > 0; j_end-- {
		traceback_left(q_idx, d_idx, &q_aligned, &d_aligned, q, d_seg, direction)
		d_idx -= direction
	}
	for ; i_end > 0; i_end-- {
		traceback_top(q_idx, d_idx, &q_aligned, &d_aligned, q, d_seg, direction)
		q_idx -= direction
	}
	return &q_aligned, &d_aligned
}

func traceback_diag(q_idx, d_idx int, q_aligned, d_aligned, q, d_seg *string, direction int) {
	if direction >= 0 {
		*q_aligned = string((*q)[q_idx]) + *q_aligned
		*d_aligned = string((*d_seg)[d_idx]) + *d_aligned
	} else {
		*q_aligned += string((*q)[q_idx])
		*d_aligned += string((*d_seg)[d_idx])
	}
}

func traceback_left(q_idx, d_idx int, q_aligned, d_aligned, q, d_seg *string, direction int) {
	if direction >= 0 {
		*q_aligned = "-" + *q_aligned
		*d_aligned = string((*d_seg)[d_idx]) + *d_aligned
	} else {
		*q_aligned += "-"
		*d_aligned += string((*d_seg)[d_idx])
	}
}

func traceback_top(q_idx, d_idx int, q_aligned, d_aligned, q, d_seg *string, direction int) {
	if direction >= 0 {
		*q_aligned = string((*q)[q_idx]) + *q_aligned
		*d_aligned = "-" + *d_aligned
	} else {
		*q_aligned += string((*q)[q_idx])
		*d_aligned += "-"
	}
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