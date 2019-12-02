package ungapped_extension

import (
	"math"
	"github.com/Mal-Jovi/561_Project/utils"
)

func ProbLeft(q_idx, d_idx int, q * string, d *[][] float64, w int, S *[] string,
	           hit_thres, delta float64) (* []int, float64) {

	return _prob_extend(q_idx, d_idx, q, d, w, S, hit_thres, delta, -1)
}

func ProbRight(q_idx, d_idx int, q * string, d *[][] float64, w int, S *[] string,
	            hit_thres, delta float64) (* []int, float64) {

	return _prob_extend(q_idx, d_idx, q, d, w, S, hit_thres, delta, 1)
}

func _prob_extend(q_idx, d_idx int, q * string, d *[][] float64, w int, S *[] string,
	  hit_thres, delta float64, step int) (* []int, float64) {

	max_score := math.Inf(-1)
	max_q_idx := q_idx
	max_d_idx := d_idx
	score := 0.

	for 0 < q_idx && q_idx < len(* q) - 1 && 0 < d_idx && d_idx < len((* d)[0]) - 1  {
		q_idx += step
		d_idx += step

		score += (* d)[ utils.SliceIndex(S, string((* q)[q_idx])) ][d_idx] - hit_thres
		if score > max_score {
			max_score = score
			max_q_idx = q_idx
			max_d_idx = d_idx

		} else if max_score - score > delta {
			break
		}
	}
	return & []int{max_q_idx, max_d_idx}, max_score
	// return & []int{max_q_idx, max_d_idx}, score
}

