package ungapped_extension

import (
	"math"
	"github.com/Mal-Jovi/561_Project/utils"
)

func extend(q_idx, d_idx int,
			q *string, d *[][]float64, w int, S *[]string,
			hit_thres, delta float64,
			step int) (*[]int, float64) {

	has_exceeded_delta := false
	max_score := math.Inf(-1)
	max_q_idx := q_idx
	max_d_idx := d_idx
	score := 0.

	for 0 < q_idx && q_idx < len(* q) - 1 && 0 < d_idx && d_idx < len((* d)[0]) - 1  {
		q_idx += step
		d_idx += step

		// add diff between prob and hit_thres
		// score += (* d)[ utils.SliceIndex(S, string((* q)[q_idx])) ][d_idx] - hit_thres
		
		// score = score + prob if prob > hit_thres else score - (hit_thres - prob)
		prob := (*d)[ utils.SliceIndex(S, string((*q)[q_idx])) ][ d_idx ]
		if prob > hit_thres {
			score += prob
		} else {
			score -= hit_thres - prob
		}

		if score >= max_score {
			max_score = score
			max_q_idx = q_idx
			max_d_idx = d_idx

		} else if math.Abs(max_score - score) > delta {
			has_exceeded_delta = true
			break
		}
	}
	if has_exceeded_delta {
		return &[]int{max_q_idx, max_d_idx}, max_score
	}
	return &[]int{q_idx, d_idx}, score
// return & []int{max_q_idx, max_d_idx}, score
}

func Left(q_idx, d_idx int,
		  q *string,
		  d *[][]float64,
		  w int, S *[]string,
		  hit_thres, delta float64) (*[]int, float64) {

	return extend(q_idx, d_idx, q, d, w, S, hit_thres, delta, -1)
}

func Right(q_idx, d_idx int,
		   q *string,
		   d *[][]float64,
		   w int,
		   S *[]string,
		   hit_thres, delta float64) (*[]int, float64) {

	return extend(q_idx + w - 1, d_idx + w - 1, q, d, w, S, hit_thres, delta, 1)
}