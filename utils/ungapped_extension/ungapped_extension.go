package ungapped_extension

import (
	// "fmt"
	"math"
	"github.com/Mal-Jovi/561_Project/utils"
	"github.com/Mal-Jovi/561_Project/utils/gapped_extension"
	. "github.com/Mal-Jovi/561_Project/utils/structs"
)

func extend(q_idx, d_idx int, q *string, d *[][]float64, params *Params, direction int) (*[]int, float64) {
	has_exceeded_delta := false
	// max_score := math.Inf(-1)
	max_score := 0.
	max_q_idx := q_idx
	max_d_idx := d_idx
	score := 0.

	if direction >= 0 {
		direction = 1
	} else {
		direction = -1
	}

	S_idx := gapped_extension.SIdx(&params.S)

	for 0 < q_idx && q_idx < len(*q) - 1 && 0 < d_idx && d_idx < len((*d)[0]) - 1  {
		q_idx += direction
		d_idx += direction
		
		// score += prob if prob > hit_thres else score -= (hit_thres - prob)
		char := string((*q)[q_idx])
		char_idx := (*S_idx)[char]
		// prob := (*d)[ utils.SliceIndex(&params.S, string((*q)[q_idx])) ][ d_idx ]
		prob := (*d)[ char_idx ][ d_idx ]
		// fmt.Println("char:", char)
		if prob >= params.HitThres {
			score += prob
		} else {
			score -= params.HitThres - prob
		}

		if score >= max_score {
			max_score = score
			max_q_idx = q_idx
			max_d_idx = d_idx

		} else if math.Abs(max_score - score) > params.Delta {
			has_exceeded_delta = true
			break
		}
	}
	if has_exceeded_delta {
		return &[]int{max_q_idx, max_d_idx}, max_score
	}
	return &[]int{q_idx, d_idx}, score
}

func Left(q_idx, d_idx int, q *string, d *[][]float64, params *Params) (*[]int, float64) {
	return extend(q_idx, d_idx, q, d, params, -1)
}

func Right(q_idx, d_idx int, q *string, d *[][]float64, params *Params) (*[]int, float64) {
	return extend(q_idx + params.W - 1, d_idx + params.W - 1, q, d, params, 1)
}

func ScoreMiddle(q_idx, d_idx int, q *string, d *[][]float64, params *Params) float64 {
	score := 0.
	for i := 0; i < params.W; i++ {
		score += CharScore(q_idx+i, d_idx+i, q, d, params)
	}
	return score
}

func CharScore(q_idx, d_idx int, q *string, d *[][]float64, params *Params) float64 {
	char := string((*q)[q_idx])
	char_idx := utils.SliceIndex(&params.S, char)
	prob := (*d)[ char_idx ][ d_idx ]
	
	if prob >= params.HitThres {
		return prob
	} else {
		return - (params.HitThres - prob)
	}
}