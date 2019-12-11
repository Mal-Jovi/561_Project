package indexing

import (
	"math"
	"strings"
	"github.com/Mal-Jovi/561_Project/utils"
	. "github.com/Mal-Jovi/561_Project/utils/structs"
)

func ProbIndexTable(index *map[string][]int,
					start, end int,
					seeds *[][]string,
					d *[][]float64,
					params *Params) {

	for i := start; i < end; i++ {
		seed := strings.Join((*seeds)[i], "")
		indices := prob_find(&seed, d, params)
		(*index)[seed] = indices
	}
}

func prob_find(seed *string, d *[][]float64, params *Params) []int {

	indices := []int{}
	thres := float64(params.W) * math.Log(params.HitThres)

	for i := 0; i < len((*d)[0]) - params.W + 1; i++ {
		sum_log_prob := 0.
		for j := 0; j < params.W; j++ {
			sum_log_prob += math.Log( (* d)[ utils.SliceIndex(&params.S, string((*seed)[j])) ][i+j])
		}
		if sum_log_prob > thres {
			indices = append(indices, i)
		}
	}
	return indices
}

func GenSeeds(S *[]string, w int) [][] string {
	// seeds := make([][]string, len(S)^w)
	indexes := make([]int, w)
	seeds := [][]string{}

	for indexes != nil {
		seed := make([]string, w)
		for i, x := range indexes {
			seed[i] = (*S)[x]
		}

		for i := len(indexes) - 1; i >= 0; i-- {
			indexes[i]++
			if indexes[i] < len(*S) {
				break
			}
			indexes[i] = 0
			if i <= 0 {
				indexes = nil
				break
			}
		}
		seeds = append(seeds, seed)
	}
	return seeds
}