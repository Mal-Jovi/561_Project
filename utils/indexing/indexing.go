package indexing

import (
	"math"
	"github.com/Mal-Jovi/561_Project/utils"
)

func GenSeeds(S *[] string, w int) [][] string {
	// seeds := make([][]string, len(S)^w)
	indexes := make([]int, w)
	seeds := [][]string{}

	for indexes != nil {
		seed := make([]string, w)
		for i, x := range indexes {
			seed[i] = (* S)[x]
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

func ProbFind(seed * string, d *[][] float64, w int,
	S *[] string, hit_thres float64) [] int {

	indices := []int{}
	thres := float64(w) * math.Log(hit_thres)

	for i := 0; i < len((*d)[0]) - w + 1; i++ {
		sum_log_prob := 0.
		for j := 0; j < w; j++ {
			sum_log_prob += math.Log( (* d)[ utils.SliceIndex(S, string((*seed)[j])) ][i+j])
		}

		if sum_log_prob > thres {
			indices = append(indices, i)
		}
	}
	return indices
}