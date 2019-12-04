package main

import (
	"os"
	"fmt"
	"math"
	"strings"
	"github.com/Mal-Jovi/561_Project/utils"
	"github.com/Mal-Jovi/561_Project/utils/indexing"
	"github.com/Mal-Jovi/561_Project/utils/ungapped_extension"
	"github.com/Mal-Jovi/561_Project/utils/gapped_extension"
)

func main() {
	w := 7
	hit_thres := 0.9
	delta := 10.
	// delta := math.Inf(1)
	// hsp_thres := 0.
	// hsp_thres := 1.
	hsp_thres := 5.
	// hsp_thres := 0.5
	e_thres := 0.
	S := []string{"A", "T", "G", "C"}

	d := utils.GetProbSeq(
		"data/raw/chr22.maf.ancestors.42000000.complete.boreo.fa",
		utils.GetDConf("data/raw/chr22.maf.ancestors.42000000.complete.boreo.conf"),
		& S,
	)
	q := utils.SeqFromFasta("query.fa")

	prob_blast(q, d, w, & S, hit_thres, delta, hsp_thres, e_thres)
}

func prob_blast(q * string,
				d *[][] float64,
				w int,
				S *[] string,
		        hit_thres, delta, hsp_thres, e_thres float64) {

	// alingments := []int{}
	// index := prob_index_table(d, w, S, hit_thres)
	index := prob_index_table(d, w, S, hit_thres)	
	hsps := prob_extend(q, d, w, S, index, hit_thres, delta, hsp_thres, e_thres)
	fmt.Println(* hsps)
	prob_extend_gap(q, d, S, hsps, hit_thres, delta)
}

func prob_index_table(d *[][] float64, w int, S *[] string,
		              hit_thres float64) * map[string][]int {
	
	path := fmt.Sprint("data/processed/prob_index_table.w", w, ".hit_thres", hit_thres,".json")

	if _, err := os.Stat(path); !os.IsNotExist(err) {
		return utils.IndexFromJson(& path)
	}

	index := make(map[string] []int)
	seeds := indexing.GenSeeds(S, w)

	// Concurrently index database using 4 cores
	num_cores := 4
	for i := 0; i < num_cores - 1; i++ {
		go _prob_index_table(
			& index,
			i * (len(seeds) / num_cores), // Start index
			(i + 1) * (len(seeds) / num_cores), // End index
			& seeds, d, w, S, hit_thres,
		)
	}
	_prob_index_table(
		& index,
		(num_cores-1)*(len(seeds) / 4),
		len(seeds),
		& seeds, d, w, S, hit_thres,
	)
	utils.IndexToJson(& index, & path)

	return & index
}

func _prob_index_table(index * map[string][]int,
					   start, end int,
					   seeds *[][] string,
					   d *[][] float64,
					   w int,
					   S *[] string,
					   hit_thres float64)  {

	for i := start; i < end; i++ {
		seed := strings.Join((* seeds)[i], "")
		(* index)[seed] = indexing.ProbFind(& seed, d, w, S, hit_thres)
	}
}

func prob_extend(q *string,
				 d *[][]float64,
				 w int,
				 S *[]string,
				 index *map[string][]int,
				 hit_thres, delta, hsp_thres, e_thres float64) *[][][] int {

	hsps := [][][]int{}

	for q_idx := 0; q_idx < len(*q) - w + 1; q_idx++ {
		seed := (*q)[q_idx : q_idx + w]

		for _, d_idx := range (* index)[seed] {
			var hsp_score float64

			hsp_left, score_left := ungapped_extension.Left(q_idx, d_idx, q, d, w, S, hit_thres, delta)
			hsp_right, score_right := ungapped_extension.Right(q_idx, d_idx, q, d, w, S, hit_thres, delta)

			if score_left == math.Inf(-1) && score_right == math.Inf(-1) {
				fmt.Println("both")
				continue
			
			} else if score_left == math.Inf(-1) {
				fmt.Println("left")
				hsp_score = score_right
			
			} else if score_right == math.Inf(-1) {
				fmt.Println("right")
				hsp_score = score_left
			
			} else {
				fmt.Println("none")
				hsp_score = score_left + score_right
			}

			if hsp_score > hsp_thres {
				hsps = append(hsps, [][]int{*hsp_left, *hsp_right})
			}
		}
	}
	return &hsps
}

func prob_extend_gap(q *string, d *[][]float64, S *[]string, hsps *[][][]int, hit_thres, delta float64) {

	// substitution_matrix := gapped_extension.SubstitutionMatrix(S)
	S_idx := gapped_extension.SIdx(S)

	// fmt.Println(substitution_matrix)

	// for i, hsp := range *hsps {
	for i := 0; i < len(*hsps); i++ {
		// hsp := &(*hsps)[i]
		// fmt.Println(*hsp)
		// gapped_extension.NeedlemanWunsch(substitution_matrix)

		// gapped_extension.Left(q, d, &(*hsps)[i], hit_thres, delta, substitution_matrix)
		// gapped_extension.Right(q, d, &(*hsps)[i], hit_thres, delta, substitution_matrix)

		// gapped_extension.Extend(q, d, &(*hsps)[i], hit_thres, delta, substitution_matrix)
		// gapped_extension.Extend(len(*q) - 3, len((*d)[0]) - 6, q, d, S, S_idx, hit_thres, substitution_matrix, 1)
		
		gapped_extension.Extend(len(*q) - 3, len((*d)[0]) - 6, q, d, S, S_idx, hit_thres, 1)
		// gapped_extension.Extend(3, 6, q, d, S, S_idx, hit_thres, -1)
		// gapped_extension.Test()
		
		break
	}
	fmt.Println(len(*hsps))
}