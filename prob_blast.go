package main

import (
	"os"
	"fmt"
	"math"
	"strings"
	"github.com/Mal-Jovi/561_Project/utils"
	"github.com/Mal-Jovi/561_Project/blast/index"
	"github.com/Mal-Jovi/561_Project/blast/ungapped_extension"
)

func main() {
	w := 4
	hit_thres := 0.9
	delta := 10.
	// delta := math.Inf(1)
	hsp_thres := 0.
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
	prob_extend_gap(q, d, w, S, hsps)
}

func prob_index_table(d *[][] float64, w int, S *[] string,
		              hit_thres float64) * map[string][]int {
	
	path := fmt.Sprint("data/processed/prob_index_table.w", w, ".hit_thres", hit_thres,".json")

	if _, err := os.Stat(path); !os.IsNotExist(err) {
		return utils.IndexFromJson(& path)
	}

	index := make(map[string] []int)
	seeds := index.GenSeeds(S, w)

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
		(* index)[seed] = index.ProbFind(& seed, d, w, S, hit_thres)
	}
}

func prob_extend(q * string, d *[][] float64, w int, S *[] string,
                 index * map[string][]int, hit_thres, delta, hsp_thres, e_thres float64) *[][][] int {

	hsps := [][][]int{}

	for q_idx := 0; q_idx < len(* q) - w + 1; q_idx++ {
		seed := (* q)[q_idx:q_idx+w]

		for _, d_idx := range (* index)[seed] {
			var hsp_score float64

			hsp_left, score_left := ungapped_extension.ProbLeft(q_idx, d_idx, q, d, w, S, hit_thres, delta)
			hsp_right, score_right := ungapped_extension.ProbRight(q_idx, d_idx, q, d, w, S, hit_thres, delta)

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
				hsps = append(hsps, [][]int{* hsp_left, * hsp_right})
			}
		}
	}
	return & hsps
}

func prob_extend_gap(q * string, d *[][] float64, w int, S *[] string, hsps *[][][] int) {

	// substitution_matrix := substitution_matrix(S)
	// needleman_wunsch(substitution_matrix)
}

func needleman_wunsch(seq1, seq2 * string,
	                  S *[] string,
					  gap_open_penaty, gap_extend_penalty int,
					  substitution_matrix *[][] int) {

	N := init_nw(seq1, seq2, gap_open_penaty, gap_extend_penalty)
	S_idx := S_idx(S)

	for i := 1; i < len(* seq1); i++ {
		for j := 1; j < len(* seq2); j++ {
			(* N)[i][j] = utils.Max(
				(* N)[i-1][j-1] + float64((* substitution_matrix)[(* S_idx)[string((* seq1)[i])]][(* S_idx)[string((* seq2)[j])]]),
				(* N)[i][j-1] + float64(gap_extend_penalty),
				(* N)[i-1][j] + float64(gap_extend_penalty),
			)
		}
	}
}

func S_idx(S *[] string) * map[string]int {
	S_idx := make(map[string]int)

	for i, word := range * S {
		S_idx[word] = i
	}
	return & S_idx
}

func init_nw(seq1, seq2 * string, gap_open_penaty, gap_extend_penalty int) *[][] float64 {
	N := utils.Mat(len(* seq1), len(* seq2))
	
	
	return N
}

func substitution_matrix(S *[] string) *[][] int {
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

func backtrace(N *[][] float64) {

}