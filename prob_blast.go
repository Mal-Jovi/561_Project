package main

import (
	"os"
	"fmt"
	"math"
	"sort"
	"github.com/kr/pretty"
	"github.com/Mal-Jovi/561_Project/utils"
	"github.com/Mal-Jovi/561_Project/utils/indexing"
	"github.com/Mal-Jovi/561_Project/utils/ungapped_extension"
	"github.com/Mal-Jovi/561_Project/utils/gapped_extension"
	. "github.com/Mal-Jovi/561_Project/utils/structs"
)

func main() {
	// Default config file path config.json
	config := "config.json"
	if len(os.Args[1:]) == 1 {
		// Parse command line argument for config file path
		config = os.Args[1]

	} else if len(os.Args[1:]) > 1 {
		// > 1 command line arguments are invalid. Stop
		fmt.Println("Invalid arguments. Usage: ./prob_blast [config file]")
		return
	}

	// Parse parameters from config file
	params := utils.ParamsFromJson(&config)
	pretty.Println(params)

	// Load query sequence q, probabilistic sequence d
	q := utils.SeqFromFasta(params.Q)
	d := utils.GetProbSeq(params.D, utils.GetDConf(params.DConf), &params.S)

	// Run probabilistic BLAST
	alignments := prob_blast(q, d, params)
	pretty.Println(*alignments)
	utils.SaveAlignments(alignments, params)
}

// prob_blast
//
// Probabilistic BLAST
//
func prob_blast(q * string, d *[][] float64, params *Params) *[]*Alignment {
	index := prob_index_table(d, params)
	hsps := prob_extend(q, d, index, params)
	alignments := prob_extend_gap(q, d, hsps, params)
	
	return alignments
}

// prob_index_table
//
// Index database, then export the index table
//
func prob_index_table(d *[][]float64, params *Params) *map[string][]int {
	path := fmt.Sprintf("data/processed/prob_index_table.w.%d.hit_thres.%.2f.json", params.W, params.HitThres)

	if _, err := os.Stat(path); !os.IsNotExist(err) {
		return utils.IndexFromJson(&path)
	}

	fmt.Println("Indexing database with word size", params.W, "and hit_thres", params.HitThres, "...")

	index := make(map[string][]int)
	seeds := indexing.GenSeeds(&params.S, params.W)

	// Concurrently index database using num_cores cores
	var start int
	var end int

	for i := 0; i < params.NumCores - 1; i++ {
		start = i * (len(seeds) / params.NumCores) // Start index
		end = (i + 1) * (len(seeds) / params.NumCores) // End index

		// Spawn go routine for concurrency
		go indexing.ProbIndexTable(&index, start, end, &seeds, d, params)
	}
	start = (params.NumCores - 1) * (len(seeds) / 4)
	end = len(seeds)
	indexing.ProbIndexTable(&index, start, end, &seeds, d, params)
	
	// Save index table to path
	utils.ExportToJson(&index, &path)
	fmt.Println("Saved index to", path)

	return &index
}

// prob_extend
//
// Ungapped extension. Return list of high scoring pairs
//
func prob_extend(q *string, d *[][]float64, index *map[string][]int, params *Params) *[]*Hsp {
	hsps := []*Hsp{}

	for q_idx := 0; q_idx < len(*q) - params.W + 1; q_idx++ {
		seed := (*q)[q_idx : q_idx + params.W]

		for _, d_idx := range (*index)[seed] {
			hsp_score := ungapped_extension.ScoreMiddle(q_idx, d_idx, q, d, params)

			// dummy := ""
			// dummy_S_idx := gapped_extension.SIdx(&params.S)

			// fmt.Println("seed hit acc:", hsp_score/utils.SumProbSeq(d_idx, d_idx+params.W-1, &dummy, d, dummy_S_idx, params.HitThres, params))

			hsp_left, score_left := ungapped_extension.Left(q_idx, d_idx, q, d, params)
			hsp_right, score_right := ungapped_extension.Right(q_idx, d_idx, q, d, params)

			if score_left == math.Inf(-1) && score_right == math.Inf(-1) {
				continue
			}
			if score_left > math.Inf(-1) {
				hsp_score += score_left
			}
			if score_right > math.Inf(-1) {
				hsp_score += score_right
			}

			// E-value of HSP
			e_val := 0.711 * float64(len(*q) * len((*d)[0])) * math.Exp(-1.37 * hsp_score)

			// if hsp_score >= params.HspThres {
			if e_val <= params.EThres {
				var hsp Hsp
				hsp.QIdxLeft = (*hsp_left)[0]
				hsp.QIdxRight = (*hsp_right)[0]
				hsp.DIdxLeft = (*hsp_left)[1]
				hsp.DIdxRight = (*hsp_right)[1]
				hsp.Score = hsp_score
				hsp.EVal = e_val
				hsps = append(hsps, &hsp)

				// hsp_acc := hsp.Score/utils.SumProbSeq(hsp.DIdxLeft, hsp.DIdxRight, &dummy, d, dummy_S_idx, params.HitThres, params)
			}
		}
	}
	return &hsps
}

// prob_extend_gap
//
// Gapped extension. Further extend HSPs using Needleman-Wunsch
// Return list of resulting alignments
//
func prob_extend_gap(q *string, d *[][]float64, hsps *[]*Hsp, params *Params) *[]*Alignment {
	S_idx := gapped_extension.SIdx(&params.S)
	alignments := []*Alignment{}

	for i := 0; i < len(*hsps); i++ {
		alignment := gapped_extension.Extend((*hsps)[i], q, d, S_idx, params)
		alignments = append(alignments, &alignment)
	}
	sort.Slice(alignments, func(i, j int) bool {
		// return alignments[i].Accuracy > alignments[j].Accuracy
		return alignments[i].EVal < alignments[j].EVal
	})
	return &alignments
}