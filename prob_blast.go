package main

import (
	"os"
	"fmt"
	"math"
	"strings"
	"github.com/kr/pretty"
	"github.com/Mal-Jovi/561_Project/utils"
	"github.com/Mal-Jovi/561_Project/utils/indexing"
	"github.com/Mal-Jovi/561_Project/utils/ungapped_extension"
	"github.com/Mal-Jovi/561_Project/utils/gapped_extension"
	. "github.com/Mal-Jovi/561_Project/utils/structs"
)

func main() {
	var params Params
	params.W = 4 // Word size
	params.HitThres = 0.9
	params.Delta = 2.5
	// params.HspThres = 5.
	params.HspThres = 11.
	params.S = []string{"A", "T", "G", "C"} // Alphabet
	params.EThres = 0.

	d := utils.GetProbSeq(
		"data/raw/chr22.maf.ancestors.42000000.complete.boreo.fa",
		utils.GetDConf("data/raw/chr22.maf.ancestors.42000000.complete.boreo.conf"),
		&params.S,
	)
	q := utils.SeqFromFasta("query.fa")

	// prob_blast(q, d, w, &S, hit_thres, delta, hsp_thres, e_thres)
	alignments := prob_blast(q, d, &params)
	utils.SaveAlignments(alignments, &params)
}

func prob_blast(q * string, d *[][] float64, params *Params) *[]*Alignment {

	// alingments := []int{}
	// index := prob_index_table(d, w, S, hit_thres)
	index := prob_index_table(d, params)
	// hsps := prob_extend(q, d, w, S, index, hit_thres, delta, hsp_thres, e_thres)
	hsps := prob_extend(q, d, index, params)
	// pretty.Println(*hsps)
	// prob_extend_gap(q, d, S, hsps, hit_thres, delta)
	alignments := prob_extend_gap(q, d, hsps, params)
	pretty.Println(*alignments)
	return alignments
}

func prob_index_table(d *[][]float64, params *Params) *map[string][]int {
	
	path := fmt.Sprintf("data/processed/prob_index_table.w%d.hit_thres%.2f.json", params.W, params.HitThres)

	fmt.Println(path)

	if _, err := os.Stat(path); !os.IsNotExist(err) {
		return utils.IndexFromJson(&path)
	}

	fmt.Println("Indexing database with word size", params.W, "and hit_thres", params.HitThres, "...")

	index := make(map[string][]int)
	seeds := indexing.GenSeeds(&params.S, params.W)

	// Concurrently index database using 4 cores
	num_cores := 4
	for i := 0; i < num_cores - 1; i++ {
		go _prob_index_table(
			&index,
			i * (len(seeds) / num_cores), // Start index
			(i + 1) * (len(seeds) / num_cores), // End index
			&seeds, d,
			params,
		)
	}
	_prob_index_table(
		&index,
		(num_cores-1)*(len(seeds) / 4),
		len(seeds),
		&seeds, d,
		params,
	)
	utils.ExportToJson(&index, &path)

	return &index
}

func _prob_index_table(index *map[string][]int,
					   start, end int,
					   seeds *[][]string,
					   d *[][]float64,
					   params *Params)  {

	for i := start; i < end; i++ {
		seed := strings.Join((*seeds)[i], "")
		(*index)[seed] = indexing.ProbFind(&seed, d, params)
	}
}

func prob_extend(q *string, d *[][]float64, index *map[string][]int, params *Params) *[]*Hsp {
	// hsps := [][][]int{}
	hsps := []*Hsp{}

	for q_idx := 0; q_idx < len(*q) - params.W + 1; q_idx++ {
		seed := (*q)[q_idx : q_idx + params.W]

		for _, d_idx := range (*index)[seed] {
			hsp_score := ungapped_extension.ScoreMiddle(q_idx, d_idx, q, d, params)

			hsp_left, score_left := ungapped_extension.Left(q_idx, d_idx, q, d, params)
			hsp_right, score_right := ungapped_extension.Right(q_idx, d_idx, q, d, params)

			if score_left == math.Inf(-1) && score_right == math.Inf(-1) {
				// fmt.Println("both")
				continue
			}
			if score_left > math.Inf(-1) {
				// fmt.Println("left")
				hsp_score += score_left
			}
			if score_right > math.Inf(-1) {
				// fmt.Println("right")
				hsp_score += score_right
			}

			if hsp_score >= params.HspThres {
				// hsps = append(hsps, [][]int{*hsp_left, *hsp_right})
				hsps = append(hsps, &Hsp{(*hsp_left)[0], (*hsp_right)[0], (*hsp_left)[1], (*hsp_right)[1], hsp_score})
			}
		}
	}
	return &hsps
}

func prob_extend_gap(q *string, d *[][]float64, hsps *[]*Hsp, params *Params) *[]*Alignment {

	// substitution_matrix := gapped_extension.SubstitutionMatrix(S)
	S_idx := gapped_extension.SIdx(&params.S)

	// fmt.Println(substitution_matrix)

	alignments := []*Alignment{}
	
	// fmt.Println("num hsps:", len(*hsps))

	// for i, hsp := range *hsps {
	for i := 0; i < len(*hsps); i++ {
		// fmt.Println("hsp", i, ":")

		alignment := gapped_extension.Extend((*hsps)[i], q, d, S_idx, params)
		// fmt.Println(alignment.QAligned)
		// fmt.Println(alignment.DAligned)
		// fmt.Println("score:", alignment.Score)
		alignments = append(alignments, &alignment)
	}
	return &alignments
}