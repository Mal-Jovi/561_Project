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
	config := "config.json"
	if len(os.Args[1:]) == 1 {
		config = os.Args[1]

	} else if len(os.Args[1:]) > 1 {
		fmt.Println("Invalid arguments. Usage: ./prob_blast [config file]")
	}
	params := utils.ParamsFromJson(&config)
	pretty.Println(params)

	d := utils.GetProbSeq(params.D, utils.GetDConf(params.DConf), &params.S)
	q := utils.SeqFromFasta("query.fa")

	alignments := prob_blast(q, d, params)
	pretty.Println(*alignments)
	utils.SaveAlignments(alignments, params)
}

func prob_blast(q * string, d *[][] float64, params *Params) *[]*Alignment {
	index := prob_index_table(d, params)
	hsps := prob_extend(q, d, index, params)
	alignments := prob_extend_gap(q, d, hsps, params)
	
	return alignments
}

func prob_index_table(d *[][]float64, params *Params) *map[string][]int {
	path := fmt.Sprintf("data/processed/prob_index_table.w%d.hit_thres%.2f.json", params.W, params.HitThres)

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
	fmt.Println("Saved index to", path)

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
	hsps := []*Hsp{}

	for q_idx := 0; q_idx < len(*q) - params.W + 1; q_idx++ {
		seed := (*q)[q_idx : q_idx + params.W]

		for _, d_idx := range (*index)[seed] {
			hsp_score := ungapped_extension.ScoreMiddle(q_idx, d_idx, q, d, params)

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

			if hsp_score >= params.HspThres {
				hsps = append(hsps, &Hsp{(*hsp_left)[0], (*hsp_right)[0], (*hsp_left)[1], (*hsp_right)[1], hsp_score})
			}
		}
	}
	return &hsps
}

func prob_extend_gap(q *string, d *[][]float64, hsps *[]*Hsp, params *Params) *[]*Alignment {
	S_idx := gapped_extension.SIdx(&params.S)
	alignments := []*Alignment{}

	for i := 0; i < len(*hsps); i++ {
		alignment := gapped_extension.Extend((*hsps)[i], q, d, S_idx, params)
		alignments = append(alignments, &alignment)
	}
	return &alignments
}