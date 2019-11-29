package main

import (
	"os"
	"fmt"
	"math"
	"bufio"
	"strings"
	"strconv"
	"io/ioutil"
	"encoding/json"
)

func main() {
	w := 4
	hit_thres := 0.9
	delta := 10.
	// delta := math.Inf(1)
	hsp_thres := 0.
	e_thres := 0.
	S := []string{"A", "T", "G", "C"}

	d := get_prob_seq(
		"data/raw/chr22.maf.ancestors.42000000.complete.boreo.fa",
		get_d_conf("data/raw/chr22.maf.ancestors.42000000.complete.boreo.conf"),
		& S,
	)
	q := seq_from_fasta("query.fa")

	prob_blast(q, d, w, & S, hit_thres, delta, hsp_thres, e_thres)
}

func prob_blast(q * string, d *[][] float64, w int, S *[] string,
		        hit_thres, delta, hsp_thres, e_thres float64) {

	// alingments := []int{}
	// index := prob_index_table(d, w, S, hit_thres)
	index := prob_index_table(d, w, S, hit_thres)	
	hsps := prob_extend(q, d, w, S, index, hit_thres, delta, hsp_thres, e_thres)
	fmt.Println(* hsps)
	prob_extend_gap(hsps)
}

func prob_index_table(d *[][] float64, w int, S *[] string,
		              hit_thres float64) * map[string][]int {
	
	path := fmt.Sprint("data/processed/prob_index_table.w", w, ".hit_thres", hit_thres,".json")

	if _, err := os.Stat(path); !os.IsNotExist(err) {
		return index_from_json(& path)
	}

	index := make(map[string] []int)
	seeds := gen_seeds(S, w)

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
	index_to_json(& index, & path)

	return & index
}

func _prob_index_table(index * map[string][]int, start, end int, seeds *[][] string,
                       d *[][] float64, w int, S *[] string, hit_thres float64)  {

	for i := start; i < end; i++ {
		seed := strings.Join((* seeds)[i], "")
		(* index)[seed] = prob_find(& seed, d, w, S, hit_thres)
	}
}

func prob_extend(q * string, d *[][] float64, w int, S *[] string,
                 index * map[string][]int, hit_thres, delta, hsp_thres, e_thres float64) *[][][] int {

	hsps := [][][]int{}
	for q_idx := 0; q_idx < len(* q) - w + 1; q_idx++ {
		seed := (* q)[q_idx:q_idx+w]

		for _, d_idx := range (* index)[seed] {
			var hsp_score float64

			hsp_left, score_left := prob_left(q_idx, d_idx, q, d, w, S, hit_thres, delta)
			hsp_right, score_right := prob_right(q_idx, d_idx, q, d, w, S, hit_thres, delta)

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

func prob_left(q_idx, d_idx int, q * string, d *[][] float64, w int, S *[] string,
	           hit_thres, delta float64) (* []int, float64) {

	return _prob_extend(q_idx, d_idx, q, d, w, S, hit_thres, delta, -1)
}

func prob_right(q_idx, d_idx int, q * string, d *[][] float64, w int, S *[] string,
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

		score += (* d)[ slice_index(S, string((* q)[q_idx])) ][d_idx] - hit_thres
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

func prob_extend_gap(hsps *[][][] int) {

}

func prob_find(
		seed * string, d *[][] float64, w int,
		S *[] string, hit_thres float64) [] int {

	indices := []int{}
	thres := float64(w) * math.Log(hit_thres)

	for i := 0; i < len((*d)[0]) - w + 1; i++ {
		sum_log_prob := 0.
		for j := 0; j < w; j++ {
			sum_log_prob += math.Log( (* d)[ slice_index(S, string((*seed)[j])) ][i+j])
		}

		if sum_log_prob > thres {
			indices = append(indices, i)
		}
	}
	return indices
}

func slice_index(S *[] string, el string) int {
	for i := 0; i < len(* S); i++ {
		if (* S)[i] == el { return i }
	}
	return -1
}

func gen_seeds(S *[] string, w int) [][] string {
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

func get_prob_seq(fasta string, conf *[] float64, S *[] string) *[][] float64 {
	seq := seq_from_fasta(fasta)
	
	prob_seq := make([][]float64, len(* S))
	for i := range prob_seq {
		prob_seq[i] = make([]float64, len(* seq))
	}

	for j := 0; j < len(* seq); j++ {
		for i := 0; i < len(* S); i++ {
			if string((* seq)[j]) == (* S)[i] {
				prob_seq[i][j] = (* conf)[j]
			
			} else {
				prob_seq[i][j] = (1 - (* conf)[j]) / float64((len(* S) - 1))
			}
		}
	}
	return & prob_seq
}

func seq_from_fasta(fasta string) * string {
	// dat, err := ioutil.ReadFile(fasta)
	fin, err := os.Open(fasta)
	check(err)
	defer fin.Close()

	br := bufio.NewReader(fin)

	br.ReadString('\n')
	seq, err := br.ReadString('\n')
	check(err)

	seq = strings.TrimSpace(seq)
	return & seq
}

func index_from_json(path * string) * map[string][]int {
	index := make(map[string][]int)
	fin, err := os.Open(*path)
	check(err)
	defer fin.Close()

	data, _ := ioutil.ReadAll(fin)
	fin.Close()
	json.Unmarshal([]byte(data), & index)
	return & index
}

func index_to_json(index * map[string][]int, path * string) {
	data, err := json.Marshal(index)
	check(err)

	fout, err := os.Create(* path)
	check(err)
	defer fout.Close()

	fout.Write(data)
	fout.Close()
}

func get_d_conf(d_conf string) *[] float64 {
	data, err := ioutil.ReadFile(d_conf)
	check(err)

	str_conf := strings.Split(strings.TrimSpace(string(data)), " ")
	float_conf := make([] float64, len(str_conf))

	for i, conf := range str_conf {
		if f, err := strconv.ParseFloat(conf, 64); err == nil {
			float_conf[i] = f
		}
	}
	return & float_conf
}

func check(e error) {
    if e != nil {
        panic(e)
	}
}