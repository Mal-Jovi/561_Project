package main

import (
	"os"
	"fmt"
	"math"
	"bufio"
	"strings"
	"strconv"
	"io/ioutil"
)

func main() {
	w := 4
	hit_thres := 0.9
	delta := 10.
	hsp_thres := 0.
	e_thres := 0.
	S := []string{"A", "T", "G", "C"}

	d := get_prob_seq(
		"data/raw/chr22.maf.ancestors.42000000.complete.boreo.fa",
		get_d_conf("data/raw/chr22.maf.ancestors.42000000.complete.boreo.conf"),
		& S,
	)
	q := seq_from_fasta("query.fa")
	
	prob_blast(& q, & d, w, & S, hit_thres, delta, hsp_thres, e_thres)
}

func prob_blast(
		q * string, d *[][] float64, w int, S *[] string,
		hit_thres, delta, hsp_thres, e_thres float64) {

	// alingments := []int{}
	index := prob_index_table(d, w, S, hit_thres)
	// prob_index_table(d, w, S, hit_thres)

	fmt.Println(index)
}

func prob_index_table(
		d *[][] float64, w int, S *[] string,
		hit_thres float64) map[string][]int {

	index := make(map[string][]int)
	
	seeds := gen_seeds(S, w)

	// Concurrently index database using 4 cores
	num_cores := 4
	for i := 0; i < num_cores - 1; i++ {
		go _prob_index_table(
			& index,
			i * (len(seeds) / num_cores), // Start index
			(i + 1) * (len(seeds) / num_cores), // End index
			& seeds, d, w, S, hit_thres)
	}
	_prob_index_table(& index, (num_cores-1)*(len(seeds) / 4), len(seeds), & seeds, d, w, S, hit_thres)

	return index
}

func _prob_index_table(
		index *map[string][]int, start, end int, seeds *[][] string,
		d *[][] float64, w int, S *[]string, hit_thres float64) {

	for i := start; i < end; i++ {
		seed := strings.Join((* seeds)[i], "")
		(*index)[seed] = prob_find(& seed, d, w, S, hit_thres)
	}
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

func get_prob_seq(fasta string, conf [] float64, S *[] string) [][] float64 {
	seq := seq_from_fasta(fasta)
	
	prob_seq := make([][]float64, len(* S))
	for i := range prob_seq {
		prob_seq[i] = make([]float64, len(seq))
	}

	for j := 0; j < len(seq); j++ {
		for i := 0; i < len(* S); i++ {
			if string(seq[j]) == (* S)[i] {
				prob_seq[i][j] = conf[j]
			
			} else {
				prob_seq[i][j] = (1 - conf[j]) / float64((len(* S) - 1))
			}
		}
	}
	return prob_seq
}

func seq_from_fasta(fasta string) string {
	// dat, err := ioutil.ReadFile(fasta)
	file, err := os.Open(fasta)
	check(err)
	defer file.Close()

	br := bufio.NewReader(file)

	br.ReadString('\n')
	seq, err := br.ReadString('\n')
	check(err)

	return strings.TrimSpace(seq)
}

func get_d_conf(d_conf string) [] float64 {
	data, err := ioutil.ReadFile(d_conf)
	check(err)

	str_conf := strings.Split(strings.TrimSpace(string(data)), " ")
	float_conf := make([] float64, len(str_conf))

	for i, conf := range str_conf {
		if f, err := strconv.ParseFloat(conf, 64); err == nil {
			float_conf[i] = f
		}
	}
	return float_conf
}

func check(e error) {
    if e != nil {
        panic(e)
	}
}