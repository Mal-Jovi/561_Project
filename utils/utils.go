package utils

import (
	"os"
	"fmt"
	"math"
	"bufio"
	"strings"
	"strconv"
	"io/ioutil"
	"path/filepath"
	"encoding/json"
	. "github.com/Mal-Jovi/561_Project/utils/structs"
)

func SumProbSeq(start int, seg_aligned *string, prob_seq *[][]float64, S_idx *map[string]int, hit_thres float64) float64 {
	sum := 0.
	i := start
	
	for j := 0; j < len(*seg_aligned); j++ {
		char := string((*seg_aligned)[j])
		if char == "-" {
			sum -= 1
			continue
		}
		conf := (*prob_seq)[ (*S_idx)[char] ][i]
		if conf >= hit_thres {
			sum += conf
		} else {
			sum -= hit_thres - conf
		}
		i++
	}
	return sum
}

func SaveAlignments(alignments *[]*Alignment, params *Params) {
	path := fmt.Sprintf("output/alignments.q.%s.w.%d.hit_thres.%.2f.delta.%.2f.e_thres.%.2f.json",
		params.Q, params.W, params.HitThres, params.Delta, params.EThres)
	
	ExportToJson(alignments, &path)
	fmt.Println("Saved alignments to", path)
}

func ExpandMat(mat *[][]float64, expand_by int) *[][]float64 {
	ExpandMatBottom(mat, expand_by)
	ExpandMatRight(mat, expand_by)
	return mat
}

func ExpandMatBottom(mat *[][]float64, expand_by int) *[][]float64 {
	*mat = append(*mat, *Mat(expand_by, len((*mat)[0]))...)
	return mat
}

func ExpandMatRight(mat *[][]float64, expand_by int) *[][]float64 {
	for i := 0; i < len(*mat); i++ {
		(*mat)[i] = append((*mat)[i], make([]float64, expand_by)...)
	}
	return mat
}

func ExpandMatInt(mat *[][]int, expand_by int) *[][]int {
	ExpandMatIntBottom(mat, expand_by)
	ExpandMatIntRight(mat, expand_by)
	return mat
}

func ExpandMatIntBottom(mat *[][]int, expand_by int) *[][]int {
	*mat = append(*mat, *MatInt(expand_by, len((*mat)[0]))...)
	return mat
}

func ExpandMatIntRight(mat *[][]int, expand_by int) *[][]int {
	for i := 0; i < len(*mat); i++ {
		(*mat)[i] = append((*mat)[i], make([]int, expand_by)...)
	}
	return mat
}

func SliceIndex(S *[]string, el string) int {
	for i := 0; i < len(*S); i++ {
		if (*S)[i] == el {
			return i
		}
	}
	return -1
}

func PrettyProbSeg(prob_seq *[][]float64, start, end int, S *[]string, hit_thres float64) *string {
	s := ""

	for j := start; j < end; j++ {
		max_val := math.Inf(-1)
		i_max := -1

		for i := 0; i < len(*prob_seq); i++ {
			if (*prob_seq)[i][j] >= max_val {
				max_val = (*prob_seq)[i][j]
				i_max = i
			}
		}

		if (*prob_seq)[i_max][j] >= hit_thres {
			s += (*S)[i_max]
		} else {
			s += "/"
		}
	}
	return &s
}

func GetProbSeq(fasta string, conf *[]float64, S *[]string) *[][]float64 {
	seq := SeqFromFasta(fasta)
	
	// prob_seq := make([][]float64, len(* S))
	// for i := range prob_seq {
	// 	prob_seq[i] = make([]float64, len(* seq))
	// }
	prob_seq := Mat(len(*S), len(*seq))

	for j := 0; j < len(*seq); j++ {
		for i := 0; i < len(*S); i++ {
			if string((*seq)[j]) == (*S)[i] {
				(*prob_seq)[i][j] = (*conf)[j]
			
			} else {
				(*prob_seq)[i][j] = (1 - (*conf)[j]) / float64(len(*S) - 1)
			}
		}
	}
	return prob_seq
}

func SeqFromFasta(fasta string) * string {
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

func GetDConf(d_conf string) *[] float64 {
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

func Mat(num_rows, num_cols int) *[][]float64 {
	// Make num_rows x num_cols float64 array
	mat := make([][]float64, num_rows)
	
	for i := range mat {
		mat[i] = make([]float64, num_cols)
	}
	// mat := [num_rows][num_cols]float64{}
	return &mat
}

func MatInt(num_rows, num_cols int) *[][]int {
	mat := make([][]int, num_rows)
	
	for i := range mat {
		mat[i] = make([]int, num_cols)
	}
	return &mat
}

func ParamsFromJson(path *string) *Params {
	var params Params
	load_json(path, &params)
	return &params
}

func IndexFromJson(path *string) *map[string][]int {
	index := make(map[string][]int)
	load_json(path, &index)
	return &index
}

func load_json(path *string, target interface{}) {
	fin, err := os.Open(*path)
	check(err)
	defer fin.Close()

	data, _ := ioutil.ReadAll(fin)
	fin.Close()
	err = json.Unmarshal([]byte(data), target)
	check(err)
}

// func IndexToJson(index *map[string][]int, path *string) {
func ExportToJson(data interface{}, path *string) {
	dat, err := json.MarshalIndent(data, "", "    ")
	check(err)

	os.MkdirAll(filepath.Dir(*path), 0700)
	fout, err := os.Create(*path)
	check(err)
	defer fout.Close()

	fout.Write(dat)
	fout.Close()
}

func Max(nums ...float64) (float64, int) {
	// return math.Max(math.Max(a, b), c)
	max := math.Inf(-1)
	argmax := -1

	for i, num := range nums {
		if num > max {
			max = num
			argmax = i
		}
	}
	return max, argmax
}

func check(e error) {
    if e != nil {
        panic(e)
	}
}