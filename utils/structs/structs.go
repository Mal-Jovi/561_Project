package structs

type Params struct {
	W int `json:w` // word size
	HitThres float64 `json:"hit_thres"`
	Delta float64 `json:"delta"`
	HspThres float64 `json:"hsp_thres"`
	S []string `json:"S"`
	EThres float64 `json:"e_thres"` // not used right now
	D string `json:"d"`
	DConf string `json:"d_conf"`
	Q string `json:"q"`
	NumCores int `json:"num_cores"`
}

type Alignment struct {
	QAligned string
	DAligned string
	Score float64
	Accuracy float64
	EVal float64
	QIndices []int
	DIndices []int
	Hsp Hsp
}

type Hsp struct {
	QIdxLeft int
	QIdxRight int
	DIdxLeft int
	DIdxRight int
	Score float64
	EVal float64
}