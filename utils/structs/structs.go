package structs

type Alignment struct {
	QAligned string
	DAligned string
	Hsp Hsp
	Score float64
}

type Params struct {
	W int // word size
	HitThres float64
	Delta float64
	HspThres float64
	S []string

	EThres float64
}

type Hsp struct {
	QIdxLeft int
	QIdxRight int
	DIdxLeft int
	DIdxRight int
	Score float64
}