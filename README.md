# COMP 561 Project: Probabilistic BLAST

## What is it
An implementation of BLAST adapted for querying a probabilistic database sequence.

## Instructions

### Go version
Parameters such as:
* `w`: word size for seed
* `hit_thres`: threshold for a seed hit
* `delta`: for ungapped extension
* `hsp_thres`: threshold for HSP score
* `e_thres`: threshold for e-value
* `S`: alphabet

are defined and can be modified in the main Go script `prob_blast.go`.  
There's actually only one Go script, as importing from different files in Go is very tedious.  

To run the script, enter the following command in the terminal:
```
go run prob_blast.go
```
