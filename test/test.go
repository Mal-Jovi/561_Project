package main

import "fmt"

func main() {
	s := str()
	fmt.Println(* s)
}

func str() * string {
	s := "Hello"
	return & s
}