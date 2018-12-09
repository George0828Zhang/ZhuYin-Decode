# ZhuYin-Decode
### Environment
Ubuntu 18.04, 64-bit
### How to compile
To compile on a i686-m64 system with srilm installed at `$(dir)`, run
```sh
make MACHINE_TYPE=i686-m64 SRIPATH=$(dir) all
```
### How to execute
To execute mapping program, run
```
make map
```
To execute the task in the assignment, run
```
make run
```
To execute the task but with order 3 (trigram) and a diffrent model `ngram.lm`, run
```
make run ORD=3 LM=ngram.lm
```
Alternatively, to decode a single file named `file.txt`, using map `m.map`, language model `gram.lm`, and order 3, run
```
./mydisambig -text file.txt -map m.map -lm gram.lm -order 3
```
The result is printed to stdout this way.
### What is done
1. `mapping.py`: Create a Zhu-yin to Character map from the given Big5-ZhuYin.map
2. `result1/`: Sequences decoded using SRILM's disambig are stored here (1~10.txt).
2. `mydisambig.cpp`: Implemented Viterbi algorithm to decode zhuyin-mixed sequences.
3. bonus: In addition to the default order 2 solution, I added support for tri-gram decoding.
