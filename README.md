# MPHFcan

Minimal perfect hashing of canonical k-mers

## Compilation

```
gcc -O3 -D_FILE_OFFSET_BITS=64 -pthread -mbmi -o canonical fasta.c canonical.c -lm
```

## Usage


```
./canonical <fasta> <k> <t>
```

where 

* `fasta` is a (multiple) fasta file
* `k` is the k-mer length, optional, default 5, maximum 31
* `t` is the number of buckets, optional, default 4

Output: distribution of k-mers to buckets


## 2-bit encoding

To switch to standard 2-bit encoding, (uncomment) the following lines:

```
// process_string(seq,k,threads,t)
   process_string_std(seq,k,threads,t)
```


## License

* Licensed under the [GNU general public license](https://gitlab.ub.uni-bielefeld.de/gi/sans/blob/master/LICENSE).

