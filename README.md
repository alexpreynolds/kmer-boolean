# kmer-boolean

This utility tests if a specified kmer is or is not in a set of FASTA sequences provided on standard input, for a given *k*, returning the according "true" or "false" result. Alternatively, the binary can instead report kmers that are found (`--present`) or not found (`--absent`), or either case (`--all`).

## Memory usage

Internally, this test keeps an array of bits to minimize the memory overhead of storing per-kmer presence or absence state. This requires at least 2<sup>2k-3</sup> bytes to store the bitarray. Querying 16mers, for example, will require 537 MB of memory.

For the C++ binary, an additional 8 MB buffer is reserved for storing intermediate sequence data that streams in from the input FASTA file. If the `--read-in-all-at-once` option is used, the sequence data is read into memory all at once. It is strongly recommended to instead use the default streaming option to minimize memory usage and runtime overhead required to read all sequences into memory.

## Runtime

To explore runtime characteristics, we used `hg38` assembly data, looking for kmers that are absent in this genome build from <em>k</em>=2 up to 12, to get a general sense of how runtimes trend.

This test is naturally parallelized by chromosome, as any kmer absent in one chromosome is absent over all chromosomes, by definition. We include scripts compatible with a Slurm job scheduler for local testing.

While we include scripts to parallelize the generation of runtime data for all k-by-chromosome combinations, plotting uses the maximum runtime among all chromosomes as the natural bound for querying kmers for the given value of *k*.

![graph](https://user-images.githubusercontent.com/33584/79310219-6ca02580-7eb0-11ea-9e5b-bf443aebf51b.png)

Some work is still left to do to optimize data structures used for bitsetting kmers.

## Notes regarding FASTA input

Note that searches are not performed on "canonical" DNA kmers, but on all unique kmers. For example, a search for the `AG` 2mer will be treated separately from the reverse complement `CT` 2mer.

The input FASTA file may contain one or more records (so-called "multi-FASTA"). Each FASTA record in a multi-FASTA file is scanned separately, but each record contributes to the overall kmer query report. Split a multi-FASTA file if you want to query kmer distributions for individual records.

Kmers containing hard-masked bases (`N`) are ignored. Soft-masked bases (lowercase bases) are included in queries.

An 8 MB buffer is kept of sequence data streaming in from the input FASTA. It is also possible to use the `--read-in-all-at-once` option to read the entire FASTA records into memory. Using this option would not be recommended for genome-scale input files.

## C++

This C++ implementation includes a custom bitset container which can be sized at runtime. The STL `bitset` library can only be sized at compile time and thus is not used here.

### Compilation

Run `make` to build the `kmer-boolean` binary. Run `make clean` to clean up temporary files.

### Usage

#### Options

Run `kmer-boolean --help` to get a listing of options:

```
kmer-boolean
  version: 1.1
  author:  Alex Reynolds

  Usage:

  $ kmer-boolean [arguments] < input

  Test if specified kmer is or is not in a set of
  FASTA sequences provided on standard input, for a
  given k, returning the according 'true' or 'false'
  result.

  General Options:

  --k=n                          K-value for kmer length (integer)
  --query-kmer=s                 Test or query kmer (string)
 [--present | --absent | --all]  If query kmer is omitted, print all
                                 present or absent kmers, or all kmers
  --read-in-all-at-once          Read all sequence data into memory
                                 before processing (not recommended for
                                 very large FASTA inputs; default is to
                                 stream through input in chunks)
                                 
  Process Flags:

  --help                         Show this usage message
  --version                      Show binary version

```

#### Examples

You can return a list of all kmers found in the input via `--present`:

```
$ echo -e ">foo\nCATTCTC\nGGGAC\n>bar\nTTATAT\n>baz\nTTTATTAG\nACCTCT" | ./kmer-boolean --k=2 --present
AC found
AT found
AG found
CA found
CC found
CT found
CG found
TA found
TC found
TT found
GA found
GG found
```

Conversely, you can return a list of kmers not found in the input via `--absent`:

```
$ echo -e ">foo\nCATTCTC\nGGGAC\n>bar\nTTATAT\n>baz\nTTTATTAG\nACCTCT" | ./kmer-boolean --k=2 --absent
AA not found
TG not found
GC not found
GT not found
```

You can also use `--all` to return a list of all kmers, whether found or not found in the input FASTA:

```
$ echo -e ">foo\nCATTCTC\nGGGAC\n>bar\nTTATAT\n>baz\nTTTATTAG\nACCTCT" | ./kmer-boolean --k=3 --all
AAA not found
AAC not found
AAT not found
AAG not found
ACA not found
ACC found
ACT not found
...
GTT not found
GTG not found
GGA found
GGC not found
GGT not found
GGG found
```

If a filter option (`--present`, etc.) is not specified, `kmer-boolean` will return a list of all present kmers.

You may also query for a specific kmer with the `--query-kmer=<kmer>` option:

```
$ echo -e ">foo\nCATTCTC\nGGGAC\n>bar\nTTATAT\n>baz\nTTTATTAG\nACCTCT" | ./kmer-boolean --k=5 --query-kmer="ACGGC"
ACGGC not found
```

If the queried kmer is found, the exit code for the application is zero ("success"). If the kmer is not found, the exit code is non-zero ("failure").

## Python

The `kmer-boolean.py` script uses the Python `bitarray` library to set up the required number of booleans for the specified *k*, and tests the boolean value for the specified kmer, or else returns those kmers which are observed in the input.

Note that searches are performed not on "canonical" kmers, but all unique kmers. In other words, for example, a search for the `AG` 2mer will be treated separately from the reverse complement `CT` 2mer.

### Testing

This has been tested under Python 3.6.3 on OS X. Use of `conda` or `pip` may be required for installation of the `bitstring` library. 
The `makefile` offers a `python_test` target to run expected "true" and "false" tests of searching for a pair of 2mers from some ad-hoc sequences piped in via standard input to the `kmer-boolean.py` script.
