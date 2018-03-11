# kmer-boolean

This utility tests if a specified kmer is or is not in a set of FASTA sequences provided on standard input, for a given *k*, returning the according "true" or "false" result. Alternatively, the binary will return all kmers and whether they are found or not found. Internally, this test keeps an array of bits to minimize the memory overhead of storing per-kmer presence or absence state.

## C++

This C++ implementation includes a custom bitset container which can be sized at runtime. The STL `bitset` library can only be sized at compile time and thus cannot be used here.

Note that searches are performed not on "canonical" DNA kmers, but all unique kmers. In other words, for example, a search for the `AG` 2mer will be treated separately from the reverse complement `CT` 2mer.

### Compilation

Run `make kmer-boolean` to build the `kmer-boolean` binary. Run `make clean` to clean up temporary files.

### Usage

#### Options

Run `kmer-boolean --help` to get a listing of options:

```
kmer-boolean
  version: 1.0
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
