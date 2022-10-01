# C++ Bloom Filter implementation.  
  
Has been made for an internship application.

### Requirements:  
- boost library  
(installation example : `sudo apt install libboost-all-dev`)
- c++ compiler
- fasta file

### ARGS:  
- `--help`: produces help message
- `--file`: fasta file to map into the bloom filter
- `--k`: size of kmers <= 31
- `--n`: size of bloom filter in bits <= 2^34 (bits)
- `--nf`: number of applied hash functions <= 64  
- `--r`: number of random `is_present()` tests

### Use examples:  
- ##### Compilation:  
`g++ bf.cpp -lboost_program_options -o bf`
- ###### Execution:
    - `./bf --help`
    - `./bf --file cov.fasta --k 10 --n 1048576 --nf 5 --r 100000`
    - `./bf --file cov.fasta --k 31 --n 17179869184 --nf 64 --r 200` (max values for k, n, and nf)

### Note:  
It is possible to test a specific DNA Sequence in the `main()` by using the `compute_specific_requests(bf, string kmer);` function.  
There is already an example of size k=10 in `bf.cpp` at line `80`.  