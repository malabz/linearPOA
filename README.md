# Memory-efficient Partial Order Alignment with linearPOA

linearPOA is a library/program written in C++17 for applying Hirschberg algorithm to Partial Order Alignment (POA). It runs on Linux and Windows.

## Test Environment 

linearPOA is a standalone program. In order to test our methods, we need to compile the program first:

```bash
cd src
make -j16 DEBUG=0
```

You will find `linearPOA` program in `src` folder.

## Usage

```bash
Available options:
        --in        FILE      sequence file name (Required)
        --out       FILE      output file name (Required)
        --threads   N         use N threads (N >= 1, default: 1)
        --open      O         gap open penalty (default: 3)
        --ext       E         gap extension penalty (default: 1)
        --match     M         match score (default: 0)
        --mismatch  X         mismatch score (default: 2)
        --nolinear            do not use linear method (default: disabled)
        --genmode   M         generate mode, 1: generate MSA, 2: generate consensus, 3: generate MSA+consensus (default: 1)
        --help                print help message
        --version             show program version
Example:
        ./linearpoa --in seq.fasta --out seq_out.fasta
```


## Test dataset and compiled program

The test dataset and compiled program are stored at [`https://doi.org/10.5281/zenodo.15637837`](https://zenodo.org/uploads/15637837). You can use the data for testing the program.

## Citation

## Contacts

If you find any bug, welcome to contact us on the [issues page](https://github.com/malabz/FORAlign/issues) or [email us](mailto:wym6912@outlook.com).

More tools and infomation can visit our [github](https://github.com/malabz).