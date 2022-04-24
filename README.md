# Multiple Sequence Alignment - Simulated Annealing (MSASA)

## Introduction
### About this project
This is the source code related to the Ms. Bioinformatics thesis work made by Adrián Diaz under the supervision of Gabriela Minetti as his thesis director at the Universidad Nacional de Quilmes in Buenos Aires, Argentina.

### Abstract
The main goal of this work is to compare the performance of Simulated Annelaing (SA), a metaheuristic method, performing a Multiple Sequence Alignment (MSA) versus different available tools in the bioinformatics world.

#### About the selected tools
After a research of the current state of the art in the bioinformatics scientific field, the selected tools to compare are:

- Clustal Omega (ClustalO)
- KAlign
- MAFFT
- Muscle

##### Clustal Omega
Clustal Omega is a new multiple sequence alignment program that uses seeded guide trees and HMM profile-profile techniques to generate alignments between three or more sequences. For the alignment of two sequences please instead use our pairwise sequence alignment tools.

Publications:

> Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology 7:539 doi:10.1038/msb.2011.75

> Sievers F, Higgins DG (2018) Clustal Omega for making accurate alignments of many protein sequences. Protein Sci 27:135-145

> Sievers F, Barton GJ, Higgins DG (2020) Multiple Sequence Alignment. Bioinformatics 227, pp 227-250, AD Baxevanis, GD Bader, DS Wishart (Eds)

##### KAlign
A fast and accurate multiple sequence alignment algorithm for DNA and protein sequences. This tool uses the Wu-Manber string-matching algorithm, to improve both the accuracy and speed of multiple sequence alignment. It is a fast and robust alignment method, especially well suited for the increasingly important task of aligning large numbers of sequences. Since Kalign is very fast, a user can align sequences, verify the alignment through visual inspection and re-run Kalign (if necessary) with different parameters. This interactive process can provide more insight into the evolutionary relationships of the sequences to be aligned than when a single alignment is used.

Publications:

> Kalign2: high-performance multiple alignment of protein and nucleotide sequences allowing external features. Lassmann T., Frings, O. and Erik L.L. Sonnhammer (2009). Nucleic Acids Research, 37:858-865

> Kalign, Kalignvu and Mumsa: web servers for multiple sequence alignment. Lassmann T. and Erik L.L. Sonnhammer (2006). Nucleic Acids Research, 34:W596-W599

> Kalign - an accurate and fast multiple sequence alignment algorithm. Lassmann T. and Erik L.L. Sonnhammer (2005). BMC Bioinformatics, 6:298

> Automatic assessment of alignment quality. Lassmann T. and Erik L.L. Sonnhammer (2005). Nucleic Acids Research, 33(22):7120-7128

##### MAFFT
MAFFT is a multiple sequence alignment program for unix-like operating systems.  It offers a range of multiple alignment methods, L-INS-i (accurate; for alignment of <∼200 sequences), FFT-NS-2 (fast; for alignment of <∼30,000 sequences).

Publications:
> Katoh, Misawa, Kuma, Miyata 2002 (Nucleic Acids Res. 30:3059-3066). MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. (describes the FFT-NS-1, FFT-NS-2 and FFT-NS-i strategies)

##### Muscle
MUSCLE stands for MUltiple Sequence Comparison by Log- Expectation. MUSCLE is claimed to achieve both better average accuracy and better speed than ClustalW2 or T-Coffee, depending on the chosen options.
Compared to previous versions, Muscle v5 is much more accurate, is often faster, and scales to much larger datasets. At the time of writing (late 2021), Muscle v5 has the highest scores on multiple alignment benchmarks including Balibase, Bralibase, Prefab and Balifam. It can align tens of thousands of sequences with high accuracy on a low-cost commodity computer (say, an 8-core Intel CPU with 32 Gb RAM). On large datasets, Muscle v5 is 20-30% more accurate than MAFFT and Clustal-Omega.

Publications:

> Edgar, RC (2021), MUSCLE v5 enables improved estimates of phylogenetic tree confidence by ensemble bootstrapping, bioRxiv 2021.06.20.449169. https://doi.org/10.1101/2021.06.20.449169.

#### About the selected test scenarios

There are two dimensions of comparison: length of the sequences and the number of sequences to align at the same time. On the first hand, the length of the sequences may be short or long, while on the other hand the number of sequences may be small quantity or large quantity.

For short-length sequences we consider ones with less than 200 residues, and for small-quantity of sequences we consider up to 100 sequences to align.

An additional dimension to take into account would be the known percentage of similarity (20%, 30%, et cetera)

So far, this project has four groups of interest:

1. Short-length, small-quantity sequences
1. Short-length, large-quantity sequences
1. Long-length, small-quantity sequences
1. Long-length, large-quantity sequences

From different sources of benchmarking sequences and alignments, these are the groups based on the two mentioned dimensions:

|             | Sequence quantity | Max. sequence length (residues) | Residue type  |
|-------------|-------------------|---------------------------------|---------------|
| Short-Small | 4                 | 91                              | Protein (AAs) |
| Short-Large | 5                 | 1192                            | Protein (AAs) |
| Long-Small  | 3133              | 79                              | Protein (AAs) |
| Long-Large  | 113               | 1061                            | Protein (AAs) |

## Installation
The dependencies in use are listed inside the `requirements.txt` file, you could install them using a conda environment.

### Clone the repository to your working environment
To get a copy of this source code you could just clone it to your computer:

```console
$ git clone https://github.com/agdiaz/msasa.git
$ cd msasa
```

### Directory structure

Inside the root directory you will find these subdirectories and files.

#### Root folder (`./`)
This `README.md` is located inside the root folder as well as the Makefile to build a compiled version of this code in C language. However, the performance improvement was not the expected taking into account the excessive time required to convert the code in C. Anyway, you can compile the code by `$ make build`.
Requirements are listed inside the `requirements.txt` file.

#### Scripts directory (`/scripts`)
This directory contains bash scripts to run each MSA tool. The main file inside this directory is `./scripts/batch.sh` which has the responsibility of run all the tools (including MSASA) for a target file.
Given the probabilistic nature of all these tools, it is required to run each tool a number of times to get a confident result overall the differences that may occur on individual executions. Therefore, all scripts have a variable `EXECUTIONS_PER_TOOL` to set the number of runs. The minimum recommended number is 30.

#### Source code directory (`/src`)
This directory contains the Python files of the implementation of SA to resolve the MSA problem. The entry file is `./src/msa.py` where the parser of the command line parameters triggers the execution of the algorithm.

#### Test directory (`/test`)
This directory contains several examples in FASTA file (`.fa` or `.fasta`) to test the tool.

### Creation of the conda environment
Once you have installed anaconda or miniconda on your working environment, please follow these steps to create the conda environment.

```console
$ conda create --name msasa-dev python=3 --yes
$ conda activate msasa-dev
$ conda install --file requirements.txt --yes
```

## Usage
This project has been implemented using two approaches: DataFrames from Pandas or arrays from NumPy, both aiming to get a better performance than the original Python data structures.

### Parameters
This piece of software supports the following parameters defined in `src/msa.py`:

- `--input`: the path to the sequence files to align
- `--output`: the path to the output file
- `--comparer`: the method to compare two sequences (`global_ms`, `global_ms_min`, `blosum`, `matching`)
- `--optimization`: whether the algorithm will perform a maximization or a minimization of the global energy (`min`, `max`) (default=`min`)
- `--n-iterations`: number of iterations of the simulated annealing main loop (default=50)
- `--temperature`: the initial temperature of the simulated model (default=10)
- `--engine`: whether `pandas` or `numpy` data structures to be used by the simulated annealing algorithm to store the states (default=`pandas`)
- `--debug`: be extra verbose printing out status changes to the standard output (default=`False`)

Given this software generates plots to visualize the performance of the execution, these parameters allow you to define the path of the .png plot files:

- `--output-best-plot`: path to the best results plot file
- `--output-temp-plot`: path to the temperatures plot file

#### About the sequences comparer parameter
Independently of the comparer function chosen, the process will be the same for any MSA to analyze: There will be a combination of all the sequences and for each pair of sequences the selected method will be run giving a number value. The total energy is the sum of all this pair values. For any pair of sequences, there will be an iteration position by position comparing the two residues.

##### Method `global_ms`
Identical characters are given 5 points, 4 point is deducted for each non-identical character, 3 points are deducted when opening a gap, and 0.1 points are deducted when extending it

##### Method `global_ms_min`
Identical characters are given 5 points, 4 point is deducted for each non-identical character, 3 points are deducted when opening a gap, and 0.1 points are deducted when extending it.

##### Method `blosum`
Using the BLOSUM62 matrix to get the score of the pair of residues.

##### Method `matching`
Identical characters are given -2 points, 0 point is deducted for each non-identical character, 2 points are deducted when opening a gap, and 1 points are deducted when extending it

### Quick start
You could use the rule `test` from the Makefile:

```console
$ make test
```

Replace the environment variables with your own files (the recommended engine is `numpy`)

```console
$ cd msasa
$ python ./src/msa.py --input $INPUT_FASTA --output $OUTPUT_MSA_FILENAME --comparer global_ms_min --n-iterations 50 --output-best-plot $PLOT_BEST_FILENAME --output-temp-plot $PLOT_TEMP_FILENAME --optimization min --engine numpy >> $LOG_FILENAME
```

An example of the output is:

```console
i       current_temp    best_eval       curr_eval       candidate_eval  diff
00000   10.000000       970.000000      970.000000      970.000000      -6.000000
00001   5.000000        970.000000      970.000000      970.000000      +0.000000
00002   3.333333        970.000000      970.000000      970.000000      +0.000000
00003   2.500000        970.000000      970.000000      1047.000000     +77.000000
```

### Code structure

At the beginning, an initial state is generated together its energy. In the scope of the MSA problem, a state is a MSA that is a list of sequences merged with GAPs (represented by the symbol "-") in order to make all the sequences to have the same length. Therefore, the energy of a MSA is a number that represents the quality of the multiple alignment.
The idea of this solution to the problem of generating a MSA is to introduce changes in the state and compare the energies while performing a reduction of the temperature until the number of iterations was reached. As long as the energy value is better than the current one, the state considered the best is replaced with this new one. Otherwise, if this new state is not better than the current best one, there is a sightly probability to replace it only if the Metropolis condition is met. This condition depends on the temperature and provides the algorithm a mechanism to explore nearby states preventing getting stuck in local minimum or maximum regions of the search space of solutions.

#### Workflow

`src/MSA.py` is the entry point of this software, once it parses the parameters, it will instantiate the MSA engine depending on the `engine` parameter (either Pandas or NumPy). The engine and the parsed parameters are used to instantiate an instance of `Runner` and invoke its `start()` method.

The `start` method will execute the SA loop coded inside the engine instance. As soon as the loop is done, the results will be stored into the output file plotting the internal data if the user sent the parameters related to plotting.


#### About the different implementations

##### Pandas
The states of the process are handled with DataFrames.

##### Numpy
This version seems to be faster than the implementation that used DataFrame.

## Next steps

- [] Parallelize the execution iterations
- [] Implement Hidden Markov Model to improve the performance