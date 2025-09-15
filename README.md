# (JA)TI - (Joint Alignment and) Tree Inference

(JA)TI is a command-line tool for maximum likelihood phylogenetic tree reconstruction. Currently it supports tree search with SPR moves on DNA and protein sequences, and allows two gap handling strategies: treating gaps in alignments as missing data or using the Poisson Indel Process (PIP) ([paper]( https://www.pnas.org/doi/10.1073/pnas.1220450110 )) model.

## Features

- **Multiple substitution models**: JC69, K80, HKY85, TN93, GTR for DNA; WAG, HIVB, BLOSUM for proteins;
- **Gap handling**: PIP model or treating gaps as missing data;
- **Character frequency optimisation**: fixed or empirical frequencies (ML estimates to be added);
- **Tree optimisation**: SPR-based topology and branch length optimisation;
- **Parallel computation**: optional multi-threading support for faster tree search;
- **Comprehensive logging**: debug and info level logging with timestamps;
- **Reproducible results**: optional PRNG seed for deterministic runs.

### Optimisation process

#### Starting conditions

- **Input alignment**: JATI requires aligned sequences in FASTA format as input.
- **Starting tree**: You can provide a starting tree in Newick format using the `--tree-file` option. If no tree is provided, JATI will construct a starting tree using the Neighbor-Joining (NJ) algorithm.
- **Model parameters**: You can optionally provide starting values for model parameters (e.g., λ and μ for PIP, substitution rate matrix values, kappa, frequencies) using the `--params` and `--freqs` options. These values will be used as initial guesses and will be optimised during the run.

#### Iterative optimisation

JATI uses an iterative optimisation approach:

1. **Model parameter optimisation**: maximum likelihood estimation of substitution model parameters;
2. **Tree topology optimisation**: SPR (Subtree Pruning and Regrafting) moves to improve the current tree topology;
3. **Branch length optimisation**: optimises all branch lengths given the current topology.

The process continues until convergence (change in log-likelihood < epsilon) or maximum iterations reached.

### Performance

- PIP model likelihood computation is more intensive and will take longer than treating gaps as missing data, as more computation needs to be done per site;
- Tree topology optimisation with SPR moves scales quadratically with the number of taxa ($n$) and linearly with the number of sites ($m$), making the overall complexity $O(n^2m)$.

### Parallelisation

When built with the `parallel` feature, JATI can use multiple CPU cores to speed up tree search by parallelising the evaluation of possible regrafting positions for each pruning location during SPR moves.

#### Performance considerations

- **Benefits**: Speedup on multi-core systems during tree topology optimisation, especially for larger trees;
- **Memory usage**: Parallel execution may increase memory consumption due to thread-local data structures;
- **Thread count**: JATI automatically detects and uses all available CPU cores. You can control the number of threads using the `RAYON_NUM_THREADS` environment variable:

```bash
# Limit to 4 threads
RAYON_NUM_THREADS=4 ./target/release/jati --seq-file data/dna_alignment.fasta --model JC69

# Use all available cores (default behavior)
./target/release/jati --seq-file data/dna_alignment.fasta --model JC69
```

#### Reproducibility with parallel execution

When using the `parallel` feature, results remain deterministic when a PRNG seed is specified, despite the parallel execution. The parallelisation is designed to maintain reproducible results across runs and the final tree and log-likelihood should be indentical to a run without the `parallel` feature given the same seed.

## Installation

###

Pre-compiled binaries will be available shortly.

### Building from source

#### Prerequisites

- **Minimum Supported Rust Version**: 1.84.0 (detected using [`cargo-msrv`]( https://github.com/foresterre/cargo-msrv ));
- Git.

#### Building instructions

1. Clone the repository:
```bash
git clone https://github.com/acg-team/JATI.git
cd JATI
```

2. Build the project:
**Standard build** (single-threaded):
```bash
cargo build --release
```

**Parallel build** (multi-threaded support):
```bash
cargo build --release --features parallel
```

The compiled binary will be available at `target/release/jati`.


#### Build features

- **`parallel`**: Enables multi-threading support for improved performance on multi-core systems. This feature uses Rayon for parallel computation of likelihood calculations and tree search operations.

## Usage

### Basic usage

```bash
./target/release/jati -s <sequence_file> -m <model>
```

### Example commands

```bash
# Basic DNA analysis with JC69 model and gaps as missing data
./target/release/jati \
  --seq-file data/dna_alignment.fasta \
  --model JC69 \
  --gap-handling missing

# Protein analysis with input tree, WAG substitution model and PIP
./target/release/jati \
  --seq-file data/protein_alignment.fasta \
  --tree-file data/protein_starting_tree.newick \
  --model WAG \
  --gap-handling pip \
  --out-folder results/
```

### Parameters

#### Required parameters

| Parameter    | Short | Description                                     |
| ------------ | ----- | ----------------------------------------------- |
| `--seq-file` | `-s`  | Input sequence file in FASTA format             |
| `--model`    | `-m`  | Evolutionary model (see supported models below) |


#### Optional parameters

| Parameter          | Short | Default     | Description                                                                                                   |
| ------------------ | ----- | ----------- | ------------------------------------------------------------------------------------------------------------- |
| `--out-folder`     | `-d`  | `.`         | Output directory for results                                                                                  |
| `--tree-file`      | `-t`  | None        | Input tree file in Newick format (if not provided, NJ tree will be built)                                     |
| `--run-name`       | `-r`  | None        | Custom identifier for the run                                                                                 |
| `--max-iterations` | `-x`  | `5`         | Maximum number of optimisation iterations                                                                     |
| `--params`         | `-p`  | None        | Model-specific parameters (space-separated), when using PIP the first two parameters will be used for λ and μ |
| `--freqs`          | `-f`  | None        | Stationary frequencies: `π_T π_C π_A π_G` (DNA) or custom (protein)                                           |
| `--freq-opt`       | `-o`  | `empirical` | Frequency optimisation: `fixed` or `empirical`                                                                |
| `--gap-handling`   | `-g`  | `pip`       | Gap handling strategy: `pip` or `missing`                                                                     |
| `--epsilon`        | `-e`  | `1e-5`      | Convergence threshold for optimisation                                                                        |
| `--seed`           |       | None        | PRNG seed for reproducible results                                                                            |

### Supported substitution models

#### DNA Models
- **JC69**: Jukes-Cantor model (equal rates, equal frequencies);
- **K80**: Kimura 2-parameter model (equal frequencies):
  - Parameters: `α` (transition/transversion ratio), `π_T π_C π_A π_G` (stationary frequencies);
- **HKY85/HKY**: Hasegawa-Kishino-Yano model:
  - Parameters: `α` (transition/transversion ratio), `π_T π_C π_A π_G` (stationary frequencies);
- **TN93**: Tamura-Nei model:
  - Parameters: `α₁` `α₂` (transition rates), `π_T π_C π_A π_G` (stationary frequencies);
- **GTR**: General Time Reversible model:
  - Parameters: `r_TC` `r_TA` `r_TG` `r_CA` `r_CG` (rate matrix values, `r_AG` fixed to 1.0), `π_T π_C π_A π_G` (stationary frequencies).

#### Protein Models
- **WAG**: Whelan and Goldman model;
- **HIVB**: HIV between-host model;
- **BLOSUM**: BLOSUM62-based model.

### Gap Handling Strategies

- **PIP**: Uses the Poisson Indel Process model to handle gaps as evolutionary events;
- **Missing**: Treats gaps as missing data during likelihood calculation.

### Frequency Optimisation

- **fixed**: Use provided frequencies without optimisation;
- **empirical**: Calculate frequencies from the input sequences.

## Output

JATI creates an output folder with the following structure:

```
<run_id>_out/
├── <run_id>_start_tree.newick  # Starting tree (input or NJ tree)
├── <run_id>_tree.newick        # Final optimised phylogenetic tree
├── <run_id>_logl.out           # Final log-likelihood value
└── <run_id>.log                # Detailed execution log
```

### Output Files

- **`*_tree.newick`**: The final optimised phylogenetic tree in Newick format with branch lengths;
- **`*_start_tree.newick`**: The starting tree used for optimisation;
- **`*_logl.out`**: The final log-likelihood value of the optimised model and tree;
- **`*.log`**: Comprehensive log file with:
  - Run configuration and parameters;
  - Model setup and gap handling strategy;
  - Iteration-by-iteration optimisation progress;
  - Final parameter values and frequencies;
  - Timing and convergence information.

### Run Identification

If no `--run-name` is provided, runs are identified by timestamp. Output folder format:
- With run name: `<run_name>_<timestamp>_out`
- Without run name: `<timestamp>_out`

## Examples

### DNA analysis with JC69 model

```bash
./target/release/jati \
  --seq-file data/wickd3b_7705.processed.fasta \
  --model JC69 \
  --gap-handling pip \
  --max-iterations 10 \
  --out-folder results/
```

### Protein analysis with WAG model

```bash
./target/release/jati \
  --seq-file data/whela7_Gene_0917.fasta \
  --model WAG \
  --gap-handling missing \
  --freq-opt empirical \
  --run-name protein_analysis
```

### DNA analysis with PIP and GTR models and custom parameters

```bash
./target/release/jati \
  --seq-file data/wickd3b_7705.processed.fasta \
  --model GTR \
  --gap-handling pip \
  --params 2.0 1.5 0.4 0.6 0.7 0.6 0.8 \
  --freqs 0.28 0.22 0.25 0.25 \
  --freq-opt fixed \
  --seed 12345 \
  --epsilon 1e-6
```

### HKY model with empirical frequencies

```bash
./target/release/jati \
  --seq-file data/wickd3b_7705.processed.fasta \
  --model HKY \
  --params 2.0 1.5 2.5 \
  --freq-opt empirical \
  --max-iterations 15
```

### Large dataset analysis with parallel execution

```bash
# Build with parallel support
cargo build --release --features parallel

# Run with thread control
RAYON_NUM_THREADS=4 ./target/release/jati \
  --seq-file data/HIV-1_env_DNA_mafft_alignment.fasta \
  --model GTR \
  --gap-handling pip \
  --max-iterations 20 \
  --run-name hiv_analysis
```

### Test alignments

The repository includes example alignments in the `data/` folder for testing purposes.

- **`wickd3b_7705.processed.fasta`**: A medium size DNA alignment from the phylogenetic benchmarks ([Zhou et al. 2018](https://academic.oup.com/mbe/article/35/2/486/4644721#113627412)).

- **`whela7_Gene_0917.fasta`**: Protein alignment from the same phylogenetic benchmarks ([Zhou et al. 2018](https://academic.oup.com/mbe/article/35/2/486/4644721#113627412)).

- **`HIV-1_env_DNA_mafft_alignment.fasta`**: Large DNA alignment of HIV-1 env gene sequences.

## Logging

At the moment JATI provides comprehensive logging with two levels:
- **Console**: Info-level messages showing run progress and major steps;
- **File**: Debug-level messages with detailed parameter information and convergence details.

Log timestamps are formatted as `HH:MM:SS` for clear timing information. The log file contains the complete run configuration at the start for reproducibility.

## PRNG Seed

For reproducible results, you can specify a PRNG seed:
- If not provided: Uses current timestamp as seed (logged for reference);
- If provided: Uses the specified seed for all random operations.

This ensures deterministic results across runs on the same system with the same input and seed.

## Dependencies

- Built with Rust 2021 edition;
- Uses the `phylo` library for phylogenetic functionality;
- Command-line parsing with `clap`;
- Logging with `log` and `ftail`;
- Timestamps with `ntimestamp`;
- Error handling with `anyhow`;
- Optional: `rayon` for parallel computation (when `parallel` feature is enabled).

## Related Projects

- [``phylo: a high-performance Rust library for phylogenetic analysis and multiple sequence alignment under maximum likelihood and parsimony optimality criteria]( https://github.com/acg-team/rust-phylo )

## Support

For questions, bug reports, or feature requests, please go to [JATI discussion page]( https://github.com/acg-team/JATI/discussions ) and/or open an issue on [GitHub]( https://github.com/acg-team/JATI/issues ).

## Licence and Attributions

This project is developed by Jūlija Pečerska ([GitHub]( https://github.com/junniest ), [email]( mailto:julija.pecerska@zhaw.ch )).

### Licence

This project is licensed under either of
- Apache Licence, Version 2.0 ([LICENSE-APACHE]( LICENSE-APACHE ) or [www.apache.org/licenses/LICENSE-2.0]( http://www.apache.org/licenses/LICENSE-2.0 )), or
- MIT Licence ([LICENSE-MIT]( LICENSE-MIT ) or [opensource.org/licenses/MIT]( http://opensource.org/licenses/MIT ))

at your option.