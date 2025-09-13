# (JA)TI - (Joint Alignment and) Tree Inference

(JA)TI is a command-line tool for maximum likelihood phylogenetic tree reconstruction. Currently it supports tree search with SPR moves on DNA and protein sequences, and allows two gap handling strategies: treating gaps in alignments as missing data or using the Poisson Indel Process (PIP) ([paper]( https://www.pnas.org/doi/10.1073/pnas.1220450110 )) model.

## Features

- **Multiple substitution models**: JC69, K80, HKY85, TN93, GTR for DNA; WAG, HIVB, BLOSUM for proteins;
- **Gap handling**: PIP model or treating gaps as missing data;
- **Character frequency optimisation**: fixed or empirical frequencies (ML estimates to be added);
- **Tree optimisation**: SPR-based topology and branch length optimisation;
- **Comprehensive logging**: debug and info level logging with timestamps;
- **Reproducible results**: optional PRNG seed for deterministic runs.

## Optimisation process

JATI uses an iterative optimisation approach:

1. **Model parameter optimisation**: maximum likelihood estimation of substitution model parameters;
2. **Tree topology optimisation**: SPR (Subtree Pruning and Regrafting) moves to improve the current tree topology;
3. **Branch length optimisation**: optimises all branch lengths given the current topology.

The process continues until convergence (change in log-likelihood < epsilon) or maximum iterations reached.


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
```bash
cargo build --release
```

The compiled binary will be available at `target/release/jati`.

## Usage

### Basic Usage

```bash
./target/release/jati -s <sequence_file> -m <model>
```

### Example Commands

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

## Parameters

### Required parameters

| Parameter    | Short | Description                                     |
| ------------ | ----- | ----------------------------------------------- |
| `--seq-file` | `-s`  | Input sequence file in FASTA format             |
| `--model`    | `-m`  | Evolutionary model (see supported models below) |


### Optional parameters

| Parameter          | Short | Default     | Description                                                               |
| ------------------ | ----- | ----------- | ------------------------------------------------------------------------- |
| `--out-folder`     | `-d`  | `.`         | Output directory for results                                              |
| `--tree-file`      | `-t`  | None        | Input tree file in Newick format (if not provided, NJ tree will be built) |
| `--run-name`       | `-r`  | None        | Custom identifier for the run                                             |
| `--max-iterations` | `-x`  | `5`         | Maximum number of optimisation iterations                                 |
| `--params`         | `-p`  | None        | Model-specific parameters (space-separated)                               |
| `--freqs`          | `-f`  | None        | Stationary frequencies: `π_T π_C π_A π_G` (DNA) or custom (protein)       |
| `--freq-opt`       | `-o`  | `empirical` | Frequency optimisation: `fixed`, `empirical`, or `estimated`              |
| `--gap-handling`   | `-g`  | `pip`       | Gap handling strategy: `pip` or `missing`                                 |
| `--epsilon`        | `-e`  | `1e-5`      | Convergence threshold for optimisation                                    |
| `--seed`           |       | None        | PRNG seed for reproducible results                                        |

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

### DNA analysis with JC69

```bash
./target/release/jati \
  --seq-file data/dna_alignment.fasta \
  --model JC69 \
  --gap-handling pip \
  --max-iterations 10 \
  --out-folder results/
```

### Protein analysis with WAG

```bash
./target/release/jati \
  --seq-file data/protein_alignment.fasta \
  --model WAG \
  --gap-handling missing \
  --tree-file data/protein_starting_tree.newick \
  --freq-opt estimated \
  --run-name protein_analysis
```

### DNA analysis with PIP and GTR with custom starting parameters and reproducible results

```bash
./target/release/jati \
  --seq-file data/dna_aligment.fasta \
  --model GTR \
  --gap-handling pip \
  --params 2.0 1.5 0.7 0.6 0.8 1.2 3.0 \
  --freqs 0.28 0.22 0.25 0.25 \
  --freq-opt fixed \
  --seed 12345 \
  --epsilon 1e-6
```

### HKY Model with empirical frequencies

```bash
./target/release/jati \
  --seq-file data/dna_alignment.fasta \
  --model HKY \
  --params 0.1 0.4 2.5 \
  --freq-opt empirical \
  --max-iterations 15
```

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

## Performance Notes

- PIP model computation is more intensive than treating gaps as missing data;
- Tree topology optimisation scales with the number of taxa.

## Dependencies

- Built with Rust 2021 edition;
- Uses the `phylo` library for phylogenetic functionality;
- Command-line parsing with `clap`;
- Logging with `log` and `ftail`;
- Timestamps with `ntimestamp`;
- Error handling with `anyhow`.

## License

This project is developed by Jūlija Pečerska (julija.pecerska@zhaw.ch).