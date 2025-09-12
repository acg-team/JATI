# (JA)TI - (Joint Alignment and) Tree Inference

(JA)TI is a command-line tool for (for now) phylogenetic tree reconstruction. 

that performs joint alignment and tree inference using various evolutionary models. It supports both DNA and protein sequences and provides multiple gap handling strategies including the PIP (Poisson Indel Process) model.

## Features

- **Multiple evolutionary models**: JC69, K80, HKY85, TN93, GTR for DNA; WAG, HIVB, BLOSUM for proteins
- **Gap handling**: PIP model or treating gaps as missing data
- **Flexible frequency optimization**: Fixed, empirical, or estimated frequencies
- **Tree optimization**: Topology and branch length optimization
- **Comprehensive logging**: Debug and info level logging with timestamps

## Installation

### Prerequisites

- Rust (latest stable version)
- Git

### Building from Source

1. Clone the repository:
```bash
git clone <repository-url>
cd JATI
```

2. Build the project:
```bash
cargo build --release
```

The compiled binary will be available at `target/release/jati`.

### Running Tests

```bash
cargo test
```

## Usage

### Basic Usage

```bash
./target/release/jati -s <sequence_file> -m <model> -g <gap_handling>
```

### Example

```bash
./target/release/jati \
  --seq-file data/sequences.fasta \
  --tree-file data/tree.newick \
  --model JC69 \
  --gap-handling PIP \
  --out-folder results/
```

## Parameters

### Required Parameters

| Parameter        | Short | Description                                     |
| ---------------- | ----- | ----------------------------------------------- |
| `--seq-file`     | `-s`  | Input sequence file in FASTA format             |
| `--model`        | `-m`  | Evolutionary model (see supported models below) |
| `--gap-handling` | `-g`  | Gap handling strategy: `PIP` or `Missing`       |

### Optional Parameters

| Parameter          | Short | Default     | Description                                                               |
| ------------------ | ----- | ----------- | ------------------------------------------------------------------------- |
| `--out-folder`     | `-d`  | `.`         | Output directory for results                                              |
| `--tree-file`      | `-t`  | None        | Input tree file in Newick format (if not provided, NJ tree will be built) |
| `--run-name`       | `-r`  | None        | Custom identifier for the run                                             |
| `--max-iterations` | `-x`  | `5`         | Maximum number of optimization iterations                                 |
| `--params`         | `-p`  | None        | Model-specific parameters (space-separated)                               |
| `--freqs`          | `-f`  | None        | Stationary frequencies: π_T π_C π_A π_G                                   |
| `--freq-opt`       | `-o`  | `empirical` | Frequency optimization: `fixed`, `empirical`, or `estimated`              |
| `--epsilon`        | `-e`  | `1e-5`      | Convergence threshold for optimization                                    |

### Supported Models

#### DNA Models
- **JC69**: Jukes-Cantor model (equal rates, equal frequencies)
- **K80**: Kimura 2-parameter model
  - Parameters: `α` (transition/transversion ratio)
- **HKY85/HKY**: Hasegawa-Kishino-Yano model
  - Parameters: `α` (transition/transversion ratio)
- **TN93**: Tamura-Nei model
  - Parameters: `α₁` `α₂` (transition rates)
- **GTR**: General Time Reversible model
  - Parameters: `r_TC` `r_TA` `r_TG` `r_CA` `r_CG` `r_AG` (rate matrix values)

#### Protein Models
- **WAG**: Whelan and Goldman model
- **HIVB**: HIV Between-host model
- **BLOSUM**: BLOSUM62-based model

### Gap Handling Strategies

- **PIP**: Uses the Poisson Indel Process model to handle gaps as evolutionary events
- **Missing**: Treats gaps as missing data during likelihood calculation

### Frequency Optimization

- **fixed**: Use provided frequencies without optimization
- **empirical**: Calculate frequencies from the input sequences

## Output

JATI creates an output folder with the following structure:

```
<run_id>_out/
├── <run_id>_tree.newick      # Optimized phylogenetic tree
├── <run_id>_logl.out         # Final log-likelihood value
└── <run_id>.log              # Detailed execution log
```

### Output Files

- **`*_tree.newick`**: The optimized phylogenetic tree in Newick format
- **`*_logl.out`**: Contains the final log-likelihood value of the optimized model
- **`*.log`**: Comprehensive log file with:
  - Run configuration and parameters
  - Iteration-by-iteration optimization progress
  - Final parameter values and frequencies
  - Timing information

### Run Identification

If no `--run-name` is provided, runs are identified by timestamp. Output folder format:
- With run name: `<run_name>_<timestamp>_out`
- Without run name: `<timestamp>_out`

## Examples

### DNA Analysis with JC69 Model

```bash
./target/release/jati \
  --seq-file data/dna_sequences.fasta \
  --model JC69 \
  --gap-handling PIP \
  --max-iterations 10 \
  --out-folder results/
```

### Protein Analysis with WAG Model

```bash
./target/release/jati \
  --seq-file data/protein_sequences.fasta \
  --model WAG \
  --gap-handling Missing \
  --tree-file data/starting_tree.newick \
  --freq-opt estimated \
  --run-name protein_analysis
```

### GTR Model with Custom Parameters

```bash
./target/release/jati \
  --seq-file data/sequences.fasta \
  --model GTR \
  --gap-handling PIP \
  --params 1.0 2.0 1.0 1.0 2.0 1.0 \
  --freqs 0.25 0.25 0.25 0.25 \
  --freq-opt fixed
```

## Logging

JATI provides comprehensive logging with two levels:
- **Console**: Info-level messages showing progress
- **File**: Debug-level messages with detailed parameter information

Log timestamps are formatted as `HH:MM:SS.mmm` for precise timing information.

## License

This project is developed by Jūlija Pečerska (julija.pecerska@zhaw.ch).

## Dependencies

- Built with Rust 2021 edition
- Uses the `phylo` library for phylogenetic computations
- Command-line parsing with `clap`
- Logging with `log` and `ftail`
