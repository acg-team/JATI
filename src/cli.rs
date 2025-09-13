use std::fmt::{self, Display};
use std::path::PathBuf;

use clap::Parser;
use ftail::Ftail;
use log::LevelFilter;
use ntimestamp::Timestamp;

use phylo::evolutionary_models::FrequencyOptimisation;

use crate::Result;

#[derive(Clone, clap::ValueEnum, Copy, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub(super) enum GapHandling {
    PIP,
    Missing,
}

#[derive(Clone, clap::ValueEnum, Copy, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub(super) enum SubstModelId {
    WAG,
    HIVB,
    BLOSUM,
    JC69,
    K80,
    HKY,
    HKY85,
    TN93,
    GTR,
}

impl Display for SubstModelId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub(super) struct Cli {
    /// Output folder
    #[arg(short = 'd', long, value_name = "OUTPUT_FOLDER", default_value = ".")]
    pub(super) out_folder: PathBuf,

    /// Run identifier for output
    #[arg(short, long, value_name = "RUN_NAME")]
    pub(super) run_name: Option<String>,

    /// Max iterations for the optimisation
    #[arg(short = 'x', long, value_name = "MAX_ITERATIONS", default_value = "5")]
    pub(super) max_iterations: usize,

    /// Sequence file in fasta format
    #[arg(short, long, value_name = "SEQ_FILE")]
    pub(super) seq_file: PathBuf,

    /// Tree file in newick format
    #[arg(short, long, value_name = "TREE_FILE")]
    pub(super) tree_file: Option<PathBuf>,

    /// Sequence evolution model
    #[arg(short, long, value_name = "MODEL", ignore_case = true)]
    pub(super) model: SubstModelId,

    /// Sequence evolution model parameters, e.g. alpha for k80 and
    /// r_tc r_ta r_tg r_ca r_cg r_ag for GTR (in this order)
    #[arg(short = 'p', long, value_name = "MODEL_PARAMS", num_args = 0..)]
    pub(super) params: Vec<f64>,

    /// Sequence evolution model stationary frequencies
    /// pi_t pi_c pi_a pi_g (in this order)
    #[arg(short = 'f', long, value_name = "FREQUENCIES", num_args = 0..)]
    pub(super) freqs: Vec<f64>,

    /// Frequency optimisation method: fixed, empirical, estimated
    #[arg(
        short = 'o',
        long,
        value_name = "FREQ_OPTIMISATION",
        default_value = "empirical"
    )]
    pub(super) freq_opt: FrequencyOptimisation,

    /// Gap handling method: using the PIP model or treating gaps as missing data
    #[arg(
        short,
        long,
        value_name = "GAP_HANDLING",
        ignore_case = true,
        default_value = "pip"
    )]
    pub(super) gap_handling: GapHandling,

    /// Epsilon value for numerical optimisation
    #[arg(short, long, value_name = "EPSILON", default_value = "1e-5")]
    pub(super) epsilon: f64,

    /// PRNG seed that can be fixed for reproducible results
    #[arg(long = "seed", value_name = "PRNG_SEED")]
    pub(super) prng_seed: Option<u64>,
}

pub struct ConfigBuilder {
    pub timestamp: Timestamp,
    pub out_path: PathBuf,
    pub run_name: Option<String>,
    pub max_iters: usize,
    pub seq_file: PathBuf,
    pub input_tree: Option<PathBuf>,
    pub model: SubstModelId,
    pub params: Vec<f64>,
    pub freqs: Vec<f64>,
    pub freq_opt: FrequencyOptimisation,
    pub gap_handling: GapHandling,
    pub epsilon: f64,
    pub prng_seed: Option<u64>,
}

impl From<Cli> for ConfigBuilder {
    fn from(cli: Cli) -> Self {
        ConfigBuilder {
            timestamp: Timestamp::now(),
            out_path: cli.out_folder,
            run_name: cli.run_name,
            max_iters: cli.max_iterations,
            seq_file: cli.seq_file,
            input_tree: cli.tree_file,
            model: cli.model,
            params: cli.params,
            freqs: cli.freqs,
            freq_opt: cli.freq_opt,
            gap_handling: cli.gap_handling,
            epsilon: cli.epsilon,
            prng_seed: cli.prng_seed,
        }
    }
}

pub struct Config {
    pub timestamp: Timestamp,
    pub out_fldr: PathBuf,
    pub out_tree: PathBuf,
    pub start_tree: PathBuf,
    pub out_logl: PathBuf,
    pub run_id: String,
    pub max_iters: usize,
    pub seq_file: PathBuf,
    pub input_tree: Option<PathBuf>,
    pub model: SubstModelId,
    pub params: Vec<f64>,
    pub freqs: Vec<f64>,
    pub freq_opt: FrequencyOptimisation,
    pub gap_handling: GapHandling,
    pub epsilon: f64,
    pub prng_seed: Option<u64>,
}

impl Display for Config {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Run start time: {}", self.timestamp)?;
        writeln!(f, "Run ID: {}", self.run_id)?;

        writeln!(f, "Input sequence file: {}", self.seq_file.display())?;
        match self.input_tree {
            Some(ref tree_file) => writeln!(f, "Input tree file: {}", tree_file.display())?,
            None => writeln!(f, "No input tree file provided.")?,
        }
        writeln!(f, "Output folder: {}", self.out_fldr.display())?;

        let overmodel = match self.gap_handling {
            GapHandling::PIP => "PIP",
            GapHandling::Missing => "Substitution",
        };
        writeln!(f, "Model setup: {overmodel} model with {} Q ", self.model)?;
        writeln!(f, "Model parameters: {:?}", self.params)?;
        writeln!(f, "Model frequencies: {:?}", self.freqs)?;

        writeln!(
            f,
            "Optimisation setup: frequencies: {:#?}, max iterations: {}, epsilon: {}",
            self.freq_opt, self.max_iters, self.epsilon
        )
    }
}

impl ConfigBuilder {
    pub(crate) fn setup(self) -> Result<Config> {
        let run_id = self.run_name.map_or_else(
            || format!("{}", self.timestamp.as_u64()),
            |name| format!("{name}_{}", self.timestamp.as_u64()),
        );

        let out_fldr = self.out_path.join(format!("{run_id}_out"));
        std::fs::create_dir_all(&out_fldr)?;

        Ftail::new()
            .datetime_format("%H:%M:%S")
            .console(LevelFilter::Info)
            .single_file(
                out_fldr.join(format!("{run_id}.log")).to_str().unwrap(),
                true,
                LevelFilter::Debug,
            )
            .init()?;

        let out_tree = out_fldr.join(format!("{run_id}_tree.newick"));
        let start_tree = out_fldr.join(format!("{run_id}_start_tree.newick"));
        let out_logl = out_fldr.join(format!("{run_id}_logl.out"));

        Ok(Config {
            timestamp: self.timestamp,
            out_fldr,
            out_tree,
            start_tree,
            out_logl,
            run_id,
            max_iters: self.max_iters,
            seq_file: self.seq_file,
            input_tree: self.input_tree,
            model: self.model,
            params: self.params,
            freqs: self.freqs,
            freq_opt: self.freq_opt,
            gap_handling: self.gap_handling,
            epsilon: self.epsilon,
            prng_seed: self.prng_seed,
        })
    }
}
