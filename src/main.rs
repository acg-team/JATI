use std::fmt::Display;
use std::fs::File;
use std::io::Write;
use std::result::Result::Ok;

use anyhow::{bail, Error};
use clap::Parser;

use log::{debug, error, info};

use phylo::alphabets::{dna_alphabet, protein_alphabet};
use phylo::evolutionary_models::FrequencyOptimisation;
use phylo::io::write_newick_to_file;
use phylo::likelihood::{ModelSearchCost, TreeSearchCost};
use phylo::optimisers::{
    Compatible, ModelOptimiser, MoveOptimiser, SprOptimiser, TopologyOptimiser,
};
use phylo::phylo_info::PhyloInfoBuilder;
use phylo::pip_model::{PIPCostBuilder, PIPModel};
use phylo::random::DefaultGenerator;
use phylo::substitution_models::{
    dna_models::*, protein_models::*, SubstModel, SubstitutionCostBuilder,
};
use phylo::tree::Tree;

mod cli;
use crate::cli::{Cli, ConfigBuilder, GapHandling};

type Result<T> = std::result::Result<T, Error>;

macro_rules! run_pip_optimisation {
    ($model:ty, $cfg:expr, $info:expr) => {
        run_optimisation::<SprOptimiser>(
            PIPCostBuilder::new(PIPModel::<$model>::new(&$cfg.freqs, &$cfg.params), $info)
                .build()?,
            $cfg.freq_opt,
            $cfg.max_iters,
            $cfg.epsilon,
        )?
    };
}

macro_rules! run_subst_optimisation {
    ($model:ty, $cfg:expr, $info:expr) => {
        run_optimisation::<SprOptimiser>(
            SubstitutionCostBuilder::new(
                SubstModel::<$model>::new(&$cfg.freqs, &$cfg.params),
                $info,
            )
            .build()?,
            $cfg.freq_opt,
            $cfg.max_iters,
            $cfg.epsilon,
        )?
    };
}

fn main() -> Result<()> {
    let cli = match Cli::try_parse() {
        Ok(cli) => {
            info!("Successfully parsed the command line parameters");
            cli
        }
        Err(error) => {
            bail!("Unable to parse command line arguments: \n {}", error)
        }
    };
    let cfg_build: ConfigBuilder = cli.into();
    let cfg = cfg_build.setup()?;

    info!("JATI run started.");
    info!("{}", cfg);

    info!("Running on sequences from {}.", cfg.seq_file.display());
    match &cfg.input_tree {
        Some(tree_file) => info!("Using tree from {}.", tree_file.display()),
        None => info!("No input tree provided, building NJ tree from sequences."),
    }

    let alphabet = match cfg.model.to_uppercase().as_str() {
        "JC69" | "K80" | "HKY85" | "HKY" | "TN93" | "GTR" => {
            info!("Assuming DNA sequences");
            Some(dna_alphabet())
        }
        "WAG" | "HIVB" | "BLOSUM" => {
            info!("Assuming protein sequences");
            Some(protein_alphabet())
        }
        _ => None,
    };

    let info = PhyloInfoBuilder::new(cfg.seq_file)
        .tree_file(cfg.input_tree)
        .alphabet(alphabet)
        .build()?;

    let (cost, tree) = match cfg.gap_handling {
        GapHandling::PIP => {
            info!("Gap handling: PIP.");
            match cfg.model.as_str() {
                "JC69" => run_pip_optimisation!(JC69, cfg, info),
                "K80" => run_pip_optimisation!(K80, cfg, info),
                "HKY85" | "HKY" => run_pip_optimisation!(HKY, cfg, info),
                "TN93" => run_pip_optimisation!(TN93, cfg, info),
                "GTR" => run_pip_optimisation!(GTR, cfg, info),
                "WAG" => run_pip_optimisation!(WAG, cfg, info),
                "HIVB" => run_pip_optimisation!(HIVB, cfg, info),
                "BLOSUM" => run_pip_optimisation!(BLOSUM, cfg, info),
                _ => bail!("Unknown model: {}", cfg.model),
            }
        }
        GapHandling::Missing => {
            info!("Gap handling: as missing data.");
            match cfg.model.as_str() {
                "JC69" => run_subst_optimisation!(JC69, cfg, info),
                "K80" => run_subst_optimisation!(K80, cfg, info),
                "HKY85" | "HKY" => run_subst_optimisation!(HKY, cfg, info),
                "TN93" => run_subst_optimisation!(TN93, cfg, info),
                "GTR" => run_subst_optimisation!(GTR, cfg, info),
                "WAG" => run_subst_optimisation!(WAG, cfg, info),
                "HIVB" => run_subst_optimisation!(HIVB, cfg, info),
                "BLOSUM" => run_subst_optimisation!(BLOSUM, cfg, info),
                _ => bail!("Unknown model: {}", cfg.model),
            }
        }
    };

    info!("Putting resulting tree in {}", cfg.out_tree.display());

    write_newick_to_file(std::slice::from_ref(&tree), cfg.out_tree)?;
    let out_logl = File::create(cfg.out_logl);

    info!("Final log-likelihood: {}", cost);
    match out_logl {
        Ok(mut file) => {
            writeln!(file, "{cost}")?;
        }
        Err(error) => {
            error!("Error creating log file: {}", error);
        }
    }
    Ok(())
}

fn run_optimisation<MO>(
    cost: impl TreeSearchCost + ModelSearchCost + Display + Clone + Send + Compatible<MO>,
    freq_opt: FrequencyOptimisation,
    max_iterations: usize,
    epsilon: f64,
) -> Result<(f64, Tree)>
where
    MO: MoveOptimiser + Default,
{
    let mut cost = cost;
    let mut prev_cost = f64::NEG_INFINITY;
    let mut final_cost = TreeSearchCost::cost(&cost);

    let rng = DefaultGenerator::default();

    let mut iterations = 0;
    while final_cost - prev_cost > epsilon && iterations < max_iterations {
        iterations += 1;
        info!("Iteration: {}", iterations);

        let move_optimiser = MO::default();

        prev_cost = final_cost;
        let model_optimiser = ModelOptimiser::new(cost, freq_opt);
        let o = TopologyOptimiser::new(model_optimiser.run()?.cost, move_optimiser, &rng)
            .run()
            .unwrap();
        final_cost = o.final_cost;
        cost = o.cost;
    }

    info!("Final cost after {} iterations: {}", iterations, final_cost);
    debug!("Final parameters: {:?}", cost.params());
    debug!("Final frequencies: {:?}", cost.freqs());
    debug!("Final tree: {}", cost.tree());
    Ok((final_cost, cost.tree().clone()))
}
