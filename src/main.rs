use std::fmt::Display;
use std::fs::File;
use std::io::Write;
use std::result::Result::Ok;

use anyhow::Error;
use clap::Parser;
use log::{debug, info};

use phylo::alphabets::{dna_alphabet, protein_alphabet};
use phylo::evolutionary_models::FrequencyOptimisation;
use phylo::io::write_newick_to_file;
use phylo::likelihood::{ModelSearchCost, TreeSearchCost};
use phylo::optimisers::{
    Compatible, ModelOptimiser, MoveOptimiser, SprOptimiser, TopologyOptimiser,
    TopologyOptimiserPredicate,
};
use phylo::phylo_info::PhyloInfoBuilder;
use phylo::pip_model::{PIPCostBuilder, PIPModel};
use phylo::random::{DefaultGenerator, RandomSource};
use phylo::substitution_models::{
    dna_models::*, protein_models::*, SubstModel, SubstitutionCostBuilder,
};
use phylo::tree::Tree;

mod cli;
use crate::cli::{Cli, ConfigBuilder, GapHandling as Gap, SubstModelId as Model};

type Result<T> = std::result::Result<T, Error>;

macro_rules! pip_optimisation {
    ($optimiser:ty, $model:ty, $cfg:expr, $info:expr, $rng:expr) => {
        run_optimisation::<$optimiser>(
            PIPCostBuilder::new(PIPModel::<$model>::new(&$cfg.freqs, &$cfg.params), $info)
                .build()?,
            $cfg.freq_opt,
            $cfg.max_iters,
            $cfg.epsilon,
            $rng,
        )?
    };
}

macro_rules! subst_optimisation {
    ($optimiser:ty, $model:ty, $cfg:expr, $info:expr, $rng:expr) => {
        run_optimisation::<$optimiser>(
            SubstitutionCostBuilder::new(
                SubstModel::<$model>::new(&$cfg.freqs, &$cfg.params),
                $info,
            )
            .build()?,
            $cfg.freq_opt,
            $cfg.max_iters,
            $cfg.epsilon,
            $rng,
        )?
    };
}

fn main() -> Result<()> {
    let cfg = ConfigBuilder::from(Cli::try_parse()?).setup()?;

    info!("JATI run started.");
    info!("{}", cfg);

    let rng = setup_rng(&cfg);

    info!("Running on sequences from {}.", cfg.seq_file.display());

    let alphabet = match cfg.model {
        Model::JC69 | Model::K80 | Model::HKY85 | Model::HKY | Model::TN93 | Model::GTR => {
            info!("Assuming DNA sequences");
            dna_alphabet()
        }
        Model::WAG | Model::HIVB | Model::BLOSUM => {
            info!("Assuming protein sequences");
            protein_alphabet()
        }
    };

    match &cfg.input_tree {
        Some(tree_file) => info!("Using tree from {}.", tree_file.display()),
        None => info!("No initial tree provided, building NJ tree from sequences."),
    }

    let info = PhyloInfoBuilder::new(cfg.seq_file)
        .tree_file(cfg.input_tree)
        .alphabet(Some(alphabet))
        .build_w_rng(&rng)?;

    info!("Putting start tree in {}", cfg.start_tree.display());

    write_newick_to_file(std::slice::from_ref(&info.tree), cfg.start_tree)?;

    info!(
        "Gap handling: {}.",
        match cfg.gap_handling {
            Gap::PIP => "PIP",
            Gap::Missing => "as missing data",
        }
    );
    let (cost, tree) = match (cfg.gap_handling, cfg.model) {
        (Gap::PIP, Model::JC69) => {
            pip_optimisation!(SprOptimiser, JC69, cfg, info, &rng)
        }
        (Gap::PIP, Model::K80) => {
            pip_optimisation!(SprOptimiser, K80, cfg, info, &rng)
        }
        (Gap::PIP, Model::HKY85) | (Gap::PIP, Model::HKY) => {
            pip_optimisation!(SprOptimiser, HKY, cfg, info, &rng)
        }
        (Gap::PIP, Model::TN93) => {
            pip_optimisation!(SprOptimiser, TN93, cfg, info, &rng)
        }
        (Gap::PIP, Model::GTR) => {
            pip_optimisation!(SprOptimiser, GTR, cfg, info, &rng)
        }
        (Gap::PIP, Model::WAG) => {
            pip_optimisation!(SprOptimiser, WAG, cfg, info, &rng)
        }
        (Gap::PIP, Model::HIVB) => {
            pip_optimisation!(SprOptimiser, HIVB, cfg, info, &rng)
        }
        (Gap::PIP, Model::BLOSUM) => {
            pip_optimisation!(SprOptimiser, BLOSUM, cfg, info, &rng)
        }
        (Gap::Missing, Model::JC69) => {
            subst_optimisation!(SprOptimiser, JC69, cfg, info, &rng)
        }
        (Gap::Missing, Model::K80) => {
            subst_optimisation!(SprOptimiser, K80, cfg, info, &rng)
        }
        (Gap::Missing, Model::HKY85) | (Gap::Missing, Model::HKY) => {
            subst_optimisation!(SprOptimiser, HKY, cfg, info, &rng)
        }
        (Gap::Missing, Model::TN93) => {
            subst_optimisation!(SprOptimiser, TN93, cfg, info, &rng)
        }
        (Gap::Missing, Model::GTR) => {
            subst_optimisation!(SprOptimiser, GTR, cfg, info, &rng)
        }
        (Gap::Missing, Model::WAG) => {
            subst_optimisation!(SprOptimiser, WAG, cfg, info, &rng)
        }
        (Gap::Missing, Model::HIVB) => {
            subst_optimisation!(SprOptimiser, HIVB, cfg, info, &rng)
        }
        (Gap::Missing, Model::BLOSUM) => {
            subst_optimisation!(SprOptimiser, BLOSUM, cfg, info, &rng)
        }
    };

    info!("Putting resulting tree in {}", cfg.out_tree.display());
    write_newick_to_file(std::slice::from_ref(&tree), cfg.out_tree)?;

    info!("Final log-likelihood: {cost}");
    let mut out_logl = File::create(cfg.out_logl)?;
    writeln!(out_logl, "{cost}")?;

    Ok(())
}

fn setup_rng(cfg: &cli::Config) -> DefaultGenerator {
    let seed = match cfg.prng_seed {
        None => {
            let seed = cfg.timestamp.as_u64();
            info!("Using current timestamp in milliseconds as the PRNG seed: {seed}");
            seed
        }
        Some(seed) => {
            info!("Using provided PRNG seed: {seed}");
            seed
        }
    };
    DefaultGenerator::new(seed)
}

fn run_optimisation<MO>(
    cost: impl TreeSearchCost + ModelSearchCost + Display + Clone + Send + Compatible<MO>,
    freq_opt: FrequencyOptimisation,
    max_iterations: usize,
    epsilon: f64,
    rng: &impl RandomSource,
) -> Result<(f64, Tree)>
where
    MO: MoveOptimiser + Default,
{
    let mut cost = cost;
    let mut prev_cost = f64::NEG_INFINITY;
    let mut final_cost = TreeSearchCost::cost(&cost);

    let mut iterations = 0;
    while final_cost - prev_cost > epsilon && iterations < max_iterations {
        iterations += 1;
        info!("Iteration: {}", iterations);

        let move_optimiser = MO::default();

        prev_cost = final_cost;
        let model_optimiser = ModelOptimiser::new(cost, freq_opt);
        let o = TopologyOptimiser::new_with_pred(
            model_optimiser.run()?.cost,
            move_optimiser,
            rng,
            TopologyOptimiserPredicate::gt_epsilon(epsilon),
        )
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
