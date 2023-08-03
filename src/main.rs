#[macro_use]
extern crate assert_float_eq;

use crate::alignment::Alignment;
use crate::cli::Cli;
use crate::parsimony_alignment::pars_align_on_tree;
use crate::parsimony_alignment::parsimony_costs::parsimony_costs_model::DNAParsCosts;
use crate::parsimony_alignment::parsimony_costs::parsimony_costs_model::ProteinParsCosts;
use crate::phylo_info::setup_phylogenetic_info;
use crate::phylo_info::PhyloInfo;
use crate::sequences::get_sequence_type;
use anyhow::Error;
use clap::Parser;
use log::LevelFilter;
use log::{error, info};
use pretty_env_logger::env_logger::Builder;
use std::path::PathBuf;
use std::result::Result::Ok;
use tree::get_percentiles;

mod alignment;
mod cli;
mod io;
mod parsimony_alignment;
mod phylo_info;
mod sequences;
mod substitution_models;
mod tree;

type Result<T> = std::result::Result<T, Error>;

#[allow(non_camel_case_types)]
type f64_h = ordered_float::OrderedFloat<f64>;

fn cmp_f64() -> impl Fn(&f64, &f64) -> std::cmp::Ordering {
    |a, b| a.partial_cmp(b).unwrap()
}

fn indel_map_align_dna(
    info: &PhyloInfo,
    model_name: String,
    model_params: Vec<f64>,
    go: f64,
    ge: f64,
    categories: u32,
) -> Result<(Vec<Alignment>, Vec<f64>)> {
    let scoring = DNAParsCosts::new(
        &model_name,
        &model_params,
        go,
        ge,
        &get_percentiles(&info.tree.get_all_branch_lengths(), categories),
        false,
        false,
    )?;
    Ok(pars_align_on_tree(&Box::new(&scoring), info))
}

fn indel_map_align_protein(
    info: &PhyloInfo,
    model_name: String,
    _: Vec<f64>,
    go: f64,
    ge: f64,
    categories: u32,
) -> Result<(Vec<Alignment>, Vec<f64>)> {
    let scoring = ProteinParsCosts::new(
        &model_name,
        go,
        ge,
        &get_percentiles(&info.tree.get_all_branch_lengths(), categories),
        false,
        false,
    )?;
    Ok(pars_align_on_tree(&Box::new(&scoring), info))
}

fn main() -> Result<()> {
    Builder::new()
        .filter_level(LevelFilter::Info)
        .format_timestamp_secs()
        .format_module_path(false)
        .init();
    info!("JATI run started");
    let cli = Cli::try_parse()?;
    info!("Successfully parsed the command line parameters");
    let info = setup_phylogenetic_info(cli.seq_file, cli.tree_file);
    match info {
        Ok(info) => {
            match cli.command {
                Some(command) => match command {
                    cli::Commands::IndelMAP { go, ge, categories } => {
                        let (alignment, scores) = match get_sequence_type(&info.sequences) {
                            crate::sequences::SequenceType::DNA => {
                                info!("Working on DNA data -- please ensure that data type is inferred correctly.");
                                indel_map_align_dna(
                                    &info,
                                    cli.model,
                                    cli.model_params,
                                    go,
                                    ge,
                                    categories,
                                )?
                            }
                            crate::sequences::SequenceType::Protein => {
                                info!("Working on protein data -- please ensure that data type is inferred correctly.");
                                indel_map_align_protein(
                                    &info,
                                    cli.model,
                                    cli.model_params,
                                    go,
                                    ge,
                                    categories,
                                )?
                            }
                        };
                        info!("Final alignment scores are: \n{:?}", scores);
                        let out_msa_path = match cli.output_msa_file {
                            Some(path) => path,
                            None => {
                                let path = PathBuf::from("msa.fasta");
                                info!(
                                    "No output file name provided, writing MSA to default file {}.",
                                    path.display()
                                );
                                path
                            }
                        };
                        io::write_sequences_to_file(
                            &alignment::compile_alignment_representation(&info, &alignment, None),
                            out_msa_path,
                        )?;
                        info!("IndelMAP alignment done, quitting.");
                    }
                },
                None => {
                    info!("Can only do IndelMAP alignments at the moment, quitting.");
                    return Ok(());
                }
            };
        }
        Err(error) => {
            error!("Error in input data: {}", error);
            return Ok(());
        }
    }
    Ok(())
}
