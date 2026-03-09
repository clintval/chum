//! In silico hybrid selection bait toolkit.
use std::path::PathBuf;
use std::process;

use anyhow::{Error, Result};
use clap::{CommandFactory, FromArgMatches, Parser, Subcommand};
use env_logger::Env;
use log::*;

use chumlib::ScoreRunArgs;
use chumlib::rnafold::ParFile;

use clap::builder::styling::{AnsiColor, Effects, Style, Styles};

pub(crate) const HEADER: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
pub(crate) const USAGE: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
pub(crate) const LITERAL: Style = AnsiColor::Cyan.on_default().effects(Effects::BOLD);
pub(crate) const PLACEHOLDER: Style = AnsiColor::Cyan.on_default();
pub(crate) const ERROR: Style = AnsiColor::Red.on_default().effects(Effects::BOLD);
pub(crate) const VALID: Style = AnsiColor::Cyan.on_default().effects(Effects::BOLD);
pub(crate) const INVALID: Style = AnsiColor::Yellow.on_default().effects(Effects::BOLD);

/// Cargo's color style
/// [source](https://github.com/crate-ci/clap-cargo/blob/master/src/style.rs)
pub(crate) const CARGO_STYLING: Styles = Styles::styled()
    .header(HEADER)
    .usage(USAGE)
    .literal(LITERAL)
    .placeholder(PLACEHOLDER)
    .error(ERROR)
    .valid(VALID)
    .invalid(INVALID);

/// Long help text for the `score` subcommand.
const SCORE_LONG_ABOUT: &str = "\
Evaluate the effectiveness of baits in a hybrid selection panel.

The primary input is a set of bait intervals in BED, Interval List, or FASTA
format. BED and Interval List inputs require a reference FASTA (--ref) to
extract sequences for computing sequence-based metrics. FASTA input embeds the
sequence directly; coordinates are read from `chrom:start-end` tokens in the
header when present. The primary output is a tab-separated metrics file with
one row per bait in the same order as the input.

TARGET EVALUATION
-----------------

When a target interval file is supplied (`--targets`), `chum score` outputs
per-target aggregate metrics including minimum, mean, and maximum of each
numeric bait field across all overlapping (or nearby) baits. Targets can be
optionally padded which helps with linking near-baits to targets on either
side of the original target interval with `--target-padding`. The aggregate
output file is specified with `--per-target`, which must be supplied
whenever `--targets` is set.

When interpreting mean metric values, note that baits with no value for a
field are excluded from the mean, which may skew results. For example, if a
target has three baits with `blast_hits` of `None`, `Some(2)`, and `None`, the
mean `blast_hits` is reported as `2`, not `0.67`.

BLAST ALIGNMENT
---------------

When `blastn` is on the system path and a BLAST database is supplied
(`--blast-db`), each bait sequence is aligned to the database and summary
statistics from the top two BLAST hits are included in the output. BLAST
metrics help assess how specific a bait is to its intended target in the genome.

A high BLAST hit count does not necessarily indicate poor specificity. A user
must examine the percent identity and e-value of the first and second hits to
judge whether off-target baiting is a concern. Baits from FASTA input without
embedded coordinates will have mappability and RepBase metrics computed from
their BLAST top-hit coordinates when available.

BLAST QUERY COMPLEXITY MASKING
------------------------------

BLAST uses the DUST program by default to mask low-complexity regions in query
sequences before seeding alignments, which limits statistically valid but
biologically uninteresting hits. Because `chum score` evaluates baits against
the whole genome including repeat-rich and low-complexity loci, DUST masking is
disabled by default. Enable it with `--blast-dust` when you want to suppress
alignments driven purely by sequence composition.

MAPPABILITY
-----------

When a block-compressed and tabix-indexed bedGraph of mappability scores is
supplied (`--mappability`), per-base mappability statistics are included in
the output. Mappability reflects how uniquely short sub-sequences of the bait
map back to the reference genome and predicts how well sequenced fragments
from the target region can be unambiguously placed. It is recommended to use
Umap bedGraph files derived from multi-read measure and a _k_-mer size of 36.

REPEAT MASKER / REPBASE
-----------------------

When a block-compressed and tabix-indexed BED file of RepBase features is
supplied (`--rep-base`), each bait is intersected with the feature set and
overlapping repeat names are included in the output. Repeat overlaps are
important for predicting bait synthesis efficacy and the likelihood of
alignment artifacts. For example, simple repeats such as `(GGGGGC)n` can
reduce hybrid selection efficiency, while elements like Alu or LINE/L1 can
confound captures due to high endogenous copy number.

OLIGO SECONDARY STRUCTURE
-------------------------

When `--oligo-fold` is set and `RNAfold` (ViennaRNA) is on the system path,
each bait sequence is folded at the temperature specified by
`--oligo-fold-temp` (default 65 °C). The minimum free energy (MFE, kcal/mol)
and dot-bracket secondary structure string are included in the output. A
highly stable secondary structure (very negative MFE) may reduce
hybridization efficiency. Although `RNAFold` was originally designed to fold
longer RNA oligonucleotides, `chum score` is parameterized to use a DNA
configuration by default for ssDNA baits.";

/// Hybrid selection bait toolkit.
#[derive(Debug, Parser)]
#[command(
    author,
    version,
    color = clap::ColorChoice::Always,
    term_width = 80
)]
#[clap(styles = CARGO_STYLING)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Score hybrid selection baits in silico.
    Score(Score),
}

/// Arguments for the `score` subcommand.
#[derive(Debug, Parser)]
#[command(about, long_about = SCORE_LONG_ABOUT, rename_all = "kebab-case")]
struct Score {
    /// Baits in BED, Interval List, or FASTA format.
    #[arg(short = 'b', long)]
    baits: PathBuf,

    /// Output per-bait metrics TSV. If unset, uses stdout.
    #[arg(short = 'o', long)]
    per_bait: Option<PathBuf>,

    /// Indexed FASTA reference. Required with BED or Interval List input.
    #[arg(short = 'r', long)]
    reference: Option<PathBuf>,

    /// Targets in BED or Interval List format.
    #[arg(short = 't', long)]
    targets: Option<PathBuf>,

    /// Output per-target aggregated metrics TSV.
    #[arg(short = 'g', long)]
    per_target: Option<PathBuf>,

    /// Bases to pad targets when matching baits.
    #[arg(long, default_value_t = 0)]
    target_padding: u32,

    /// BLASTn database name.
    #[arg(long)]
    blast_db: Option<String>,

    /// Directory containing BLAST databases. If unset, uses $BLASTDB.
    #[arg(long)]
    blast_db_path: Option<PathBuf>,

    /// CPU threads per blastn subprocess.
    /// With `--threads N > 1`, N parallel subprocesses already run; set this to 1
    /// unless running single threaded, where the full core count may be beneficial.
    #[arg(long, default_value_t = 1)]
    blast_threads: usize,

    /// Enable BLAST DUST query complexity masking.
    #[arg(long, default_value_t = false)]
    blast_dust: bool,

    /// Block-compressed indexed bedGraph of mappability scores.
    #[arg(long)]
    mappability: Option<PathBuf>,

    /// Block-compressed indexed BED of RepBase features.
    #[arg(long)]
    rep_base: Option<PathBuf>,

    /// Enable oligo secondary structure prediction.
    #[arg(long, default_value_t = false)]
    oligo_fold: bool,

    /// Temperature in °C for oligo folding.
    #[arg(long, default_value_t = 65.0)]
    oligo_fold_temp: f64,

    /// ViennaRNA parameters for oligo folding.
    #[arg(long, default_value_t = ParFile::DnaMathews2004)]
    oligo_fold_par: ParFile,

    /// Parallel worker threads.
    #[arg(long, default_value_t = 1)]
    threads: usize,

    /// How many baits to analyze per batch.
    #[arg(long, default_value_t = 50)]
    batch_size: usize,
}

/// Main binary entrypoint.
#[cfg(not(tarpaulin_include))]
fn main() -> Result<(), Error> {
    let env = Env::default().default_filter_or("info");
    env_logger::Builder::from_env(env).init();

    let matches = Cli::command().term_width(80).get_matches();
    let cli = Cli::from_arg_matches(&matches).unwrap_or_else(|e| e.exit());

    match cli.command {
        Commands::Score(score) => {
            match &score.per_bait {
                Some(output) if output.to_str().unwrap_or("") != "-" => {
                    info!("Output file: {output:?}")
                }
                _ => info!("Output stream: STDOUT"),
            }

            let args = ScoreRunArgs {
                baits: score.baits,
                per_bait: score.per_bait,
                reference: score.reference,
                targets: score.targets,
                per_target: score.per_target,
                target_padding: score.target_padding,
                blast_db: score.blast_db,
                blast_db_path: score.blast_db_path,
                blast_dust: score.blast_dust,
                blast_threads: score.blast_threads,
                mappability: score.mappability,
                rep_base: score.rep_base,
                oligo_fold: score.oligo_fold,
                oligo_fold_temp: score.oligo_fold_temp,
                oligo_fold_par: score.oligo_fold_par,
                threads: score.threads,
                batch_size: score.batch_size,
            };

            match chumlib::run_score(args) {
                Ok(()) => process::exit(0),
                Err(e) => {
                    error!("{e:#}");
                    process::exit(1);
                }
            }
        }
    }
}
