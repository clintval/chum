//! Interactive RNAfold subprocess wrapper.
//!
//! Spawns a single `RNAfold --noPS --temp=<T>` process and keeps it alive across
//! all bait sequences to amortize process startup cost. Each bait sequence is
//! sent to stdin; two lines are read back: the dot-bracket structure and the MFE
//! annotation line `" (-12.34)"`.
//!
//! # Parameter file
//!
//! By default [`RnaFoldProcess::spawn`] passes the vendored `dna_mathews2004.par`
//! parameter file via `--paramFile`, which configures RNAfold for DNA thermodynamics.

use std::io::{BufRead, BufReader, Write};
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};

use anyhow::{Context, Result, bail};
use tempfile::NamedTempFile;

/// The vendored DNA Mathews 1999 parameter file, embedded at compile time.
pub const DNA_MATHEWS1999_PAR: &[u8] = include_bytes!("data/dna_mathews1999.par");

/// The vendored DNA Mathews 2004 parameter file, embedded at compile time.
pub const DNA_MATHEWS2004_PAR: &[u8] = include_bytes!("data/dna_mathews2004.par");

/// The vendored RNA Andronescu 2007 parameter file, embedded at compile time.
pub const RNA_ANDRONESCU2007_PAR: &[u8] = include_bytes!("data/rna_andronescu2007.par");

/// The vendored RNA Langdon 2018 parameter file, embedded at compile time.
pub const RNA_LANGDON2018_PAR: &[u8] = include_bytes!("data/rna_langdon2018.par");

/// The vendored RNA Turner 1999 parameter file, embedded at compile time.
pub const RNA_TURNER1999_PAR: &[u8] = include_bytes!("data/rna_turner1999.par");

/// The vendored RNA Turner 2004 parameter file, embedded at compile time.
pub const RNA_TURNER2004_PAR: &[u8] = include_bytes!("data/rna_turner2004.par");

/// A bundled ViennaRNA thermodynamic parameter file.
///
/// Each variant corresponds to a PAR file vendored into the binary. Use
/// [`ParFile::bytes`] to obtain the raw bytes for passing to
/// [`RnaFoldProcess::spawn`], or [`ParFile::file_name`] to get the canonical
/// file name for display purposes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, clap::ValueEnum)]
pub enum ParFile {
    /// DNA Mathews 1999 parameters.
    DnaMathews1999,
    /// DNA Mathews 2004 parameters.
    #[default]
    DnaMathews2004,
    /// RNA Andronescu 2007 parameters.
    RnaAndronescu2007,
    /// RNA Langdon 2018 parameters.
    RnaLangdon2018,
    /// RNA Turner 1999 parameters.
    RnaTurner1999,
    /// RNA Turner 2004 parameters.
    RnaTurner2004,
    /// No parameter file; ViennaRNA's RNA defaults.
    ViennaDefault,
}

impl ParFile {
    /// Return the embedded bytes for this parameter file, or `None` for [`ParFile::ViennaDefault`].
    pub fn bytes(self) -> Option<&'static [u8]> {
        match self {
            ParFile::DnaMathews1999 => Some(DNA_MATHEWS1999_PAR),
            ParFile::DnaMathews2004 => Some(DNA_MATHEWS2004_PAR),
            ParFile::RnaAndronescu2007 => Some(RNA_ANDRONESCU2007_PAR),
            ParFile::RnaLangdon2018 => Some(RNA_LANGDON2018_PAR),
            ParFile::RnaTurner1999 => Some(RNA_TURNER1999_PAR),
            ParFile::RnaTurner2004 => Some(RNA_TURNER2004_PAR),
            ParFile::ViennaDefault => None,
        }
    }

    /// Return the canonical file name for this parameter file, or `None` for [`ParFile::ViennaDefault`].
    pub fn file_name(self) -> Option<&'static str> {
        match self {
            ParFile::DnaMathews1999 => Some("dna_mathews1999.par"),
            ParFile::DnaMathews2004 => Some("dna_mathews2004.par"),
            ParFile::RnaAndronescu2007 => Some("rna_andronescu2007.par"),
            ParFile::RnaLangdon2018 => Some("rna_langdon2018.par"),
            ParFile::RnaTurner1999 => Some("rna_turner1999.par"),
            ParFile::RnaTurner2004 => Some("rna_turner2004.par"),
            ParFile::ViennaDefault => None,
        }
    }
}

impl std::fmt::Display for ParFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Delegate to the clap-generated possible-value name for consistency.
        write!(
            f,
            "{}",
            <Self as clap::ValueEnum>::to_possible_value(self)
                .map(|v| v.get_name().to_string())
                .unwrap_or_else(|| format!("{:?}", self))
        )
    }
}

/// An active RNAfold subprocess.
pub struct RnaFoldProcess {
    child: Child,
    /// `Option` so `Drop` can take ownership to close stdin before killing the process.
    stdin: Option<ChildStdin>,
    stdout: BufReader<ChildStdout>,
    temp: f64,
    /// Keeps the temporary param file alive for the lifetime of the process.
    _param_file: Option<NamedTempFile>,
}

impl RnaFoldProcess {
    /// Spawn an `RNAfold` process at the given temperature (°C).
    ///
    /// `param_content` — raw bytes of a ViennaRNA parameter file to pass via
    /// `--paramFile`. Use [`DNA_MATHEWS2004_PAR`] for DNA thermodynamics (the
    /// recommended default) or `None` to rely on RNAfold's built-in RNA defaults.
    pub fn spawn(temp: f64, param_content: Option<&[u8]>) -> Result<Self> {
        let mut cmd = Command::new("RNAfold");
        cmd.arg("--noPS").arg(format!("--temp={}", temp));

        // Write param bytes to a temp file and hand the path to RNAfold.
        // The NamedTempFile is kept alive in the struct so it isn't deleted
        // until the process is dropped.
        let param_file: Option<NamedTempFile> = match param_content {
            Some(bytes) => {
                let mut f = NamedTempFile::new()
                    .with_context(|| "Cannot create temp file for RNAfold param")?;
                f.write_all(bytes)
                    .with_context(|| "Cannot write RNAfold param file")?;
                cmd.arg("--paramFile").arg(f.path());
                Some(f)
            }
            None => None,
        };

        let mut child = cmd
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .with_context(|| "Failed to spawn `RNAfold`; ensure it is on $PATH")?;

        let stdin = child.stdin.take().context("Could not take RNAfold stdin")?;
        let stdout = child
            .stdout
            .take()
            .context("Could not take RNAfold stdout")?;

        Ok(RnaFoldProcess {
            child,
            stdin: Some(stdin),
            stdout: BufReader::new(stdout),
            temp,
            _param_file: param_file,
        })
    }

    /// Fold a single nucleotide sequence. Returns the dot-bracket structure and
    /// minimum free energy (kcal/mol).
    ///
    /// Write `sequence\n` then read two lines:
    ///   1. The input sequence echoed back (or the header)
    ///   2. The structure + MFE: `"((...)). (-7.40)"`
    pub fn fold(&mut self, sequence: &str) -> Result<FoldResult> {
        let stdin = self
            .stdin
            .as_mut()
            .context("RNAfold stdin is already closed")?;
        writeln!(stdin, "{}", sequence)?;
        stdin.flush()?;

        let mut echo_line = String::new();
        if self.stdout.read_line(&mut echo_line)? == 0 {
            bail!("RNAfold exited unexpectedly before echoing the sequence");
        }

        let mut result_line = String::new();
        if self.stdout.read_line(&mut result_line)? == 0 {
            bail!("RNAfold exited unexpectedly before writing the structure");
        }

        parse_rnafold_line(&result_line)
    }

    /// Return the temperature the process was started with.
    pub fn temp(&self) -> f64 {
        self.temp
    }
}

impl Drop for RnaFoldProcess {
    fn drop(&mut self) {
        // Close stdin first so RNAfold receives EOF and can exit on its own.
        drop(self.stdin.take());
        // Kill the process as a fallback (idempotent if it already exited).
        let _ = self.child.kill();
        // Reap the process to avoid zombies.
        let _ = self.child.wait();
    }
}

/// The result of folding a single nucleotide sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct FoldResult {
    /// Dot-bracket secondary-structure notation.
    pub structure: String,
    /// Minimum Gibbs free energy (kcal/mol); more negative = more stable folding.
    pub min_free_energy: f64,
}

/// Parse a single RNAfold output line `"structure (mfe)"` without spawning a process.
pub fn parse_rnafold_line(line: &str) -> Result<FoldResult> {
    let line = line.trim();
    let last_space = line
        .rfind(" (")
        .with_context(|| format!("Unexpected RNAfold line: {:?}", line))?;
    let structure = line[..last_space].trim().to_string();
    let mfe_str = line[last_space + 2..].trim_end_matches(')').trim();
    let mfe: f64 = mfe_str
        .parse()
        .with_context(|| format!("Cannot parse MFE from {:?}", line))?;
    Ok(FoldResult {
        structure,
        min_free_energy: mfe,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::skip_if_missing;

    #[test]
    fn test_parse_rnafold_line_unfolded() {
        // Sequence with no structure
        let result = parse_rnafold_line("............ ( 0.00)").unwrap();
        assert_eq!(result.structure, "............");
        assert!((result.min_free_energy - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_parse_rnafold_line_with_structure() {
        let result = parse_rnafold_line(".((....)).. (-7.40)").unwrap();
        assert_eq!(result.structure, ".((....))..");
        assert!((result.min_free_energy - (-7.40)).abs() < 1e-3);
    }

    #[test]
    fn test_parse_rnafold_line_negative_mfe() {
        let result = parse_rnafold_line("(((.....))) (-12.34)").unwrap();
        assert!((result.min_free_energy - (-12.34)).abs() < 1e-3);
    }

    #[test]
    fn test_parse_rnafold_line_invalid_format_returns_error() {
        assert!(parse_rnafold_line("nospace").is_err());
    }

    #[test]
    fn test_parse_rnafold_line_structure_length() {
        let seq = "ACGTTTTCAAAA"; // 12 bases
        let result = parse_rnafold_line(&format!("{} ( 0.00)", ".".repeat(seq.len()))).unwrap();
        assert_eq!(result.structure.len(), seq.len());
        assert!((result.min_free_energy - 0.0).abs() < 1e-9);
    }

    /// Run only when RNAfold is on $PATH.
    #[test]
    fn test_rnafold_process_spawn_and_fold_if_available() {
        skip_if_missing!("RNAfold");
        let mut proc =
            RnaFoldProcess::spawn(65.0, Some(DNA_MATHEWS2004_PAR)).expect("spawn failed");
        let seq = "ACGTTTTCAAAA";
        let result = proc.fold(seq).expect("fold failed");
        assert_eq!(result.structure.len(), seq.len());
        assert!(
            result.min_free_energy <= 1e-3,
            "AT-rich sequence should have near-zero MFE"
        );
        assert_eq!(proc.temp(), 65.0);
        // Structure must contain only valid dot-bracket characters.
        assert!(
            result
                .structure
                .chars()
                .all(|c| matches!(c, '.' | '(' | ')')),
            "unexpected character in structure: {}",
            result.structure
        );
        // Parentheses must be balanced.
        let open = result.structure.chars().filter(|&c| c == '(').count();
        let close = result.structure.chars().filter(|&c| c == ')').count();
        assert_eq!(open, close, "unbalanced parentheses in structure");
    }

    /// Run only when RNAfold is on $PATH. Verifies a GC-rich stem-loop folds with negative MFE.
    #[test]
    fn test_rnafold_hairpin_folds_if_available() {
        skip_if_missing!("RNAfold");
        let mut proc =
            RnaFoldProcess::spawn(65.0, Some(DNA_MATHEWS2004_PAR)).expect("spawn failed");
        // Six-base GC stem flanking a TTTT loop: expected hairpin with negative MFE.
        let seq = "GCGCGCTTTTGCGCGC";
        let result = proc.fold(seq).expect("fold failed");
        assert_eq!(result.structure.len(), seq.len());
        assert!(
            result
                .structure
                .chars()
                .all(|c| matches!(c, '.' | '(' | ')')),
            "unexpected character in structure: {}",
            result.structure
        );
        let open = result.structure.chars().filter(|&c| c == '(').count();
        let close = result.structure.chars().filter(|&c| c == ')').count();
        assert_eq!(open, close, "unbalanced parentheses");
        assert!(
            result.min_free_energy < 0.0,
            "GC-rich hairpin should have negative MFE, got {}",
            result.min_free_energy
        );
        assert!(
            open >= 3,
            "expected at least 3 base pairs in GC hairpin, got {open}"
        );
    }

    #[test]
    fn test_par_file_display_all_variants() {
        // Every variant must produce a non-empty display string and not panic.
        for variant in [
            ParFile::DnaMathews1999,
            ParFile::DnaMathews2004,
            ParFile::RnaAndronescu2007,
            ParFile::RnaLangdon2018,
            ParFile::RnaTurner1999,
            ParFile::RnaTurner2004,
            ParFile::ViennaDefault,
        ] {
            assert!(
                !variant.to_string().is_empty(),
                "Display was empty for {variant:?}"
            );
        }
    }

    #[test]
    fn test_par_file_display_matches_expected_names() {
        assert_eq!(ParFile::DnaMathews1999.to_string(), "dna-mathews1999");
        assert_eq!(ParFile::DnaMathews2004.to_string(), "dna-mathews2004");
        assert_eq!(ParFile::ViennaDefault.to_string(), "vienna-default");
    }

    #[test]
    fn test_par_file_bytes_all_variants() {
        assert!(ParFile::DnaMathews1999.bytes().is_some());
        assert!(ParFile::DnaMathews2004.bytes().is_some());
        assert!(ParFile::RnaAndronescu2007.bytes().is_some());
        assert!(ParFile::RnaLangdon2018.bytes().is_some());
        assert!(ParFile::RnaTurner1999.bytes().is_some());
        assert!(ParFile::RnaTurner2004.bytes().is_some());
        assert!(ParFile::ViennaDefault.bytes().is_none());
    }

    #[test]
    fn test_par_file_file_name_all_variants() {
        assert_eq!(
            ParFile::DnaMathews1999.file_name(),
            Some("dna_mathews1999.par")
        );
        assert_eq!(
            ParFile::DnaMathews2004.file_name(),
            Some("dna_mathews2004.par")
        );
        assert_eq!(
            ParFile::RnaAndronescu2007.file_name(),
            Some("rna_andronescu2007.par")
        );
        assert_eq!(
            ParFile::RnaLangdon2018.file_name(),
            Some("rna_langdon2018.par")
        );
        assert_eq!(
            ParFile::RnaTurner1999.file_name(),
            Some("rna_turner1999.par")
        );
        assert_eq!(
            ParFile::RnaTurner2004.file_name(),
            Some("rna_turner2004.par")
        );
        assert_eq!(ParFile::ViennaDefault.file_name(), None);
    }

    #[test]
    fn test_rnafold_spawn_with_none_param_if_available() {
        skip_if_missing!("RNAfold");
        let mut proc = RnaFoldProcess::spawn(65.0, None).expect("spawn failed");
        let result = proc.fold("ACGTACGT").expect("fold failed");
        assert_eq!(result.structure.len(), 8);
    }

    /// Run only when RNAfold is on $PATH. Verifies the long-lived process handles multiple
    /// sequential sequences correctly (tests stdin/stdout framing across multiple folds).
    #[test]
    fn test_rnafold_multiple_sequences_if_available() {
        skip_if_missing!("RNAfold");
        let mut proc =
            RnaFoldProcess::spawn(65.0, Some(DNA_MATHEWS2004_PAR)).expect("spawn failed");
        let seqs = ["ACGTTTTCAAAA", "GCGCGCTTTTGCGCGC", "TTTTTTTTTTTTTTTT"];
        for seq in seqs {
            let result = proc.fold(seq).expect("fold failed");
            assert_eq!(
                result.structure.len(),
                seq.len(),
                "structure length mismatch for {seq}"
            );
            assert!(
                result
                    .structure
                    .chars()
                    .all(|c| matches!(c, '.' | '(' | ')'))
            );
            let open = result.structure.chars().filter(|&c| c == '(').count();
            let close = result.structure.chars().filter(|&c| c == ')').count();
            assert_eq!(open, close, "unbalanced parens for {seq}");
        }
    }
}
