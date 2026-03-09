//! RepBase repeat-feature annotation via tabix-indexed BED overlap.
//!
//! A bgzf-compressed, tabix-indexed BED file of RepBase features is queried
//! per bait. For each bait, overlapping feature names are collected, sorted,
//! and returned as a `Vec<String>` for the caller to apply via
//! [`crate::metrics::apply_repbase`]. No genome-wide data is held in memory.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};

use crate::intervals::Bait;
use noodles::core::Region;
use noodles::tabix;

/// A tabix-indexed bgzf-compressed BED reader for random-access RepBase lookup.
///
/// Only a `PathBuf` is stored; each call to [`overlapping_features`] opens a
/// fresh file handle so the reader is `Send + Sync` and safe for parallel use.
///
/// [`overlapping_features`]: RepBaseReader::overlapping_features
#[derive(Debug)]
pub struct RepBaseReader {
    path: PathBuf,
}

impl RepBaseReader {
    /// Open a tabix-indexed bgzf-compressed BED file.
    ///
    /// Expects a `{path}.tbi` index file alongside the compressed file.
    pub fn new(path: &Path) -> Result<Self> {
        if !path.exists() {
            bail!("RepBase file not found: {}", path.display());
        }
        let tbi = PathBuf::from(format!("{}.tbi", path.display()));
        if !tbi.exists() {
            bail!(
                "Tabix index not found: {}; run `tabix -p bed {}` to create it",
                tbi.display(),
                path.display()
            );
        }
        Ok(RepBaseReader {
            path: path.to_path_buf(),
        })
    }

    /// Return a sorted list of feature names overlapping `bait`.
    pub fn overlapping_features(&self, bait: &Bait) -> Result<Vec<String>> {
        let region_str = format!("{}:{}-{}", bait.chrom, bait.start + 1, bait.end);
        let region: Region = region_str
            .parse()
            .with_context(|| format!("Invalid region string: {region_str}"))?;

        let mut reader = tabix::io::indexed_reader::Builder::default()
            .build_from_path(&self.path)
            .with_context(|| format!("Cannot open tabix file: {}", self.path.display()))?;

        let query = reader
            .query(&region)
            .with_context(|| format!("Tabix query failed for {region_str}"))?;

        let mut names = Vec::new();
        for result in query {
            let record = result?;
            let line = record.as_ref();
            let fields: Vec<&str> = line.splitn(5, '\t').collect();
            if fields.len() < 4 {
                continue;
            }
            names.push(fields[3].trim().to_string());
        }
        names.sort();
        names.dedup();
        Ok(names)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_bait(chrom: &str, start: u64, end: u64) -> Bait {
        Bait::new(chrom, start, end, "b")
    }

    /// Return sorted feature names overlapping `bait` from in-memory records (for unit testing).
    fn features_from_records(bait: &Bait, records: &[(&str, u64, u64, &str)]) -> Vec<String> {
        let mut names = Vec::new();
        for &(chrom, rec_start, rec_end, name) in records {
            if chrom != bait.chrom {
                continue;
            }
            // Tabix-style overlap: any overlap between [rec_start, rec_end) and [bait.start, bait.end)
            if rec_start < bait.end && rec_end > bait.start {
                names.push(name.to_string());
            }
        }
        names.sort();
        names.dedup();
        names
    }

    #[test]
    fn test_no_overlap_returns_empty() {
        // TAR1 at [0,1) is 1bp left of bait [1,6) → no overlap
        let bait = make_bait("chr1", 1, 6);
        let features = features_from_records(&bait, &[("chr1", 0, 1, "TAR1")]);
        assert!(features.is_empty());
    }

    #[test]
    fn test_overlap_returns_name() {
        let bait = make_bait("chr1", 1, 6);
        let features = features_from_records(&bait, &[("chr1", 1, 6, "(TAACCC)n")]);
        assert_eq!(features, vec!["(TAACCC)n"]);
    }

    #[test]
    fn test_repbase_overlap_case() {
        // Bait [1,6), features:
        //   TAR1 [0,1) → no overlap (1 bp to the left)
        //   (TAACCC)n [1,6) → overlap
        //   L1MC5a [4,5) → overlap (within bait)
        //   MIR3 [6,7) → no overlap (1 bp to the right)
        let bait = make_bait("chr1", 1, 6);
        let records = [
            ("chr1", 0, 1, "TAR1"),
            ("chr1", 1, 6, "(TAACCC)n"),
            ("chr1", 4, 5, "L1MC5a"),
            ("chr1", 6, 7, "MIR3"),
        ];
        let features = features_from_records(&bait, &records);
        assert_eq!(features, vec!["(TAACCC)n", "L1MC5a"]);
    }

    #[test]
    fn test_features_sorted() {
        let bait = make_bait("chr1", 0, 100);
        let records = [
            ("chr1", 0, 100, "Zzz"),
            ("chr1", 0, 100, "Aaa"),
            ("chr1", 0, 100, "Mmm"),
        ];
        let features = features_from_records(&bait, &records);
        assert_eq!(features, vec!["Aaa", "Mmm", "Zzz"]);
    }

    #[test]
    fn test_wrong_chrom_returns_empty() {
        let bait = make_bait("chr2", 0, 100);
        let features = features_from_records(&bait, &[("chr1", 0, 100, "Alu")]);
        assert!(features.is_empty());
    }

    #[test]
    fn test_features_deduped() {
        let bait = make_bait("chr1", 0, 100);
        // Two records with the same name → should appear only once.
        let records = [("chr1", 0, 50, "Alu"), ("chr1", 50, 100, "Alu")];
        let features = features_from_records(&bait, &records);
        assert_eq!(features, vec!["Alu"]);
    }

    #[test]
    fn test_repbase_reader_new_file_not_found() {
        let result = RepBaseReader::new(std::path::Path::new("/nonexistent/path.bed.gz"));
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("not found"));
    }

    #[test]
    fn test_repbase_reader_new_tbi_not_found() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "placeholder").unwrap();
        let result = RepBaseReader::new(f.path());
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("Tabix index not found")
        );
    }

    /// Create a bgzf-compressed, tabix-indexed BED4 file in `dir` from `records`.
    ///
    /// Records must be sorted by chrom then start. Returns the `.gz` path, or
    /// `None` if `bgzip` or `tabix` are not on `$PATH`.
    fn make_tabix_bed(
        dir: &tempfile::TempDir,
        records: &[(&str, u64, u64, &str)],
    ) -> Option<std::path::PathBuf> {
        use std::io::Write;
        if which::which("bgzip").is_err() || which::which("tabix").is_err() {
            return None;
        }
        let raw = dir.path().join("features.bed");
        let mut f = std::fs::File::create(&raw).unwrap();
        for &(chrom, start, end, name) in records {
            writeln!(f, "{chrom}\t{start}\t{end}\t{name}").unwrap();
        }
        drop(f);
        // bgzip compresses in-place: features.bed → features.bed.gz
        let status = std::process::Command::new("bgzip")
            .arg(raw.to_str().unwrap())
            .status()
            .ok()?;
        if !status.success() {
            return None;
        }
        let gz = dir.path().join("features.bed.gz");
        // tabix -p bed creates features.bed.gz.tbi
        let status = std::process::Command::new("tabix")
            .args(["-p", "bed", gz.to_str().unwrap()])
            .status()
            .ok()?;
        if !status.success() {
            return None;
        }
        Some(gz)
    }

    #[test]
    fn test_repbase_reader_new_success_if_tools_available() {
        let dir = tempfile::TempDir::new().unwrap();
        let records = [("chr1", 0u64, 100u64, "Alu")];
        let Some(gz) = make_tabix_bed(&dir, &records) else {
            return; // bgzip or tabix not on $PATH
        };
        let reader = RepBaseReader::new(&gz);
        assert!(
            reader.is_ok(),
            "RepBaseReader::new failed: {:?}",
            reader.unwrap_err()
        );
    }

    #[test]
    fn test_overlapping_features_via_tabix_single_hit() {
        let dir = tempfile::TempDir::new().unwrap();
        let records = [("chr1", 0u64, 100u64, "Alu")];
        let Some(gz) = make_tabix_bed(&dir, &records) else {
            return;
        };
        let reader = RepBaseReader::new(&gz).unwrap();
        let bait = make_bait("chr1", 0, 100);
        let features = reader.overlapping_features(&bait).unwrap();
        assert_eq!(features, vec!["Alu"]);
    }

    #[test]
    fn test_overlapping_features_via_tabix_no_overlap() {
        let dir = tempfile::TempDir::new().unwrap();
        let records = [("chr1", 500u64, 600u64, "L1")];
        let Some(gz) = make_tabix_bed(&dir, &records) else {
            return;
        };
        let reader = RepBaseReader::new(&gz).unwrap();
        let bait = make_bait("chr1", 0, 100);
        let features = reader.overlapping_features(&bait).unwrap();
        assert!(features.is_empty());
    }

    #[test]
    fn test_overlapping_features_via_tabix_multiple_sorted_deduped() {
        let dir = tempfile::TempDir::new().unwrap();
        // Two overlapping records: one Alu appearing twice (dedup), one L1.
        let records = [
            ("chr1", 0u64, 50u64, "Alu"),
            ("chr1", 0u64, 100u64, "L1"),
            ("chr1", 50u64, 100u64, "Alu"), // duplicate name
        ];
        let Some(gz) = make_tabix_bed(&dir, &records) else {
            return;
        };
        let reader = RepBaseReader::new(&gz).unwrap();
        let bait = make_bait("chr1", 0, 100);
        let features = reader.overlapping_features(&bait).unwrap();
        // sorted + deduped → ["Alu", "L1"]
        assert_eq!(features, vec!["Alu", "L1"]);
    }
}
