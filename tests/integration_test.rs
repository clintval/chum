use assert_cmd::Command;
use tempfile::NamedTempFile;

/// 100 bp sequence that uniquely maps to `test-contig:11-110` (1-based closed) in
/// the mini-BLAST DB at `tests/data/blast/`.
const BLAST_TEST_SEQ: &str = "TCTCCCTGTCTAGGGGGGAGTGCACCCTCCTTAGGCAGTGGGGTCTGTGCTGACCGCCTG\
                               CTGACTGCCTTGCAGGTGAAATTGCCCTGTGGTCCTTGGT";

/// Name of the mini-BLAST database bundled under `tests/data/blast/`.
const BLAST_DB: &str = "hs38DH-chr3:129530791-129531030";

fn chum() -> Command {
    assert_cmd::cargo::cargo_bin_cmd!(env!("CARGO_PKG_NAME"))
}

fn score() -> Command {
    let mut cmd = chum();
    cmd.arg("score");
    cmd
}

fn data(name: &str) -> String {
    format!("tests/data/{name}")
}

/// BED input, no optional flags → only positional columns populated.
#[test]
fn test_run_with_baits_only_bed() {
    let out = NamedTempFile::new().unwrap();
    score()
        .args(["-b", &data("baits.bed"), "-o", out.path().to_str().unwrap()])
        .assert()
        .success();

    let actual = std::fs::read_to_string(out.path()).unwrap();
    let expected = std::fs::read_to_string(data("expected_per_bait.tsv")).unwrap();
    pretty_assertions::assert_eq!(
        actual,
        expected,
        "output differs from expected_per_bait.tsv"
    );
}

/// Interval List input should produce the same output as BED.
#[test]
fn test_run_with_interval_list_input() {
    let out = NamedTempFile::new().unwrap();
    score()
        .args([
            "-b",
            &data("baits.interval_list"),
            "-o",
            out.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    let actual = std::fs::read_to_string(out.path()).unwrap();
    let expected = std::fs::read_to_string(data("expected_per_bait.tsv")).unwrap();
    pretty_assertions::assert_eq!(
        actual,
        expected,
        "interval list output differs from expected_per_bait.tsv"
    );
}

/// `--targets` without `--group-output` must fail.
#[test]
fn test_targets_without_group_output_fails() {
    score()
        .args(["-b", &data("baits.bed"), "-t", &data("targets.bed")])
        .assert()
        .failure();
}

/// `--group-output` without `--targets` must fail.
#[test]
fn test_group_output_without_targets_fails() {
    let out = NamedTempFile::new().unwrap();
    let group_out = NamedTempFile::new().unwrap();
    score()
        .args([
            "-b",
            &data("baits.bed"),
            "-o",
            out.path().to_str().unwrap(),
            "--group-output",
            group_out.path().to_str().unwrap(),
        ])
        .assert()
        .failure();
}

/// `--targets` + `--group-output` together must succeed and write both output files.
#[test]
fn test_run_with_targets_and_group_output() {
    let out = NamedTempFile::new().unwrap();
    let group_out = NamedTempFile::new().unwrap();
    score()
        .args([
            "-b",
            &data("baits.bed"),
            "-o",
            out.path().to_str().unwrap(),
            "-t",
            &data("targets.bed"),
            "--group-output",
            group_out.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    // Per-bait file must have data rows.
    let per_bait = std::fs::read_to_string(out.path()).unwrap();
    let lines: Vec<&str> = per_bait.lines().collect();
    assert!(lines.len() >= 2, "expected at least header + 1 data row");

    // Per-group file must contain header + at least one row.
    let per_group = std::fs::read_to_string(group_out.path()).unwrap();
    let group_lines: Vec<&str> = per_group.lines().collect();
    assert!(
        group_lines.len() >= 2,
        "expected at least header + 1 group row"
    );

    // Header must contain num_baits column.
    assert!(
        group_lines[0].contains("num_baits"),
        "group header missing num_baits"
    );
}

/// Writing to stdout (no `-o` flag) must succeed.
#[test]
fn test_stdout_output() {
    let output = score().args(["-b", &data("baits.bed")]).output().unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(
        stdout.starts_with("bait_name\t"),
        "stdout missing TSV header"
    );
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 4, "expected header + 3 bait rows");
}

/// Full analysis integration test: exercises BLAST (required) and RNAfold (optional).
#[test]
fn test_full_analysis_with_blast() {
    if which::which("blastn").is_err() {
        return;
    }

    // FASTA header uses 1-based closed coordinates (UCSC convention), matching
    // how parse_coord_token interprets chrom:start-end tokens.
    let bait_fasta = NamedTempFile::new().unwrap();
    std::fs::write(
        bait_fasta.path(),
        format!(">bait1 test-contig:11-110\n{BLAST_TEST_SEQ}\n"),
    )
    .unwrap();

    let out = NamedTempFile::new().unwrap();
    let mut cmd = score();
    cmd.args([
        "-b",
        bait_fasta.path().to_str().unwrap(),
        "-o",
        out.path().to_str().unwrap(),
        "--blast-db",
        BLAST_DB,
        "--blast-db-path",
        &data("blast"),
    ]);
    if which::which("RNAfold").is_ok() {
        cmd.arg("--oligo-fold");
    }
    cmd.assert().success();

    let content = std::fs::read_to_string(out.path()).unwrap();
    let mut lines = content.lines();
    let headers: Vec<&str> = lines
        .next()
        .expect("missing header row")
        .split('\t')
        .collect();
    let values: Vec<&str> = lines
        .next()
        .expect("missing data row")
        .split('\t')
        .collect();
    assert!(lines.next().is_none(), "expected exactly 1 data row");

    let col = |name: &str| -> &str {
        headers
            .iter()
            .position(|&h| h == name)
            .map(|i| values[i])
            .unwrap_or_else(|| panic!("column '{name}' not found in output"))
    };

    assert_eq!(col("blast_hits"), "1", "expected exactly 1 BLAST hit");
    assert_eq!(col("blast_top_hit_interval"), "test-contig:11-110");
    assert_eq!(col("blast_top_hit_matches"), "true");
    assert_eq!(col("blast_top_hit_overlaps"), "true");

    if which::which("RNAfold").is_ok() {
        let mfe: f64 = col("min_free_energy")
            .parse()
            .expect("min_free_energy should be a float when --oligo-fold is set");
        assert!(
            mfe < 0.0,
            "expected negative MFE for GC-rich bait at 65 °C, got {mfe}"
        );
        assert!(!col("folding_structure").is_empty());
        assert_eq!(col("oligo_fold_param_file"), "dna_mathews2004.par");
    }
}
