//! CLI smoke test: the published toy example must recover labB ~ 2x and
//! labC ~ 5x (relative to labA), and must emit the label and cutoff tables.
use std::collections::HashMap;
use std::process::Command;

#[test]
fn cli_toy_recovers_2x_5x_and_emits_tables() {
    let out = std::env::temp_dir().join("phenocal_cli_smoke");
    let prefix = out.to_str().unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_phenocal"))
        .args([
            "--phenotypes",
            "examples/toy_phenotypes.tsv",
            "--pairs",
            "examples/toy_pairs.tsv",
            "--anchor",
            "labA",
            "--out",
            prefix,
        ])
        .status()
        .expect("failed to launch phenocal binary");
    assert!(status.success(), "phenocal exited with failure");

    let offsets = std::fs::read_to_string(format!("{prefix}.offsets.tsv")).unwrap();
    let mut fold = HashMap::new();
    for line in offsets.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        fold.insert(f[0].to_string(), f[3].parse::<f64>().unwrap());
    }
    assert!(
        (fold["labB"] - 2.0).abs() < 0.3,
        "labB fold {} not ~2x",
        fold["labB"]
    );
    assert!(
        (fold["labC"] - 5.0).abs() < 0.7,
        "labC fold {} not ~5x",
        fold["labC"]
    );

    assert!(std::path::Path::new(&format!("{prefix}.labels.tsv")).exists());
    assert!(std::path::Path::new(&format!("{prefix}.cutoff.tsv")).exists());
}

/// The chewBBACA-style distance matrix of the toy (upper-triangular, genomes
/// named `<sample>.fasta`) must give exactly the outputs of the long pair table.
#[test]
fn cli_chewbbaca_matrix_matches_long_pairs() {
    let run = |pairs: &str, name: &str| -> String {
        let prefix = std::env::temp_dir().join(name);
        let prefix = prefix.to_str().unwrap().to_string();
        let status = Command::new(env!("CARGO_BIN_EXE_phenocal"))
            .args([
                "--phenotypes",
                "examples/toy_phenotypes.tsv",
                "--pairs",
                pairs,
                "--anchor",
                "labA",
                "--min-support",
                "2",
                "--out",
                &prefix,
            ])
            .status()
            .expect("failed to launch phenocal binary");
        assert!(status.success(), "phenocal failed on {pairs}");
        prefix
    };
    let long = run("examples/toy_pairs.tsv", "phenocal_cli_long");
    let matrix = run("examples/toy_distance_matrix.tsv", "phenocal_cli_matrix");
    for ext in [
        "offsets",
        "sigma",
        "edges",
        "harmonised",
        "labels",
        "cutoff",
    ] {
        let a = std::fs::read_to_string(format!("{long}.{ext}.tsv")).unwrap();
        let b = std::fs::read_to_string(format!("{matrix}.{ext}.tsv")).unwrap();
        assert_eq!(a, b, "{ext}.tsv differs between long pairs and matrix");
    }
}
