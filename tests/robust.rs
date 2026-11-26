//! Regression and robustness checks for the edge options.
//!
//! 1. Defaults are unchanged: the toy example reproduces the outputs of the
//!    released version (tests/golden, written by v0.1.3).
//! 2. --clusters is reproducible: same seed, same bytes.
//! 3. Already-harmonised labs: when one outbreak clone happens to read one
//!    dilution higher in lab B, the classic pair count turns it into a 2x offset;
//!    --support-unit cluster gives one vote per clone and finds no offset.
//! 4. A real offset carried by several clones is still recovered.
//! 5. A lab with too few independent clones is reported as unlinked, not fitted.
use std::collections::HashMap;
use std::fmt::Write as _;
use std::path::PathBuf;
use std::process::Command;

fn tmp(name: &str) -> String {
    let d = std::env::temp_dir().join(format!("phenocal_robust_{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d.join(name).to_str().unwrap().to_string()
}

fn run(args: &[&str]) {
    let st = Command::new(env!("CARGO_BIN_EXE_phenocal"))
        .args(args)
        .output()
        .expect("failed to launch phenocal");
    assert!(
        st.status.success(),
        "phenocal failed: {}",
        String::from_utf8_lossy(&st.stderr)
    );
}

/// cohort -> (fold, anchor column) from an offsets table.
fn folds(prefix: &str) -> HashMap<String, (f64, String)> {
    let t = std::fs::read_to_string(format!("{prefix}.offsets.tsv")).unwrap();
    t.lines()
        .skip(1)
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (
                f[0].to_string(),
                (f[3].parse().unwrap(), f.last().unwrap().to_string()),
            )
        })
        .collect()
}

/// Same cells, numbers equal to 1e-4 (tolerates last-digit libm differences across OS).
fn same_table(a: &str, b: &str) -> bool {
    let (la, lb): (Vec<&str>, Vec<&str>) = (a.lines().collect(), b.lines().collect());
    la.len() == lb.len()
        && la.iter().zip(&lb).all(|(x, y)| {
            let (cx, cy): (Vec<&str>, Vec<&str>) =
                (x.split('\t').collect(), y.split('\t').collect());
            cx.len() == cy.len()
                && cx.iter().zip(&cy).all(|(p, q)| {
                    p == q
                        || matches!((p.parse::<f64>(), q.parse::<f64>()),
                                    (Ok(u), Ok(v)) if (u - v).abs() <= 1e-4 * (1.0 + u.abs()))
                })
        })
}

#[test]
fn defaults_reproduce_released_outputs() {
    let prefix = tmp("toy");
    run(&[
        "--phenotypes",
        "examples/toy_phenotypes.tsv",
        "--pairs",
        "examples/toy_pairs.tsv",
        "--anchor",
        "labA",
        "--out",
        &prefix,
    ]);
    for ext in [
        "offsets",
        "edges",
        "sigma",
        "harmonised",
        "labels",
        "cutoff",
    ] {
        let golden = PathBuf::from(format!("tests/golden/toy.{ext}.tsv"));
        let got = std::fs::read_to_string(format!("{prefix}.{ext}.tsv")).unwrap();
        let want = std::fs::read_to_string(&golden).unwrap();
        assert!(
            same_table(&got, &want),
            "{ext}.tsv differs from {golden:?}\n--- got\n{got}\n--- want\n{want}"
        );
    }
}

/// Synthetic study. Lab A is the anchor. `clones` = (cluster, isolates per lab,
/// log10 shift of lab B in that clone, lab C present). Isolates of one clone are
/// 1 SNP apart across labs.
fn write_study(name: &str, clones: &[(&str, usize, f64, bool)]) -> (String, String, String) {
    let (mut ph, mut pr, mut cl) = (
        String::from("sample\tvalue\tcohort\n"),
        String::from("sample_i\tsample_j\tdistance\n"),
        String::from("sample\tcluster\n"),
    );
    for (ci, &(clone, n, shift_b, with_c)) in clones.iter().enumerate() {
        let base = 0.5 + 0.25 * ci as f64; // mg/L, differs between clones
        let mut ids: Vec<(String, &str)> = Vec::new();
        for k in 0..n {
            for (lab, v) in [("labA", base), ("labB", base * 10f64.powf(shift_b))] {
                ids.push((format!("{clone}_{lab}_{k}"), lab));
                writeln!(ph, "{clone}_{lab}_{k}\t{v}\t{lab}").unwrap();
            }
            if with_c {
                ids.push((format!("{clone}_labC_{k}"), "labC"));
                writeln!(ph, "{clone}_labC_{k}\t{base}\tlabC").unwrap();
            }
        }
        for (i, (si, la)) in ids.iter().enumerate() {
            writeln!(cl, "{si}\t{clone}").unwrap();
            for (sj, lb) in &ids[i + 1..] {
                if la != lb {
                    writeln!(pr, "{si}\t{sj}\t1").unwrap();
                }
            }
        }
    }
    let (p, q, c) = (
        tmp(&format!("{name}.ph.tsv")),
        tmp(&format!("{name}.pr.tsv")),
        tmp(&format!("{name}.cl.tsv")),
    );
    std::fs::write(&p, ph).unwrap();
    std::fs::write(&q, pr).unwrap();
    std::fs::write(&c, cl).unwrap();
    (p, q, c)
}

const ROBUST: [&str; 9] = [
    "--support-unit",
    "cluster",
    "--min-support",
    "3",
    "--max-tau",
    "12",
    "--unsupported",
    "drop",
    "--unlinked",
];

#[test]
fn one_outbreak_does_not_set_an_offset() {
    // labs agree; one big outbreak clone reads one dilution higher in lab B
    let study = [
        ("outbreak", 10, 2f64.log10(), false),
        ("c1", 1, 0.0, false),
        ("c2", 1, 0.0, false),
        ("c3", 1, 0.0, false),
        ("c4", 1, 0.0, false),
    ];
    let (ph, pr, cl) = write_study("outbreak", &study);
    let classic = tmp("outbreak_classic");
    run(&[
        "--phenotypes",
        &ph,
        "--pairs",
        &pr,
        "--anchor",
        "labA",
        "--quick",
        "--out",
        &classic,
    ]);
    let robust = tmp("outbreak_robust");
    let mut a = vec![
        "--phenotypes",
        &ph,
        "--pairs",
        &pr,
        "--anchor",
        "labA",
        "--quick",
        "--clusters",
        &cl,
        "--out",
        &robust,
    ];
    a.extend(ROBUST);
    a.push("zero");
    run(&a);
    let (fc, fr) = (folds(&classic)["labB"].0, folds(&robust)["labB"].0);
    assert!(
        (fc - 2.0).abs() < 0.2,
        "classic should follow the outbreak (got {fc})"
    );
    assert!(
        (fr - 1.0).abs() < 0.1,
        "one vote per clone should find no offset (got {fr})"
    );
}

#[test]
fn real_offset_across_clones_is_recovered_and_unlinked_lab_is_flagged() {
    // lab B reads 4x higher in every clone; lab C shares only one clone
    let s4 = 4f64.log10();
    let study = [
        ("k1", 2, s4, true),
        ("k2", 3, s4, false),
        ("k3", 1, s4, false),
        ("k4", 2, s4, false),
    ];
    let (ph, pr, cl) = write_study("real", &study);
    let robust = tmp("real_robust");
    let mut a = vec![
        "--phenotypes",
        &ph,
        "--pairs",
        &pr,
        "--anchor",
        "labA",
        "--quick",
        "--clusters",
        &cl,
        "--out",
        &robust,
    ];
    a.extend(ROBUST);
    a.push("zero");
    run(&a);
    let f = folds(&robust);
    assert!(
        (f["labB"].0 - 4.0).abs() < 0.3,
        "lab B offset not recovered: {}",
        f["labB"].0
    );
    assert_eq!(
        f["labC"].1, "unlinked",
        "lab C has one clone (< 3): must be unlinked"
    );
    assert!((f["labC"].0 - 1.0).abs() < 1e-9);
    // without --unlinked zero the same data must fail loudly
    let strict = tmp("real_strict");
    let mut b = vec![
        "--phenotypes",
        &ph,
        "--pairs",
        &pr,
        "--anchor",
        "labA",
        "--quick",
        "--clusters",
        &cl,
        "--out",
        &strict,
    ];
    b.extend(&ROBUST[..8]);
    let out = Command::new(env!("CARGO_BIN_EXE_phenocal"))
        .args(&b)
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "a disconnected graph must fail by default"
    );
}

#[test]
fn clusters_bootstrap_is_reproducible() {
    let study = [
        ("k1", 3, 0.3, false),
        ("k2", 2, 0.2, false),
        ("k3", 4, 0.35, false),
        ("k4", 1, 0.25, false),
    ];
    let (ph, pr, cl) = write_study("repro", &study);
    let outs: Vec<String> = (0..2)
        .map(|i| {
            let p = tmp(&format!("repro{i}"));
            run(&[
                "--phenotypes",
                &ph,
                "--pairs",
                &pr,
                "--anchor",
                "labA",
                "--quick",
                "--clusters",
                &cl,
                "--bootstrap",
                "500",
                "--out",
                &p,
            ]);
            std::fs::read_to_string(format!("{p}.offsets.tsv")).unwrap()
        })
        .collect();
    assert_eq!(
        outs[0], outs[1],
        "same seed must give the same cluster-effective intervals"
    );
}
