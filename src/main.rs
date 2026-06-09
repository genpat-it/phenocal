//! phenocal — genome-anchored calibration of a continuous phenotype across
//! cohorts/batches, with uncertainty-propagated probabilistic labels.
//!
//! Method (Algorithm 1 of de Ruvo, Castelli, Di Pasquale & Radomski): near-clonal cross-cohort pairs are
//! genomic natural controls; a shared genome cancels the biological component,
//! so the phenotype difference between matched isolates estimates the
//! cohort/protocol offset difference. Offsets are fitted on a comparison graph
//! by weighted least squares (one node fixed as anchor), and calibration
//! uncertainty is bootstrapped into per-cohort credible intervals and
//! isolate-level probabilistic labels over a threshold grid.
//!
//! Zero external dependencies.

mod calib;
mod cutoff;
mod dashboard;
mod linalg;
mod rng;
mod trace;

use calib::{build_edges, fit, Pair, Params};
use std::collections::HashMap;
use std::f64::consts::LOG10_2;
use std::fs;
use std::process::exit;

struct Args {
    phenotypes: String,
    pairs: String,
    out_prefix: String,
    anchor: Option<String>,
    already_log: bool,
    min_support: usize,
    max_drift: f64,
    bootstrap: usize,
    seed: u64,
    thresholds: Vec<f64>,
    sigma_mode: String, // "fixed" | "empirical"
    drift_scale: u8,    // 0=doubling 1=rms 2=mean 3=max
    lambda: f64,
    robust: bool,
    nu: f64,
    dashboard: Option<String>,
    trace: Option<String>,
}

fn usage() -> ! {
    eprintln!(
        "phenocal — genome-anchored phenotype calibration

USAGE:
    phenocal --phenotypes <pheno.tsv> --pairs <pairs.tsv> --out <prefix> [OPTIONS]

REQUIRED:
    --phenotypes <FILE>   TSV with header columns: sample  value  cohort
    --pairs <FILE>        TSV with header columns: sample_i  sample_j  distance
                          (genomic distance; within- and cross-cohort allowed,
                           only cross-cohort pairs are used)
    --out <PREFIX>        Output prefix. Writes:
                            <prefix>.offsets.tsv    cohort offsets + 95% CrI
                            <prefix>.sigma.tsv      per-cohort resolution sigma_c
                            <prefix>.edges.tsv      per cohort-pair: d_min, tau, n, fold, SE, weight
                            <prefix>.harmonised.tsv per-sample raw + harmonised value
                            <prefix>.labels.tsv     per-sample P(tolerant) per threshold
                            <prefix>.cutoff.tsv     data-driven cutoff (KDE antimode + GMM crossover)

OPTIONS:
    --anchor <COHORT>     Cohort fixed to offset 0 (default: largest cohort)
    --already-log         'value' is already log10 (default: raw, log10 applied)
    --min-support <N>     Min near-clonal pairs per edge (default 5)
    --max-drift <D>       Max median drift in dilutions for tau selection (default 1.0)
    --bootstrap <B>       Bootstrap draws for credible intervals (default 10000)
    --seed <S>            RNG seed (default 20260603)
    --thresholds <list>   Comma-separated tolerance thresholds on the RAW scale
                          (default 0.75,1.0,1.25,1.5,2.0)
    --sigma <MODE>        Per-cohort noise model: 'fixed' (= log10(2), classic) or
                          'empirical' (estimated from each cohort's MIC grid spacing).
                          Default 'fixed'.
    --drift-scale <S>     Units of the tau drift tolerance: 'doubling' (log10 2),
                          'rms' (default), 'mean', or 'max' of the cohort resolutions.
                          Not load-bearing (offsets are essentially invariant).
    --lambda <L>          Biological-drift variance per unit genetic distance in the
                          edge SE (default log10(2)^2). SE_ab = sqrt((s_a^2+s_b^2+L*tau)/n).
    --robust              Robust bootstrap: sample edge perturbations from Student-t
                          instead of Normal (heavier tails, tolerates outliers).
    --nu <DF>             Degrees of freedom for --robust (default 4).
    --dashboard <FILE>    Also write a self-contained interactive HTML report
                          (the comparison graph, every twin pair and log-ratio per
                          edge with the median highlighted, offsets + credible
                          intervals, before/after distributions, per-isolate
                          probabilities). Opens in any browser, no dependencies.
    --trace <FILE>        Also write a step-by-step mathematical trace of the run
                          in Markdown (inputs, per-edge twins + tau selection,
                          WLS normal-equation matrices, bootstrap CIs, cutoff).
                          Renders directly on GitHub; no compilation needed.
    -h, --help            Show this help
"
    );
    exit(2);
}

fn parse_args() -> Args {
    let mut a = Args {
        phenotypes: String::new(),
        pairs: String::new(),
        out_prefix: String::new(),
        anchor: None,
        already_log: false,
        min_support: 5,
        max_drift: 1.0,
        bootstrap: 10000,
        seed: 20260603,
        thresholds: vec![0.75, 1.0, 1.25, 1.5, 2.0],
        sigma_mode: "fixed".to_string(),
        drift_scale: 1, // rms
        lambda: LOG10_2.powi(2), // log10(2)^2 -> reproduces the classic SE
        robust: false,
        nu: 4.0,
        dashboard: None,
        trace: None,
    };
    let mut it = std::env::args().skip(1);
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "--phenotypes" => a.phenotypes = it.next().unwrap_or_else(|| usage()),
            "--pairs" => a.pairs = it.next().unwrap_or_else(|| usage()),
            "--out" => a.out_prefix = it.next().unwrap_or_else(|| usage()),
            "--anchor" => a.anchor = Some(it.next().unwrap_or_else(|| usage())),
            "--already-log" => a.already_log = true,
            "--min-support" => {
                a.min_support = it.next().unwrap_or_else(|| usage()).parse().unwrap_or_else(|_| usage())
            }
            "--max-drift" => {
                a.max_drift = it.next().unwrap_or_else(|| usage()).parse().unwrap_or_else(|_| usage())
            }
            "--bootstrap" => {
                a.bootstrap = it.next().unwrap_or_else(|| usage()).parse().unwrap_or_else(|_| usage())
            }
            "--seed" => {
                a.seed = it.next().unwrap_or_else(|| usage()).parse().unwrap_or_else(|_| usage())
            }
            "--thresholds" => {
                a.thresholds = it
                    .next()
                    .unwrap_or_else(|| usage())
                    .split(',')
                    .map(|s| s.trim().parse().unwrap_or_else(|_| usage()))
                    .collect()
            }
            "--sigma" => a.sigma_mode = it.next().unwrap_or_else(|| usage()),
            "--drift-scale" => {
                a.drift_scale = match it.next().unwrap_or_else(|| usage()).as_str() {
                    "doubling" => 0,
                    "rms" => 1,
                    "mean" => 2,
                    "max" => 3,
                    other => {
                        eprintln!("Unknown --drift-scale '{other}' (use doubling|rms|mean|max).");
                        exit(1);
                    }
                }
            }
            "--lambda" => {
                a.lambda = it.next().unwrap_or_else(|| usage()).parse().unwrap_or_else(|_| usage())
            }
            "--robust" => a.robust = true,
            "--nu" => a.nu = it.next().unwrap_or_else(|| usage()).parse().unwrap_or_else(|_| usage()),
            "--dashboard" => a.dashboard = Some(it.next().unwrap_or_else(|| usage())),
            "--trace" => a.trace = Some(it.next().unwrap_or_else(|| usage())),
            "-h" | "--help" => usage(),
            _ => {
                eprintln!("Unknown argument: {arg}");
                usage();
            }
        }
    }
    if a.phenotypes.is_empty() || a.pairs.is_empty() || a.out_prefix.is_empty() {
        usage();
    }
    a
}

fn read_lines(path: &str) -> Vec<String> {
    fs::read_to_string(path)
        .unwrap_or_else(|e| {
            eprintln!("Cannot read {path}: {e}");
            exit(1);
        })
        .lines()
        .map(|s| s.to_string())
        .collect()
}

/// Find a column index by (case-insensitive) name among accepted aliases.
fn col(header: &[&str], aliases: &[&str]) -> usize {
    for (i, h) in header.iter().enumerate() {
        let hl = h.trim().to_lowercase();
        if aliases.iter().any(|a| *a == hl) {
            return i;
        }
    }
    eprintln!("Required column not found (one of {aliases:?}) in header {header:?}");
    exit(1);
}

struct Pheno {
    log_value: HashMap<String, f64>,
    raw_value: HashMap<String, f64>,
    cohort: HashMap<String, String>,
    cohorts: Vec<String>, // unique, stable order
    counts: HashMap<String, usize>,
}

fn read_phenotypes(path: &str, already_log: bool) -> Pheno {
    let lines = read_lines(path);
    let header: Vec<&str> = lines[0].split('\t').collect();
    let cs = col(&header, &["sample", "accession", "id"]);
    let cv = col(&header, &["value", "mic", "phenotype"]);
    let cc = col(&header, &["cohort", "author", "batch", "study"]);
    let mut p = Pheno {
        log_value: HashMap::new(),
        raw_value: HashMap::new(),
        cohort: HashMap::new(),
        cohorts: Vec::new(),
        counts: HashMap::new(),
    };
    let mut n_missing = 0usize;
    for line in &lines[1..] {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() <= cs.max(cv).max(cc) {
            continue;
        }
        let sample = f[cs].trim().to_string();
        let vstr = f[cv].trim();
        let raw: f64 = match vstr.parse() {
            Ok(v) => v,
            Err(_) => {
                // distinguish genuine missing values from a categorical column
                let low = vstr.to_lowercase();
                if vstr.is_empty()
                    || matches!(low.as_str(), "na" | "nan" | "null" | "none" | "-" | ".")
                {
                    n_missing += 1;
                    continue;
                }
                eprintln!(
                    "Error: phenotype 'value' must be numeric (a continuous measurement), \
                     but found '{vstr}' for sample '{sample}'.\n\
                     phenocal calibrates a continuous phenotype (e.g. MIC, OD, IC50); a \
                     categorical/string label (e.g. 'tolerant', 'R') cannot be calibrated \
                     directly.\n\
                     If your phenotype is categorical, calibrate the underlying continuous \
                     measurement and apply the threshold afterwards (the categorical label \
                     is an output, not an input)."
                );
                exit(1);
            }
        };
        let cohort = f[cc].trim().to_string();
        if raw <= 0.0 && !already_log {
            continue; // log10 of non-positive undefined
        }
        let logv = if already_log { raw } else { raw.log10() };
        p.raw_value
            .insert(sample.clone(), if already_log { 10f64.powf(raw) } else { raw });
        p.log_value.insert(sample.clone(), logv);
        if !p.counts.contains_key(&cohort) {
            p.cohorts.push(cohort.clone());
        }
        *p.counts.entry(cohort.clone()).or_insert(0) += 1;
        p.cohort.insert(sample, cohort);
    }
    if n_missing > 0 {
        eprintln!("note: skipped {n_missing} row(s) with missing/NA phenotype value");
    }
    p
}

fn main() {
    let args = parse_args();
    let ph = read_phenotypes(&args.phenotypes, args.already_log);
    if ph.cohorts.len() < 2 {
        eprintln!("Need at least 2 cohorts; found {}.", ph.cohorts.len());
        exit(1);
    }

    let cidx: HashMap<String, usize> = ph
        .cohorts
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, c)| (c, i))
        .collect();

    let anchor_name = args.anchor.clone().unwrap_or_else(|| {
        ph.cohorts
            .iter()
            .max_by_key(|c| ph.counts.get(*c).copied().unwrap_or(0))
            .unwrap()
            .clone()
    });
    let anchor = *cidx.get(&anchor_name).unwrap_or_else(|| {
        eprintln!("Anchor cohort '{anchor_name}' not present in data.");
        exit(1);
    });

    // read pairs -> cross-cohort Pair list
    let lines = read_lines(&args.pairs);
    let header: Vec<&str> = lines[0].split('\t').collect();
    let ci = col(&header, &["sample_i", "acc_i", "i"]);
    let cj = col(&header, &["sample_j", "acc_j", "j"]);
    let cd = col(&header, &["distance", "dist", "d"]);
    let mut pairs: Vec<Pair> = Vec::new();
    let mut n_cross = 0usize;
    for line in &lines[1..] {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() <= ci.max(cj).max(cd) {
            continue;
        }
        let si = f[ci].trim();
        let sj = f[cj].trim();
        let dist: f64 = match f[cd].trim().parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let (ca, cb) = match (ph.cohort.get(si), ph.cohort.get(sj)) {
            (Some(a), Some(b)) => (a, b),
            _ => continue,
        };
        if ca == cb {
            continue; // within-cohort: biology does not cancel cleanly
        }
        let (yi, yj) = match (ph.log_value.get(si), ph.log_value.get(sj)) {
            (Some(a), Some(b)) => (*a, *b),
            _ => continue,
        };
        pairs.push(Pair {
            a: cidx[ca],
            b: cidx[cb],
            dist,
            signed: yi - yj,
            si: si.to_string(),
            sj: sj.to_string(),
            va: yi,
            vb: yj,
        });
        n_cross += 1;
    }
    if n_cross == 0 {
        eprintln!("No usable cross-cohort pairs found.");
        exit(1);
    }

    // per-cohort noise scale sigma_c
    let sigma_fixed = LOG10_2 / std::f64::consts::SQRT_2; // log10(2)/sqrt(2)
    let sigma: Vec<f64> = match args.sigma_mode.as_str() {
        "fixed" => vec![sigma_fixed; ph.cohorts.len()],
        "empirical" => ph
            .cohorts
            .iter()
            .map(|c| empirical_sigma(&ph, c).unwrap_or(sigma_fixed))
            .collect(),
        other => {
            eprintln!("Unknown --sigma mode '{other}' (use 'fixed' or 'empirical').");
            exit(1);
        }
    };

    let params = Params {
        min_support: args.min_support,
        max_drift_dilutions: args.max_drift,
        bootstrap: args.bootstrap,
        seed: args.seed,
        sigma: sigma.clone(),
        lambda: args.lambda,
        robust: args.robust,
        nu: args.nu,
        drift_scale: args.drift_scale,
    };

    let edges = build_edges(&pairs, ph.cohorts.len(), &params);
    let sol = fit(&edges, ph.cohorts.len(), anchor, &params).unwrap_or_else(|e| {
        eprintln!("Calibration failed: {e}");
        exit(1);
    });

    // ---- offsets ----
    let mut off = String::from("cohort\tn\tdelta_log10\tfold\tfold_lo95\tfold_hi95\tanchor\n");
    for (i, c) in ph.cohorts.iter().enumerate() {
        off.push_str(&format!(
            "{}\t{}\t{:.5}\t{:.4}\t{:.4}\t{:.4}\t{}\n",
            c,
            ph.counts[c],
            sol.delta[i],
            10f64.powf(sol.delta[i]),
            10f64.powf(sol.lo95[i]),
            10f64.powf(sol.hi95[i]),
            if i == anchor { "yes" } else { "no" },
        ));
    }
    write(&format!("{}.offsets.tsv", args.out_prefix), &off);

    // ---- per-cohort measurement resolution sigma_c ----
    let mut sg = String::from("cohort\tn\tsigma_log10\tsigma_dilutions\tmode\n");
    for (i, c) in ph.cohorts.iter().enumerate() {
        sg.push_str(&format!(
            "{}\t{}\t{:.5}\t{:.3}\t{}\n",
            c,
            ph.counts[c],
            sigma[i],
            sigma[i] / LOG10_2,
            args.sigma_mode,
        ));
    }
    write(&format!("{}.sigma.tsv", args.out_prefix), &sg);

    // ---- edges (transparency: d_min / tau / n / fold) ----
    let mut et =
        String::from("cohort_a\tcohort_b\td_min\ttau\tn_pairs\tfold_a_over_b\tse_log10\tweight\n");
    for e in &edges {
        et.push_str(&format!(
            "{}\t{}\t{:.3}\t{:.3}\t{}\t{:.4}\t{:.4}\t{:.3}\n",
            ph.cohorts[e.a],
            ph.cohorts[e.b],
            e.d_min,
            e.tau,
            e.n,
            10f64.powf(e.delta),
            e.se,
            e.weight,
        ));
    }
    write(&format!("{}.edges.tsv", args.out_prefix), &et);

    // ---- harmonised values ----
    let mut samples: Vec<&String> = ph.log_value.keys().collect();
    samples.sort();
    let mut hv = String::from("sample\tcohort\tvalue_raw\tvalue_harmonised\n");
    for s in &samples {
        let c = &ph.cohort[*s];
        let i = cidx[c];
        let harm_log = ph.log_value[*s] - sol.delta[i];
        hv.push_str(&format!(
            "{}\t{}\t{:.5}\t{:.5}\n",
            s,
            c,
            ph.raw_value[*s],
            10f64.powf(harm_log),
        ));
    }
    write(&format!("{}.harmonised.tsv", args.out_prefix), &hv);

    // ---- probabilistic labels over the threshold grid ----
    // Propagate calibration uncertainty: per-cohort offset posterior approximated
    // as Normal with sd from the 95% CrI half-width; recompute label per draw.
    let mut rng = rng::Rng::new(args.seed ^ 0xABCDEF);
    let b = params.bootstrap.max(1);
    let sd: Vec<f64> = (0..ph.cohorts.len())
        .map(|i| ((sol.hi95[i] - sol.lo95[i]) / (2.0 * 1.959964)).max(0.0))
        .collect();

    let mut lab = String::from("sample\tcohort\tvalue_harmonised");
    for t in &args.thresholds {
        lab.push_str(&format!("\tP_tol_{}", fmt_thr(*t)));
    }
    lab.push('\n');
    let mut iso_rows: Vec<dashboard::Iso> = Vec::with_capacity(samples.len());
    for s in &samples {
        let c = &ph.cohort[*s];
        let i = cidx[c];
        let y = ph.log_value[*s];
        let harm = 10f64.powf(y - sol.delta[i]);
        let mut counts = vec![0usize; args.thresholds.len()];
        for _ in 0..b {
            let di = sol.delta[i] + sd[i] * rng.next_normal();
            let hvd = 10f64.powf(y - di);
            for (k, t) in args.thresholds.iter().enumerate() {
                if hvd >= *t {
                    counts[k] += 1;
                }
            }
        }
        let probs: Vec<f64> = counts.iter().map(|&c| c as f64 / b as f64).collect();
        lab.push_str(&format!("{}\t{}\t{:.5}", s, c, harm));
        for pr in &probs {
            lab.push_str(&format!("\t{:.4}", pr));
        }
        lab.push('\n');
        iso_rows.push(dashboard::Iso {
            sample: (*s).clone(),
            cohort: i,
            raw: ph.raw_value[*s],
            harm,
            probs,
        });
    }
    write(&format!("{}.labels.tsv", args.out_prefix), &lab);

    // ---- data-driven cutoff on the harmonised scale (KDE antimode + GMM crossover) ----
    let logharm: Vec<f64> = iso_rows
        .iter()
        .filter(|r| r.harm > 0.0)
        .map(|r| r.harm.log10())
        .collect();
    let cut = cutoff::estimate(&logharm, params.bootstrap, args.seed);
    let fmt = |o: Option<f64>| o.map(|x| format!("{:.4}", x)).unwrap_or_else(|| "NA".into());
    let fmt_mgl = |o: Option<f64>| o.map(|x| format!("{:.4}", 10f64.powf(x))).unwrap_or_else(|| "NA".into());
    let nan_mgl = |x: f64| if x.is_finite() { format!("{:.4}", 10f64.powf(x)) } else { "NA".into() };
    let mut ct = String::from("method\tcutoff_log10\tcutoff_mgL\tlo95_mgL\thi95_mgL\n");
    ct.push_str(&format!(
        "kde_antimode\t{}\t{}\tNA\tNA\n",
        fmt(cut.kde_antimode),
        fmt_mgl(cut.kde_antimode)
    ));
    ct.push_str(&format!(
        "gmm_crossover\t{}\t{}\t{}\t{}\n",
        fmt(cut.gmm_crossover),
        fmt_mgl(cut.gmm_crossover),
        nan_mgl(cut.gmm_lo),
        nan_mgl(cut.gmm_hi)
    ));
    write(&format!("{}.cutoff.tsv", args.out_prefix), &ct);
    eprintln!(
        "data-driven cutoff (harmonised, mg/L): KDE antimode {}, GMM crossover {} [{}, {}] \
         (mixture modes ~{:.2} vs ~{:.2})",
        fmt_mgl(cut.kde_antimode),
        fmt_mgl(cut.gmm_crossover),
        nan_mgl(cut.gmm_lo),
        nan_mgl(cut.gmm_hi),
        10f64.powf(cut.comp_lo.1),
        10f64.powf(cut.comp_hi.1)
    );

    // ---- summary ----
    eprintln!(
        "phenocal: {} cohorts, {} cross-cohort pairs, {} edges, anchor = {} | sigma={} lambda={:.4}{}",
        ph.cohorts.len(),
        n_cross,
        edges.len(),
        anchor_name,
        args.sigma_mode,
        args.lambda,
        if args.robust { format!(" robust-t(nu={})", args.nu) } else { String::new() }
    );
    eprintln!("per-cohort noise scale sigma_c (log10 units):");
    for (i, c) in ph.cohorts.iter().enumerate() {
        eprintln!("  {:<14} sigma = {:.3}", c, sigma[i]);
    }
    eprintln!(
        "weighted residual RMSE = {:.4} log10 ({:.2} dilutions)",
        sol.rmse,
        sol.rmse / LOG10_2
    );
    eprintln!("cohort offsets (fold vs anchor, 95% CrI):");
    for (i, c) in ph.cohorts.iter().enumerate() {
        eprintln!(
            "  {:<14} {:>7.2}x  [{:>6.2}, {:>6.2}]",
            c,
            10f64.powf(sol.delta[i]),
            10f64.powf(sol.lo95[i]),
            10f64.powf(sol.hi95[i])
        );
    }
    eprintln!(
        "wrote {0}.offsets.tsv {0}.sigma.tsv {0}.edges.tsv {0}.harmonised.tsv {0}.labels.tsv {0}.cutoff.tsv",
        args.out_prefix
    );

    // ---- sensitivity to our SE/drift formula choices (for the trace) ----
    // Re-fit the offsets under alternative resolution (sigma) and drift-scale
    // forms, to show these modelling choices are not load-bearing.
    let sensitivity: Vec<(String, Vec<f64>)> = if args.trace.is_some() {
        let ncoh = ph.cohorts.len();
        let s_fixed = vec![sigma_fixed; ncoh];
        let s_emp: Vec<f64> = ph
            .cohorts
            .iter()
            .map(|c| empirical_sigma(&ph, c).unwrap_or(sigma_fixed))
            .collect();
        let variants: [(&str, &Vec<f64>, u8); 5] = [
            ("sigma=fixed, drift=rms", &s_fixed, 1),
            ("sigma=empirical, drift=rms", &s_emp, 1),
            ("sigma=empirical, drift=doubling", &s_emp, 0),
            ("sigma=empirical, drift=mean", &s_emp, 2),
            ("sigma=empirical, drift=max", &s_emp, 3),
        ];
        let mut out = Vec::new();
        for (label, sig, ds) in variants {
            let p = Params {
                min_support: args.min_support,
                max_drift_dilutions: args.max_drift,
                bootstrap: args.bootstrap.min(1000),
                seed: args.seed,
                sigma: sig.clone(),
                lambda: args.lambda,
                robust: args.robust,
                nu: args.nu,
                drift_scale: ds,
            };
            let edg = build_edges(&pairs, ncoh, &p);
            if let Ok(so) = fit(&edg, ncoh, anchor, &p) {
                out.push((
                    label.to_string(),
                    (0..ncoh).map(|i| 10f64.powf(so.delta[i])).collect(),
                ));
            }
        }
        out
    } else {
        Vec::new()
    };

    // ---- optional reports (dashboard HTML and/or mathematical trace in Markdown) ----
    if args.dashboard.is_some() || args.trace.is_some() {
        let ctx = dashboard::Ctx {
            cohorts: &ph.cohorts,
            counts: &ph.counts,
            anchor,
            anchor_name: &anchor_name,
            edges: &edges,
            sol: &sol,
            sigma: &sigma,
            sigma_mode: &args.sigma_mode,
            lambda: args.lambda,
            robust: args.robust,
            nu: args.nu,
            min_support: args.min_support,
            max_drift: args.max_drift,
            bootstrap: args.bootstrap,
            n_cross,
            thresholds: &args.thresholds,
            isolates: &iso_rows,
            cutoff: &cut,
            logharm: &logharm,
            sensitivity: &sensitivity,
        };
        if let Some(path) = &args.dashboard {
            write(path, &dashboard::render(&ctx));
            eprintln!("wrote interactive dashboard: {path}");
        }
        if let Some(path) = &args.trace {
            trace::write_report(&ctx, path);
            eprintln!("wrote mathematical trace (Markdown + SVG figures): {path}");
        }
    }
}

/// Estimate a cohort's log10-MIC resolution sigma_c from its observed MIC grid:
/// take the distinct log10 values, merge rounding duplicates (within 0.05 log10),
/// and return the median spacing between consecutive rungs (floored at 0.05).
/// Returns None if fewer than two distinct rungs.
fn empirical_sigma(ph: &Pheno, cohort: &str) -> Option<f64> {
    let mut logs: Vec<f64> = ph
        .log_value
        .iter()
        .filter(|(s, _)| ph.cohort.get(*s).map(|c| c == cohort).unwrap_or(false))
        .map(|(_, v)| *v)
        .collect();
    if logs.len() < 2 {
        return None;
    }
    logs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    // merge rounding duplicates into rungs
    let tol = 0.05;
    let mut rungs: Vec<f64> = vec![logs[0]];
    for &v in &logs[1..] {
        if (v - rungs.last().unwrap()).abs() > tol {
            rungs.push(v);
        }
    }
    if rungs.len() < 2 {
        return None;
    }
    let spacings: Vec<f64> = rungs.windows(2).map(|w| w[1] - w[0]).collect();
    let s = crate::linalg::median(&spacings);
    Some(s.max(0.05))
}

fn fmt_thr(t: f64) -> String {
    format!("{t}").replace('.', "p")
}

fn write(path: &str, content: &str) {
    if let Err(e) = fs::write(path, content) {
        eprintln!("Cannot write {path}: {e}");
        exit(1);
    }
}
