//! Genome-anchored calibration: build the cohort comparison graph from
//! near-clonal cross-cohort pairs, solve weighted least squares for per-cohort
//! offsets, and bootstrap credible intervals. Mirrors Algorithm 1 of the paper.

use crate::linalg::{median, quantile, solve};
use crate::rng::Rng;
use std::collections::HashMap;

const LOG10_2: f64 = 0.3010299956639812;

/// One cross-cohort near-clonal pair (already oriented `a` vs `b`).
pub struct Pair {
    pub a: usize,        // cohort index (oriented "high" side of the signed value)
    pub b: usize,        // cohort index
    pub dist: f64,       // genomic distance between the two isolates
    pub signed: f64,     // y_a - y_b  (log10 phenotype difference, a over b)
}

/// A fitted edge between two cohorts.
pub struct Edge {
    pub a: usize,
    pub b: usize,
    pub d_min: f64,
    pub tau: f64,
    pub n: usize,
    pub delta: f64, // median signed log10 (a over b)
    pub se: f64,
    pub weight: f64,
}

pub struct Params {
    pub min_support: usize,
    pub max_drift_dilutions: f64,
    pub bootstrap: usize,
    pub seed: u64,
    /// Per-cohort measurement-noise scale sigma_c (log10 units). Length = n_cohorts.
    pub sigma: Vec<f64>,
    /// Biological-drift variance per unit genetic distance (lambda).
    pub lambda: f64,
    /// Robust bootstrap: sample edge perturbations from Student-t instead of Normal.
    pub robust: bool,
    /// Degrees of freedom for the robust t-likelihood.
    pub nu: f64,
}

/// Adaptive tau selection + edge statistics for one cohort pair.
/// `pairs` are all observed cross-cohort pairs for this (a,b), oriented a over b.
fn build_edge(
    a: usize,
    b: usize,
    sigma_a: f64,
    sigma_b: f64,
    mut pairs: Vec<(f64, f64)>,
    p: &Params,
) -> Option<Edge> {
    if pairs.is_empty() {
        return None;
    }
    // sort by distance
    pairs.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
    let d_min = pairs[0].0;
    let base_vals: Vec<f64> = pairs.iter().filter(|x| x.0 == d_min).map(|x| x.1).collect();
    let base = median(&base_vals);

    // candidate taus = the distinct distances, ascending
    let mut taus: Vec<f64> = pairs.iter().map(|x| x.0).collect();
    taus.dedup();

    let mut chosen_tau = *taus.last().unwrap();
    let mut chosen_vals: Vec<f64> = pairs.iter().map(|x| x.1).collect();
    let mut satisfied = false;
    for &tau in &taus {
        let vals: Vec<f64> = pairs.iter().filter(|x| x.0 <= tau).map(|x| x.1).collect();
        let drift = (median(&vals) - base).abs() / LOG10_2;
        if vals.len() >= p.min_support && drift <= p.max_drift_dilutions {
            chosen_tau = tau;
            chosen_vals = vals;
            satisfied = true;
            break;
        }
    }
    // If we never reached min_support, fall back to all pairs (kept, but down-weighted
    // through the large tau / small n in the SE below).
    if !satisfied {
        chosen_vals = pairs.iter().map(|x| x.1).collect();
        chosen_tau = *taus.last().unwrap();
    }

    let n = chosen_vals.len();
    let delta = median(&chosen_vals);
    // Protocol-agnostic edge SE: cohort-specific resolution sigma_a, sigma_b plus
    // a biological-drift variance growing with genetic distance (lambda * tau),
    // averaged over the n supporting pairs. With sigma = log10(2)/sqrt(2) and
    // lambda = log10(2)^2 this reduces exactly to log10(2)*sqrt((1+tau)/n).
    let se = ((sigma_a * sigma_a + sigma_b * sigma_b + p.lambda * chosen_tau) / n as f64).sqrt();
    let weight = 1.0 / (se * se);
    Some(Edge {
        a,
        b,
        d_min,
        tau: chosen_tau,
        n,
        delta,
        se,
        weight,
    })
}

/// Build all edges from the full list of cross-cohort pairs.
pub fn build_edges(pairs: &[Pair], n_cohorts: usize, p: &Params) -> Vec<Edge> {
    // group by unordered cohort pair, orient consistently (lower index = a)
    let mut groups: HashMap<(usize, usize), Vec<(f64, f64)>> = HashMap::new();
    for pr in pairs {
        let (a, b, s) = if pr.a < pr.b {
            (pr.a, pr.b, pr.signed)
        } else {
            (pr.b, pr.a, -pr.signed)
        };
        groups.entry((a, b)).or_default().push((pr.dist, s));
    }
    let mut edges = Vec::new();
    for ((a, b), v) in groups {
        if a < n_cohorts && b < n_cohorts {
            let sa = p.sigma.get(a).copied().unwrap_or(LOG10_2 / std::f64::consts::SQRT_2);
            let sb = p.sigma.get(b).copied().unwrap_or(LOG10_2 / std::f64::consts::SQRT_2);
            if let Some(e) = build_edge(a, b, sa, sb, v, p) {
                edges.push(e);
            }
        }
    }
    // stable order for reproducible output
    edges.sort_by(|x, y| (x.a, x.b).cmp(&(y.a, y.b)));
    edges
}

/// Solve weighted least squares for the free cohorts given edge deltas.
/// `anchor` is fixed to 0. `free` maps cohort index -> column index (anchor absent).
fn wls(
    edges: &[Edge],
    deltas: &[f64],
    free_col: &HashMap<usize, usize>,
    k: usize,
) -> Option<Vec<f64>> {
    // Normal equations (A^T W A) x = A^T W d, where each edge row has +1 at a, -1 at b.
    let mut ata = vec![0.0; k * k];
    let mut atb = vec![0.0; k];
    for (e, &d) in edges.iter().zip(deltas) {
        let w = e.weight;
        let ca = free_col.get(&e.a).copied();
        let cb = free_col.get(&e.b).copied();
        // row contributions
        if let Some(i) = ca {
            ata[i * k + i] += w;
            atb[i] += w * d;
            if let Some(j) = cb {
                ata[i * k + j] -= w;
            }
        }
        if let Some(j) = cb {
            ata[j * k + j] += w;
            atb[j] -= w * d;
            if let Some(i) = ca {
                ata[j * k + i] -= w;
            }
        }
    }
    solve(&ata, &atb, k)
}

pub struct Solution {
    pub delta: Vec<f64>,      // per cohort (anchor = 0)
    pub lo95: Vec<f64>,
    pub hi95: Vec<f64>,
    pub rmse: f64,
}

/// Fit offsets + bootstrap credible intervals.
pub fn fit(
    edges: &[Edge],
    n_cohorts: usize,
    anchor: usize,
    p: &Params,
) -> Result<Solution, String> {
    // free cohorts = all except anchor
    let mut free_col = HashMap::new();
    let mut free_list = Vec::new();
    for c in 0..n_cohorts {
        if c != anchor {
            free_col.insert(c, free_list.len());
            free_list.push(c);
        }
    }
    let k = free_list.len();

    let point_deltas: Vec<f64> = edges.iter().map(|e| e.delta).collect();
    let sol = wls(edges, &point_deltas, &free_col, k)
        .ok_or_else(|| "Calibration graph is singular (disconnected from anchor?).".to_string())?;

    // assemble full delta vector
    let to_full = |x: &[f64]| -> Vec<f64> {
        let mut full = vec![0.0; n_cohorts];
        for (c, &col) in free_col.iter() {
            full[*c] = x[col];
        }
        full
    };
    let delta_full = to_full(&sol);

    // weighted RMSE of residuals
    let mut num = 0.0;
    let mut den = 0.0;
    for e in edges {
        let pred = delta_full[e.a] - delta_full[e.b];
        let r = e.delta - pred;
        num += e.weight * r * r;
        den += e.weight;
    }
    let rmse = (num / den).sqrt();

    // bootstrap
    let mut rng = Rng::new(p.seed);
    let mut draws: Vec<Vec<f64>> = vec![Vec::with_capacity(p.bootstrap); k];
    for _ in 0..p.bootstrap {
        let sampled: Vec<f64> = edges
            .iter()
            .map(|e| {
                let z = if p.robust {
                    rng.next_t(p.nu)
                } else {
                    rng.next_normal()
                };
                e.delta + e.se * z
            })
            .collect();
        if let Some(x) = wls(edges, &sampled, &free_col, k) {
            for col in 0..k {
                draws[col].push(x[col]);
            }
        }
    }

    let mut lo = vec![0.0; n_cohorts];
    let mut hi = vec![0.0; n_cohorts];
    let mut med = delta_full.clone();
    for (c, &col) in free_col.iter() {
        let mut d = draws[col].clone();
        d.sort_by(|a, b| a.partial_cmp(b).unwrap());
        lo[*c] = quantile(&d, 0.025);
        hi[*c] = quantile(&d, 0.975);
        med[*c] = quantile(&d, 0.5);
    }
    Ok(Solution {
        delta: med,
        lo95: lo,
        hi95: hi,
        rmse,
    })
}
