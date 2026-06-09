//! Genome-anchored calibration: build the cohort comparison graph from
//! near-clonal cross-cohort pairs, solve weighted least squares for per-cohort
//! offsets, and bootstrap credible intervals. Mirrors Algorithm 1 of the paper.

use crate::linalg::{median, quantile, solve};
use crate::rng::Rng;
use std::collections::HashMap;

use std::f64::consts::LOG10_2;

/// One cross-cohort near-clonal pair (already oriented `a` vs `b`).
pub struct Pair {
    pub a: usize,    // cohort index of sample si
    pub b: usize,    // cohort index of sample sj
    pub dist: f64,   // genomic distance between the two isolates
    pub signed: f64, // y_a - y_b  (log10 phenotype difference, a over b)
    pub si: String,  // sample on the a side
    pub sj: String,  // sample on the b side
    pub va: f64,     // log10 phenotype of si
    pub vb: f64,     // log10 phenotype of sj
}

/// One supporting twin pair for an edge (oriented a over b), kept for inspection.
#[derive(Clone)]
pub struct EdgePair {
    pub si: String,
    pub sj: String,
    pub dist: f64,
    pub va: f64,     // log10 value, a side
    pub vb: f64,     // log10 value, b side
    pub signed: f64, // va - vb
    pub used: bool,  // within the chosen tau (entered the median)
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
    pub cand: Vec<EdgePair>, // every cross-cohort candidate (sorted by distance), used-flagged
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
    /// Drift-tolerance scale (units of the τ-drift bound):
    /// 0 = doubling (log10 2), 1 = rms of σ (default), 2 = mean of σ, 3 = max of σ.
    pub drift_scale: u8,
}

/// Adaptive tau selection + edge statistics for one cohort pair.
/// `pairs` are all observed cross-cohort pairs for this (a,b), oriented a over b.
fn build_edge(
    a: usize,
    b: usize,
    sigma_a: f64,
    sigma_b: f64,
    mut pairs: Vec<EdgePair>,
    p: &Params,
) -> Option<Edge> {
    if pairs.is_empty() {
        return None;
    }
    // sort by distance
    pairs.sort_by(|x, y| x.dist.partial_cmp(&y.dist).unwrap());
    let d_min = pairs[0].dist;
    let base_vals: Vec<f64> = pairs
        .iter()
        .filter(|x| x.dist == d_min)
        .map(|x| x.signed)
        .collect();
    let base = median(&base_vals);

    // candidate taus = the distinct distances, ascending
    let mut taus: Vec<f64> = pairs.iter().map(|x| x.dist).collect();
    taus.dedup();

    let mut chosen_tau = *taus.last().unwrap();
    let mut chosen_vals: Vec<f64> = pairs.iter().map(|x| x.signed).collect();
    let mut satisfied = false;
    // drift tolerance scale (does not assume a two-fold grid unless drift_scale=0):
    // reduces to log10(2) when both cohorts are two-fold (sigma=log10 2).
    let edge_res = match p.drift_scale {
        0 => LOG10_2,
        2 => (sigma_a + sigma_b) / 2.0,
        3 => sigma_a.max(sigma_b),
        _ => ((sigma_a * sigma_a + sigma_b * sigma_b) / 2.0).sqrt(),
    }
    .max(1e-9);
    for &tau in &taus {
        let vals: Vec<f64> = pairs
            .iter()
            .filter(|x| x.dist <= tau)
            .map(|x| x.signed)
            .collect();
        let drift = (median(&vals) - base).abs() / edge_res;
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
        chosen_vals = pairs.iter().map(|x| x.signed).collect();
        chosen_tau = *taus.last().unwrap();
    }

    // flag which candidates entered the median
    for pr in pairs.iter_mut() {
        pr.used = pr.dist <= chosen_tau;
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
        cand: pairs,
    })
}

/// Build all edges from the full list of cross-cohort pairs.
pub fn build_edges(pairs: &[Pair], n_cohorts: usize, p: &Params) -> Vec<Edge> {
    // group by unordered cohort pair, orient consistently (lower index = a)
    let mut groups: HashMap<(usize, usize), Vec<EdgePair>> = HashMap::new();
    for pr in pairs {
        let ep = if pr.a < pr.b {
            EdgePair {
                si: pr.si.clone(),
                sj: pr.sj.clone(),
                dist: pr.dist,
                va: pr.va,
                vb: pr.vb,
                signed: pr.signed,
                used: false,
            }
        } else {
            EdgePair {
                si: pr.sj.clone(),
                sj: pr.si.clone(),
                dist: pr.dist,
                va: pr.vb,
                vb: pr.va,
                signed: -pr.signed,
                used: false,
            }
        };
        let (a, b) = if pr.a < pr.b {
            (pr.a, pr.b)
        } else {
            (pr.b, pr.a)
        };
        groups.entry((a, b)).or_default().push(ep);
    }
    let mut edges = Vec::new();
    for ((a, b), v) in groups {
        if a < n_cohorts && b < n_cohorts {
            let sa = p
                .sigma
                .get(a)
                .copied()
                .unwrap_or(LOG10_2 / std::f64::consts::SQRT_2);
            let sb = p
                .sigma
                .get(b)
                .copied()
                .unwrap_or(LOG10_2 / std::f64::consts::SQRT_2);
            if let Some(e) = build_edge(a, b, sa, sb, v, p) {
                edges.push(e);
            }
        }
    }
    // stable order for reproducible output
    edges.sort_by_key(|e| (e.a, e.b));
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
    pub delta: Vec<f64>, // per cohort (anchor = 0)
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

#[cfg(test)]
mod tests {
    use super::*;

    fn test_params() -> Params {
        Params {
            min_support: 2,
            max_drift_dilutions: 1.0,
            bootstrap: 200,
            seed: 42,
            sigma: vec![LOG10_2 / std::f64::consts::SQRT_2; 4],
            lambda: LOG10_2 * LOG10_2,
            robust: false,
            nu: 4.0,
            drift_scale: 1,
        }
    }

    // Deterministic edge: se=0 so the bootstrap median equals the WLS point estimate.
    fn edge(a: usize, b: usize, delta: f64) -> Edge {
        Edge {
            a,
            b,
            d_min: 0.0,
            tau: 0.0,
            n: 5,
            delta,
            se: 0.0,
            weight: 1.0,
            cand: vec![],
        }
    }

    #[test]
    fn fit_recovers_known_offsets() {
        // true offsets: labA=0, labB=log10 2 (2x), labC=log10 5 (5x)
        let (l2, l5) = (2f64.log10(), 5f64.log10());
        let edges = vec![edge(1, 0, l2), edge(2, 0, l5), edge(2, 1, l5 - l2)];
        let s = fit(&edges, 3, 0, &test_params()).unwrap();
        assert!((10f64.powf(s.delta[0]) - 1.0).abs() < 1e-6);
        assert!((10f64.powf(s.delta[1]) - 2.0).abs() < 1e-6);
        assert!((10f64.powf(s.delta[2]) - 5.0).abs() < 1e-6);
    }

    #[test]
    fn anchor_invariance_lemma1() {
        // Re-anchoring must leave all pairwise offset differences unchanged.
        let (l2, l5) = (2f64.log10(), 5f64.log10());
        let edges = vec![edge(1, 0, l2), edge(2, 0, l5), edge(2, 1, l5 - l2)];
        let p = test_params();
        let s0 = fit(&edges, 3, 0, &p).unwrap();
        let s2 = fit(&edges, 3, 2, &p).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let d0 = s0.delta[i] - s0.delta[j];
                let d2 = s2.delta[i] - s2.delta[j];
                assert!(
                    (d0 - d2).abs() < 1e-9,
                    "anchor variance at ({i},{j}): {d0} vs {d2}"
                );
            }
        }
    }

    #[test]
    fn disconnected_graph_is_singular() {
        // cohort 3 has no edge -> graph disconnected from anchor -> Err
        let edges = vec![edge(1, 0, 0.3), edge(2, 0, 0.6)];
        assert!(fit(&edges, 4, 0, &test_params()).is_err());
    }

    #[test]
    fn build_edges_tau_median_and_empirical_sigma() {
        let mk = |dist, signed| Pair {
            a: 0,
            b: 1,
            dist,
            signed,
            si: "x".into(),
            sj: "y".into(),
            va: 0.0,
            vb: -signed,
        };
        let pairs = vec![mk(0.0, 2f64.log10()), mk(0.0, 2f64.log10()), mk(2.0, 0.5)];
        let mut p = test_params();
        p.sigma = vec![0.2, 0.3];
        p.lambda = 0.09;
        p.min_support = 2;
        let edges = build_edges(&pairs, 2, &p);
        assert_eq!(edges.len(), 1);
        let e = &edges[0];
        assert_eq!((e.a, e.b), (0, 1));
        // at tau=0 there are 2 concordant pairs (>= min_support, zero drift) -> chosen
        assert_eq!(e.n, 2);
        assert!((e.delta - 2f64.log10()).abs() < 1e-9);
        // empirical-sigma SE = sqrt((sa^2 + sb^2 + lambda*tau)/n), tau = 0 here
        let expect = ((0.2f64 * 0.2 + 0.3 * 0.3) / 2.0).sqrt();
        assert!((e.se - expect).abs() < 1e-9, "se {} vs {}", e.se, expect);
    }
}
