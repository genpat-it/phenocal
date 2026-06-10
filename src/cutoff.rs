//! Data-driven sensitive/tolerant cutoff on the harmonised log10-MIC scale,
//! by two dependency-free methods. The KDE antimode is the density valley
//! between the two modes (non-parametric). The GMM crossover is the
//! Bayes-optimal boundary of a 2-component Gaussian mixture fitted by EM, with
//! a bootstrap 95% interval (parametric). The two bracket the cutoff; phenocal
//! does not commit to a single value.

use crate::linalg::quantile;
use crate::rng::Rng;

const SQRT_2PI: f64 = 2.5066282746310002;

fn npdf(x: f64, mu: f64, sd: f64) -> f64 {
    let z = (x - mu) / sd;
    (-0.5 * z * z).exp() / (sd * SQRT_2PI)
}

fn std_dev(v: &[f64], mean: f64) -> f64 {
    if v.len() < 2 {
        return 1.0;
    }
    let s: f64 = v.iter().map(|x| (x - mean) * (x - mean)).sum();
    (s / (v.len() as f64 - 1.0)).sqrt()
}

/// KDE antimode: interior density minimum between the two highest peaks.
fn kde_antimode(sorted: &[f64]) -> Option<f64> {
    if sorted.len() < 10 {
        return None;
    }
    let mean = sorted.iter().sum::<f64>() / sorted.len() as f64;
    let sd = std_dev(sorted, mean);
    let h = (1.06 * sd * (sorted.len() as f64).powf(-0.2)).max(1e-3); // Silverman
    let (lo, hi) = (sorted[0], sorted[sorted.len() - 1]);
    let m = 512usize;
    let dens: Vec<(f64, f64)> = (0..m)
        .map(|i| {
            let x = lo + (hi - lo) * i as f64 / (m as f64 - 1.0);
            let d = sorted.iter().map(|&v| npdf(x, v, h)).sum::<f64>() / sorted.len() as f64;
            (x, d)
        })
        .collect();
    // dominant mode (usually the sensitive peak)
    let gm = (0..m)
        .max_by(|&a, &b| dens[a].1.partial_cmp(&dens[b].1).unwrap())
        .unwrap();
    // first valley to the RIGHT of the dominant mode that is followed by a
    // secondary peak (the sensitive->tolerant antimode); else search left.
    let valley = |range: &mut dyn Iterator<Item = usize>| -> Option<usize> {
        let idx: Vec<usize> = range.collect();
        for w in 1..idx.len().saturating_sub(1) {
            let i = idx[w];
            let (pi, ni) = (idx[w - 1], idx[w + 1]);
            if dens[i].1 < dens[pi].1 && dens[i].1 <= dens[ni].1 {
                let later_max = idx[w + 1..].iter().map(|&j| dens[j].1).fold(0.0, f64::max);
                if later_max > dens[i].1 * 1.05 {
                    return Some(i);
                }
            }
        }
        None
    };
    if let Some(i) = valley(&mut (gm + 1..m)) {
        return Some(dens[i].0);
    }
    if let Some(i) = valley(&mut (0..gm).rev()) {
        return Some(dens[i].0);
    }
    None
}

/// 2-component 1D Gaussian mixture by EM. Returns (w,mu,sd) per comp, sorted by mu.
fn gmm2(v: &[f64]) -> ((f64, f64, f64), (f64, f64, f64)) {
    let n = v.len() as f64;
    let mut sorted = v.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mean = v.iter().sum::<f64>() / n;
    let var = std_dev(v, mean).powi(2).max(1e-6);
    let mut mu1 = quantile(&sorted, 0.25);
    let mut mu2 = quantile(&sorted, 0.75);
    let mut sd1 = var.sqrt();
    let mut sd2 = var.sqrt();
    let mut w1 = 0.5;
    for _ in 0..200 {
        let mut sr1 = 0.0;
        let mut s_mu1 = 0.0;
        let mut s_mu2 = 0.0;
        let mut resp: Vec<f64> = Vec::with_capacity(v.len());
        for &x in v {
            let a = w1 * npdf(x, mu1, sd1.max(1e-4));
            let b = (1.0 - w1) * npdf(x, mu2, sd2.max(1e-4));
            let r1 = if a + b > 0.0 { a / (a + b) } else { 0.5 };
            resp.push(r1);
            sr1 += r1;
            s_mu1 += r1 * x;
            s_mu2 += (1.0 - r1) * x;
        }
        let sr2 = n - sr1;
        let new_mu1 = if sr1 > 1e-9 { s_mu1 / sr1 } else { mu1 };
        let new_mu2 = if sr2 > 1e-9 { s_mu2 / sr2 } else { mu2 };
        let mut t1 = 0.0;
        let mut t2 = 0.0;
        for (k, &x) in v.iter().enumerate() {
            let r1 = resp[k];
            t1 += r1 * (x - new_mu1).powi(2);
            t2 += (1.0 - r1) * (x - new_mu2).powi(2);
        }
        mu1 = new_mu1;
        mu2 = new_mu2;
        sd1 = (t1 / sr1.max(1e-9)).max(1e-6).sqrt();
        sd2 = (t2 / sr2.max(1e-9)).max(1e-6).sqrt();
        w1 = sr1 / n;
    }
    let c1 = (w1, mu1, sd1);
    let c2 = (1.0 - w1, mu2, sd2);
    if mu1 <= mu2 {
        (c1, c2)
    } else {
        (c2, c1)
    }
}

/// Crossover x in (mu1,mu2) where w1*N1(x) = w2*N2(x): the Bayes boundary.
fn crossover(c1: (f64, f64, f64), c2: (f64, f64, f64)) -> Option<f64> {
    let ((w1, m1, s1), (w2, m2, s2)) = (c1, c2);
    if m1 >= m2 {
        return None;
    }
    let steps = 1000;
    let mut prev = w1 * npdf(m1, m1, s1) - w2 * npdf(m1, m2, s2);
    for i in 1..=steps {
        let x = m1 + (m2 - m1) * i as f64 / steps as f64;
        let cur = w1 * npdf(x, m1, s1) - w2 * npdf(x, m2, s2);
        if prev.signum() != cur.signum() {
            return Some(x);
        }
        prev = cur;
    }
    None
}

pub struct Cutoff {
    pub kde_antimode: Option<f64>, // log10
    pub gmm_crossover: Option<f64>,
    pub gmm_lo: f64,
    pub gmm_hi: f64,
    pub comp_lo: (f64, f64, f64), // (weight, mu, sd) of the sensitive component
    pub comp_hi: (f64, f64, f64), // (weight, mu, sd) of the tolerant component
}

impl Cutoff {
    /// Placeholder used when cutoff estimation is skipped (`--no-cutoff`).
    pub fn empty() -> Self {
        Cutoff {
            kde_antimode: None,
            gmm_crossover: None,
            gmm_lo: f64::NAN,
            gmm_hi: f64::NAN,
            comp_lo: (f64::NAN, f64::NAN, f64::NAN),
            comp_hi: (f64::NAN, f64::NAN, f64::NAN),
        }
    }
}

/// Estimate the cutoff from harmonised log10-MIC values.
pub fn estimate(logvals: &[f64], bootstrap: usize, seed: u64) -> Cutoff {
    let mut sorted = logvals.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let kde = kde_antimode(&sorted);
    let (c1, c2) = gmm2(logvals);
    let xc = crossover(c1, c2);
    let mut rng = Rng::new(seed ^ 0x00C0FFEE);
    let mut draws: Vec<f64> = Vec::new();
    let n = logvals.len();
    for _ in 0..bootstrap.max(1) {
        let bs: Vec<f64> = (0..n).map(|_| logvals[rng.next_index(n)]).collect();
        let (b1, b2) = gmm2(&bs);
        if let Some(x) = crossover(b1, b2) {
            draws.push(x);
        }
    }
    draws.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let (lo, hi) = if draws.len() >= 20 {
        (quantile(&draws, 0.025), quantile(&draws, 0.975))
    } else {
        (f64::NAN, f64::NAN)
    };
    Cutoff {
        kde_antimode: kde,
        gmm_crossover: xc,
        gmm_lo: lo,
        gmm_hi: hi,
        comp_lo: c1,
        comp_hi: c2,
    }
}

/// One weighted component density w·N(x|mu,sd) for plotting.
pub fn component_density(x: f64, c: (f64, f64, f64)) -> f64 {
    c.0 * npdf(x, c.1, c.2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crossover_lies_between_two_modes() {
        // bimodal log10 values: a tight cluster near 0 and another near 1
        let mut v = Vec::new();
        for i in 0..60 {
            v.push(0.0 + (i as f64) * 0.002); // ~[0, 0.12]
            v.push(1.0 + (i as f64) * 0.002); // ~[1, 1.12]
        }
        let c = estimate(&v, 50, 7);
        let x = c
            .gmm_crossover
            .expect("expected a GMM crossover for clearly bimodal data");
        assert!(
            x > 0.1 && x < 1.0,
            "crossover {x} not between the two modes"
        );
        // the two fitted components must straddle the crossover
        assert!(
            c.comp_lo.1 < x && c.comp_hi.1 > x,
            "components do not straddle crossover"
        );
    }
}
