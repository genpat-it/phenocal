//! Mathematical-trace report: emits a self-contained Markdown document that
//! walks through every step of the actual run with the real numbers and the
//! formula behind it -- inputs and per-cohort resolutions, per-edge twin tables
//! and the adaptive-tau selection, the edge statistics (Delta, SE, weight), the
//! weighted least-squares normal equations as real matrices, the solution, the
//! bootstrap credible intervals, and the data-driven cutoff.
//!
//! Markdown (not LaTeX/HTML) so it renders directly on GitHub -- including the
//! `$...$`/`$$...$$` math and the tables -- with no compilation step.
//! It is the didactic worked example, auto-filled with this run's data.

use crate::dashboard::Ctx;
use crate::linalg::median;
use std::collections::HashMap;
use std::f64::consts::LOG10_2;

const MAX_TWINS: usize = 15;
const MAX_ISO: usize = 15;

/// Render a name as inline code so underscores/specials never break Markdown.
fn code(s: &str) -> String {
    format!("`{}`", s.replace('`', "'"))
}

fn fold(d: f64) -> f64 {
    10f64.powf(d)
}

pub fn render(c: &Ctx) -> String {
    let mut s = String::with_capacity(1 << 16);
    s.push_str("# phenocal — mathematical trace of this run\n\n");
    s.push_str(&format!(
        "Genome-anchored calibration of **{} cohorts** ({} cross-cohort pairs, {} edges), anchor **{}**. \
         Every step below uses the *real* numbers of this run.\n\n",
        c.cohorts.len(),
        c.n_cross,
        c.edges.len(),
        code(c.anchor_name)
    ));
    inputs_section(&mut s, c);
    model_section(&mut s);
    edges_section(&mut s, c);
    solve_section(&mut s, c);
    bootstrap_section(&mut s, c);
    cutoff_section(&mut s, c);
    labels_section(&mut s, c);
    s
}

fn inputs_section(s: &mut String, c: &Ctx) {
    s.push_str("## 1. Inputs and per-cohort resolution\n\n");
    s.push_str(&format!(
        "Parameters: $n_{{\\min}}={}$, $\\kappa={}$ (drift tolerance), $\\lambda={:.4}$, noise model `{}`, bootstrap $B={}$.\n\n",
        c.min_support, c.max_drift, c.lambda, c.sigma_mode, c.bootstrap
    ));
    s.push_str("The measurement resolution $\\sigma_c$ is the median spacing of each cohort's MIC grid (in $\\log_{10}$); it recovers $\\log_{10}2\\approx0.301$ for two-fold cohorts.\n\n");
    s.push_str("| Cohort | $n$ | $\\sigma_c$ (log10) | $\\sigma_c$ (dilutions) |\n|---|---:|---:|---:|\n");
    for (i, name) in c.cohorts.iter().enumerate() {
        let sig = c.sigma.get(i).copied().unwrap_or(0.0);
        s.push_str(&format!(
            "| {}{} | {} | {:.3} | {:.2} |\n",
            code(name),
            if i == c.anchor { " (anchor)" } else { "" },
            c.counts.get(name).copied().unwrap_or(0),
            sig,
            sig / LOG10_2
        ));
    }
    s.push('\n');
}

fn model_section(s: &mut String) {
    s.push_str("## 2. The measurement model\n\n");
    s.push_str("Observed $y_{ic}=\\log_{10}\\mathrm{MIC}_{ic}$ decomposes as\n\n");
    s.push_str("$$ y_{ic}=\\mu_i+\\delta_c+\\varepsilon_{ic} $$\n\n");
    s.push_str("with $\\mu_i$ the unknown true biology of genome $i$, $\\delta_c$ the cohort/protocol offset (what we estimate), and $\\varepsilon_{ic}$ zero-median noise. For two near-clonal isolates ($\\mu_i\\approx\\mu_j$) in cohorts $a,b$ the biology cancels:\n\n");
    s.push_str("$$ y_{ia}-y_{jb}=\\underbrace{(\\mu_i-\\mu_j)}_{\\approx 0}+(\\delta_a-\\delta_b)+(\\varepsilon_{ia}-\\varepsilon_{jb}) $$\n\n");
    s.push_str("The **median over many twin pairs** of an edge kills the noise, leaving $\\Delta_{ab}\\approx\\delta_a-\\delta_b$.\n\n");
}

fn edges_section(s: &mut String, c: &Ctx) {
    s.push_str("## 3. Per-edge construction (with the real twins)\n\n");
    for e in c.edges {
        let (a, b) = (code(&c.cohorts[e.a]), code(&c.cohorts[e.b]));
        s.push_str(&format!("### {} — {}\n\n", a, b));

        let mut cand: Vec<&crate::calib::EdgePair> = e.cand.iter().collect();
        cand.sort_by(|x, y| x.dist.partial_cmp(&y.dist).unwrap());
        let shown = cand.len().min(MAX_TWINS);
        s.push_str("Twin pairs, log-ratio $s=\\log_{10}\\frac{\\mathrm{MIC}_a}{\\mathrm{MIC}_b}$ (nearest first):\n\n");
        s.push_str("| sample $a$ | sample $b$ | $d$ | $\\mathrm{MIC}_a$ | $\\mathrm{MIC}_b$ | $s$ | used |\n|---|---|---:|---:|---:|---:|:--:|\n");
        for p in cand.iter().take(shown) {
            s.push_str(&format!(
                "| {} | {} | {:.0} | {:.3} | {:.3} | {:.3} | {} |\n",
                code(&p.si),
                code(&p.sj),
                p.dist,
                fold(p.va),
                fold(p.vb),
                p.signed,
                if p.used { "✓" } else { "" }
            ));
        }
        s.push('\n');
        if cand.len() > shown {
            s.push_str(&format!("*({} further pairs not shown.)*\n\n", cand.len() - shown));
        }

        s.push_str("**Adaptive $\\tau$ selection** — widen $\\tau$ from $d_{\\min}$ to the smallest distance with $\\ge n_{\\min}$ pairs and drift within tolerance:\n\n");
        tau_scan(s, e);

        let (sa, sb) = (c.sigma.get(e.a).copied().unwrap_or(0.0), c.sigma.get(e.b).copied().unwrap_or(0.0));
        s.push_str(&format!(
            "$$ \\Delta_{{ab}}=\\operatorname{{median}}\\{{s:d\\le\\tau\\}}={:.4}\\ (\\text{{fold }}{:.2}\\times),\\quad \\tau={:.0},\\ n={}. $$\n\n",
            e.delta, fold(e.delta), e.tau, e.n
        ));
        s.push_str(&format!(
            "$$ SE_{{ab}}=\\sqrt{{\\frac{{\\sigma_a^2+\\sigma_b^2+\\lambda\\tau}}{{n}}}}=\\sqrt{{\\frac{{{:.3}^2+{:.3}^2+{:.4}\\cdot{:.0}}}{{{}}}}}={:.4},\\quad w_{{ab}}=\\frac{{1}}{{SE^2}}={:.2}. $$\n\n",
            sa, sb, c.lambda, e.tau, e.n, e.se, e.weight
        ));
    }
}

fn tau_scan(s: &mut String, e: &crate::calib::Edge) {
    let mut pairs: Vec<(f64, f64)> = e.cand.iter().map(|p| (p.dist, p.signed)).collect();
    pairs.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
    if pairs.is_empty() {
        return;
    }
    let d_min = pairs[0].0;
    let base: Vec<f64> = pairs.iter().filter(|x| x.0 == d_min).map(|x| x.1).collect();
    let base = median(&base);
    let mut taus: Vec<f64> = pairs.iter().map(|x| x.0).collect();
    taus.dedup();
    s.push_str("| $\\tau$ | $n(\\tau)$ | median$_{\\le\\tau}$ | drift (dil.) | |\n|---:|---:|---:|---:|:--|\n");
    for &t in taus.iter().take(8) {
        let vals: Vec<f64> = pairs.iter().filter(|x| x.0 <= t).map(|x| x.1).collect();
        let med = median(&vals);
        let drift = (med - base).abs() / LOG10_2;
        let chosen = (t - e.tau).abs() < 1e-9;
        s.push_str(&format!(
            "| {:.0} | {} | {:.3} | {:.2} | {} |\n",
            t,
            vals.len(),
            med,
            drift,
            if chosen { "← chosen" } else { "" }
        ));
    }
    s.push('\n');
    if taus.len() > 8 {
        s.push_str("*(scan truncated to 8 rows.)*\n\n");
    }
}

fn solve_section(s: &mut String, c: &Ctx) {
    s.push_str("## 4. Weighted least squares (the real matrices)\n\n");
    s.push_str("One row per edge ($+1$ at $a$, $-1$ at $b$): $A\\delta\\approx\\Delta$, weights $W=\\operatorname{diag}(w_{ab})$. The solution solves the normal equations\n\n");
    s.push_str("$$ A^{\\top}WA\\,\\hat\\delta = A^{\\top}W\\Delta, $$\n\n");
    s.push_str("where $A^{\\top}WA$ is the weighted graph Laplacian; the anchor row/column is removed (fixed to $0$).\n\n");

    let free: Vec<usize> = (0..c.cohorts.len()).filter(|&i| i != c.anchor).collect();
    let col: HashMap<usize, usize> = free.iter().enumerate().map(|(j, &i)| (i, j)).collect();
    let k = free.len();
    let mut ata = vec![0.0; k * k];
    let mut atb = vec![0.0; k];
    for e in c.edges {
        let w = e.weight;
        let ca = col.get(&e.a).copied();
        let cb = col.get(&e.b).copied();
        if let Some(i) = ca {
            ata[i * k + i] += w;
            atb[i] += w * e.delta;
            if let Some(j) = cb {
                ata[i * k + j] -= w;
            }
        }
        if let Some(j) = cb {
            ata[j * k + j] += w;
            atb[j] -= w * e.delta;
            if let Some(i) = ca {
                ata[j * k + i] -= w;
            }
        }
    }
    s.push_str("Free cohorts (order): ");
    s.push_str(&free.iter().map(|&i| code(&c.cohorts[i])).collect::<Vec<_>>().join(", "));
    s.push_str(".\n\n");
    s.push_str("$$ A^{\\top}WA=\\begin{bmatrix}\n");
    for i in 0..k {
        let row: Vec<String> = (0..k).map(|j| format!("{:.2}", ata[i * k + j])).collect();
        s.push_str(&row.join(" & "));
        s.push_str(" \\\\\n");
    }
    s.push_str("\\end{bmatrix},\\quad A^{\\top}W\\Delta=\\begin{bmatrix}\n");
    for v in atb.iter() {
        s.push_str(&format!("{:.2} \\\\\n", v));
    }
    s.push_str("\\end{bmatrix}. $$\n\n");
    s.push_str("Solving gives the offsets (anchor $=0$):\n\n");
    s.push_str("| Cohort | $\\hat\\delta$ (log10) | fold |\n|---|---:|---:|\n");
    for (i, name) in c.cohorts.iter().enumerate() {
        s.push_str(&format!(
            "| {}{} | {:.4} | {:.2}× |\n",
            code(name),
            if i == c.anchor { " (anchor)" } else { "" },
            c.sol.delta[i],
            fold(c.sol.delta[i])
        ));
    }
    s.push('\n');
    s.push_str(&format!(
        "Weighted residual RMSE $={:.4}$ log10 ($={:.2}$ dilutions). Offsets are identifiable only up to an additive constant (the anchor is a gauge choice; differences are unique).\n\n",
        c.sol.rmse,
        c.sol.rmse / LOG10_2
    ));
}

fn bootstrap_section(s: &mut String, c: &Ctx) {
    s.push_str("## 5. Uncertainty: bootstrap credible intervals\n\n");
    s.push_str(&format!(
        "Each edge is resampled $\\Delta^{{*}}\\sim\\mathcal{{N}}(\\Delta,SE^2)$ ({}) and the WLS re-solved $B={}$ times; the 2.5/97.5 percentiles give the 95% interval.\n\n",
        if c.robust { format!("robust Student-$t$, $\\nu={}$", c.nu) } else { "Normal".into() },
        c.bootstrap
    ));
    s.push_str("| Cohort | fold | 95% CrI |\n|---|---:|---:|\n");
    for (i, name) in c.cohorts.iter().enumerate() {
        s.push_str(&format!(
            "| {} | {:.2}× | [{:.2}, {:.2}] |\n",
            code(name),
            fold(c.sol.delta[i]),
            fold(c.sol.lo95[i]),
            fold(c.sol.hi95[i])
        ));
    }
    s.push('\n');
}

fn cutoff_section(s: &mut String, c: &Ctx) {
    let cut = c.cutoff;
    s.push_str("## 6. Data-driven cutoff (harmonised scale)\n\n");
    s.push_str("A 2-component Gaussian mixture is fitted to the harmonised $\\log_{10}$-MIC; the KDE antimode (density valley) and the mixture crossover (Bayes-optimal boundary) bracket the cutoff.\n\n");
    let mgl = |o: Option<f64>| o.map(|x| format!("{:.3}", fold(x))).unwrap_or_else(|| "n/a".into());
    let ci = |x: f64| if x.is_finite() { format!("{:.3}", fold(x)) } else { "n/a".into() };
    s.push_str(&format!(
        "Mixture components: sensitive ≈ {:.2} mg/L (weight {:.2}); tolerant ≈ {:.2} mg/L (weight {:.2}).\n\n",
        fold(cut.comp_lo.1), cut.comp_lo.0, fold(cut.comp_hi.1), cut.comp_hi.0
    ));
    s.push_str(&format!(
        "- **KDE antimode**: {} mg/L\n- **GMM crossover**: {} mg/L  (95% CI [{}, {}])\n\n",
        mgl(cut.kde_antimode),
        mgl(cut.gmm_crossover),
        ci(cut.gmm_lo),
        ci(cut.gmm_hi)
    ));
}

fn labels_section(s: &mut String, c: &Ctx) {
    s.push_str("## 7. Per-isolate harmonised values and probabilistic labels\n\n");
    s.push_str(&format!(
        "Harmonised MIC $=\\mathrm{{MIC}}/10^{{\\delta_c}}$; $P(T)$ is the bootstrap fraction clearing each threshold. First {} of {} isolates:\n\n",
        c.isolates.len().min(MAX_ISO),
        c.isolates.len()
    ));
    s.push_str("| sample | cohort | raw | harm.");
    for t in c.thresholds {
        s.push_str(&format!(" | $P_{{\\ge{:.2}}}$", t));
    }
    s.push_str(" |\n|---|---|---:|---:|");
    for _ in c.thresholds {
        s.push_str("---:|");
    }
    s.push('\n');
    for iso in c.isolates.iter().take(MAX_ISO) {
        s.push_str(&format!(
            "| {} | {} | {:.2} | {:.2}",
            code(&iso.sample),
            code(&c.cohorts[iso.cohort]),
            iso.raw,
            iso.harm
        ));
        for p in &iso.probs {
            s.push_str(&format!(" | {:.2}", p));
        }
        s.push_str(" |\n");
    }
    s.push('\n');
}
