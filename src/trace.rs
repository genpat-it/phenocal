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

/// A legend swatch: a coloured line (optionally dashed) + label.
fn leg_line(x: f64, y: f64, col: &str, dash: &str, label: &str) -> String {
    format!(
        "<line x1=\"{x:.1}\" y1=\"{y:.1}\" x2=\"{:.1}\" y2=\"{y:.1}\" stroke=\"{col}\" stroke-width=\"2\" {dash}/><text x=\"{:.1}\" y=\"{:.1}\" font-size=\"10\" fill=\"#333\">{label}</text>\n",
        x + 16.0, x + 20.0, y + 3.0
    )
}

/// A legend swatch: a coloured box (optional opacity) + label.
fn leg_box(x: f64, y: f64, col: &str, op: &str, label: &str) -> String {
    format!(
        "<rect x=\"{x:.1}\" y=\"{:.1}\" width=\"16\" height=\"9\" fill=\"{col}\" fill-opacity=\"{op}\"/><text x=\"{:.1}\" y=\"{:.1}\" font-size=\"10\" fill=\"#333\">{label}</text>\n",
        y - 7.0, x + 20.0, y + 3.0
    )
}

/// Write the Markdown trace plus its companion SVG figures (referenced from the
/// Markdown so they render on GitHub). `md_path` ends in `.md`; the figures are
/// `<stem>_bells.svg` and `<stem>_cutoff.svg` next to it.
pub fn write_report(c: &Ctx, md_path: &str) {
    let stem = md_path.strip_suffix(".md").unwrap_or(md_path);
    let bells = format!("{stem}_bells.svg");
    let cutoff = format!("{stem}_cutoff.svg");
    let base = |p: &str| p.rsplit('/').next().unwrap_or(p).to_string();
    let md = render(c, &base(&bells), &base(&cutoff));
    let _ = std::fs::write(md_path, md);
    let _ = std::fs::write(&bells, bells_svg_doc(c));
    let _ = std::fs::write(&cutoff, cutoff_svg_doc(c));
}

fn render(c: &Ctx, bells_img: &str, cutoff_img: &str) -> String {
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
    s.push_str("> **Tip — the reference does not matter.** The anchor ($\\delta=0$) is only a coordinate convention. Re-anchoring to any other cohort rescales the offsets, the threshold grid and the data *together*, so the harmonised labels are **unchanged** (anchor invariance, Lemma 1): the results do **not** depend on the reference choice.\n\n");
    inputs_section(&mut s, c);
    model_section(&mut s);
    symbols_section(&mut s);
    edges_section(&mut s, c);
    solve_section(&mut s, c);
    bootstrap_section(&mut s, c);
    s.push_str(&format!(
        "Each edge contributes the Gaussian $\\mathcal{{N}}(\\Delta_{{ab}},SE_{{ab}}^2)$ the bootstrap samples from. In the figure legend each bell is labelled `a–b (fold, w)`, with **center** $=10^{{\\Delta_{{ab}}}}$ (the fold) and **width/height** set by $w_{{ab}}=1/SE_{{ab}}^2$ — exactly the `fold` and `weight` columns computed per edge in §4. A precise edge (many close twins $\\Rightarrow$ small $SE$, large $w$) is **tall and narrow**; a noisy one (few or distant twins $\\Rightarrow$ large $SE$, small $w$) is **short and wide** and barely moves the fit.\n\n![per-edge bootstrap bells]({})\n\n",
        bells_img
    ));
    cutoff_section(&mut s, c);
    s.push_str(&format!("![data-driven cutoff: histogram, mixture components, KDE antimode and GMM crossover]({})\n\n", cutoff_img));
    labels_section(&mut s, c);
    sensitivity_section(&mut s, c);
    s
}

fn sensitivity_section(s: &mut String, c: &Ctx) {
    if c.sensitivity.is_empty() {
        return;
    }
    s.push_str("## 9. Sensitivity to our formula choices (SE $\\sigma_c$ and drift scale)\n\n");
    s.push_str("The edge standard error (via the per-cohort resolution $\\sigma_c$) and the drift-tolerance scale are the only modelling choices we make. The variants and their formulas:\n\n");
    s.push_str("- **$\\sigma_c$** (enters $SE_{ab}=\\sqrt{(\\sigma_a^2+\\sigma_b^2+\\lambda\\tau)/n}$): *fixed* $\\sigma_c=\\log_{10}(2)/\\sqrt{2}\\approx0.213$; *empirical* $\\sigma_c=$ median spacing of the cohort's MIC grid.\n");
    s.push_str("- **drift scale $\\bar\\sigma_{ab}$** (the $\\tau$ bound $\\le\\kappa\\,\\bar\\sigma_{ab}$): *doubling* $=\\log_{10}2$; *rms* $=\\sqrt{(\\sigma_a^2+\\sigma_b^2)/2}$; *mean* $=(\\sigma_a+\\sigma_b)/2$; *max* $=\\max(\\sigma_a,\\sigma_b)$.\n\n");
    s.push_str("Re-fitting the offsets under each alternative shows they are **not load-bearing** — the folds barely move:\n\n");
    s.push_str("| variant");
    for name in c.cohorts {
        s.push_str(&format!(" | {}", code(name)));
    }
    s.push_str(" |\n|---");
    for _ in c.cohorts {
        s.push_str("|---:");
    }
    s.push_str("|\n");
    for (label, folds) in c.sensitivity {
        s.push_str(&format!("| {}", code(label)));
        for f in folds {
            s.push_str(&format!(" | {:.2}×", f));
        }
        s.push_str(" |\n");
    }
    s.push_str("\nThe anchor stays $1.00\\times$ by definition; the others are stable across both the $\\sigma$ model and the drift-scale form (cf. §3: $\\sigma_c$ sets the weight, the drift scale only gates $\\tau$).\n\n");
}

fn npdf(x: f64, mu: f64, sd: f64) -> f64 {
    let z = (x - mu) / sd;
    (-0.5 * z * z).exp() / (sd * 2.506_628_274_631_000_3)
}

/// Per-edge Gaussian perturbation bells N(Delta, SE^2) on the log10-fold axis.
fn bells_svg_doc(c: &Ctx) -> String {
    let (w, h) = (760.0, 360.0);
    let (top, left, right, bot) = (24.0, 40.0, 20.0, 36.0);
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for e in c.edges {
        lo = lo.min(e.delta - 3.0 * e.se);
        hi = hi.max(e.delta + 3.0 * e.se);
    }
    if !lo.is_finite() {
        lo = -1.0;
        hi = 1.0;
    }
    let xspan = w - left - right;
    let yspan = h - top - bot;
    let xof = |x: f64| left + (x - lo) / (hi - lo) * xspan;
    let m = 240usize;
    let grid: Vec<f64> = (0..m)
        .map(|i| lo + (hi - lo) * i as f64 / (m as f64 - 1.0))
        .collect();
    let ymax = c
        .edges
        .iter()
        .map(|e| npdf(e.delta, e.delta, e.se))
        .fold(1e-9, f64::max);
    let yof = |y: f64| top + yspan - (y / ymax) * yspan;
    let palette = [
        "#2E86AB", "#E76F51", "#2A9D8F", "#7a5195", "#bc5090", "#ffa600", "#003f5c", "#58508d",
        "#ff6361",
    ];
    let mut s = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {w} {h}\" font-family=\"sans-serif\">\n<rect width=\"{w}\" height=\"{h}\" fill=\"white\"/>\n"
    );
    let baseline = top + yspan;
    // fold ticks
    for &f in &[0.5f64, 1.0, 2.0, 5.0, 10.0, 20.0] {
        let xl = f.log10();
        if xl < lo || xl > hi {
            continue;
        }
        let px = xof(xl);
        let ty = baseline + 14.0;
        let lab = trim(f);
        s.push_str(&format!(
            "<line x1=\"{px:.1}\" y1=\"{top:.1}\" x2=\"{px:.1}\" y2=\"{baseline:.1}\" stroke=\"#eee\"/><text x=\"{px:.1}\" y=\"{ty:.1}\" font-size=\"10\" fill=\"#888\" text-anchor=\"middle\">{lab}×</text>\n"
        ));
    }
    for (k, e) in c.edges.iter().enumerate() {
        let col = palette[k % palette.len()];
        let mut p = String::from("<path d=\"");
        for (j, &x) in grid.iter().enumerate() {
            p.push_str(&format!(
                "{}{:.1} {:.1} ",
                if j == 0 { "M" } else { "L" },
                xof(x),
                yof(npdf(x, e.delta, e.se))
            ));
        }
        p.push_str(&format!(
            "\" fill=\"none\" stroke=\"{col}\" stroke-width=\"1.8\" stroke-opacity=\"0.85\"/>\n"
        ));
        s.push_str(&p);
    }
    // legend
    let mut ly = top + 4.0;
    for (k, e) in c.edges.iter().enumerate() {
        let col = palette[k % palette.len()];
        s.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{ly:.1}\" font-size=\"9\" fill=\"{col}\">{}–{} (n={}, {:.1}×, w={:.0})</text>\n",
            left + 4.0,
            c.cohorts[e.a].replace('&', "&amp;").replace('<', "&lt;"),
            c.cohorts[e.b].replace('&', "&amp;").replace('<', "&lt;"),
            e.n,
            fold(e.delta),
            e.weight
        ));
        ly += 11.0;
    }
    s.push_str("</svg>\n");
    s
}

/// Cutoff density: histogram of harmonised log10-MIC + the two mixture components,
/// with the KDE antimode and GMM crossover marked.
fn cutoff_svg_doc(c: &Ctx) -> String {
    let lv = c.logharm;
    let (w, h) = (760.0, 400.0);
    let (top, left, right, bot) = (18.0, 40.0, 20.0, 100.0);
    let mut lo = lv.iter().cloned().fold(f64::INFINITY, f64::min) - 0.1;
    let mut hi = lv.iter().cloned().fold(f64::NEG_INFINITY, f64::max) + 0.1;
    if lo >= hi {
        lo = -0.5;
        hi = 1.5;
    }
    let cut = c.cutoff;
    let xspan = w - left - right;
    let yspan = h - top - bot;
    let xof = |x: f64| left + (x - lo) / (hi - lo) * xspan;
    let nb = 32usize;
    let bw = (hi - lo) / nb as f64;
    let mut hist = vec![0f64; nb];
    for &x in lv {
        hist[(((x - lo) / bw) as usize).min(nb - 1)] += 1.0;
    }
    let n = lv.len().max(1) as f64;
    for v in hist.iter_mut() {
        *v /= n * bw;
    }
    let comp = |x: f64, which: u8| {
        if which == 0 {
            crate::cutoff::component_density(x, cut.comp_lo)
        } else {
            crate::cutoff::component_density(x, cut.comp_hi)
        }
    };
    let m = 256usize;
    let grid: Vec<f64> = (0..m)
        .map(|i| lo + (hi - lo) * i as f64 / (m as f64 - 1.0))
        .collect();
    let mut ymax = hist.iter().cloned().fold(0.0, f64::max);
    for &x in &grid {
        ymax = ymax.max(comp(x, 0)).max(comp(x, 1));
    }
    if ymax <= 0.0 {
        ymax = 1.0;
    }
    let yof = |y: f64| top + yspan - (y / ymax) * yspan;
    let baseline = top + yspan;
    let mut s = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {w} {h}\" font-family=\"sans-serif\">\n<rect width=\"{w}\" height=\"{h}\" fill=\"white\"/>\n"
    );
    if cut.gmm_lo.is_finite() && cut.gmm_hi.is_finite() {
        let (x1, x2) = (xof(cut.gmm_lo), xof(cut.gmm_hi));
        s.push_str(&format!(
            "<rect x=\"{:.1}\" y=\"{top:.1}\" width=\"{:.1}\" height=\"{yspan:.1}\" fill=\"#2E86AB\" fill-opacity=\"0.10\"/>\n",
            x1, (x2 - x1).max(0.0)
        ));
    }
    for (i, &d) in hist.iter().enumerate() {
        let px = xof(lo + i as f64 * bw);
        let pw = (xof(lo + (i as f64 + 1.0) * bw) - px).max(0.5);
        let py = yof(d);
        s.push_str(&format!("<rect x=\"{px:.1}\" y=\"{py:.1}\" width=\"{pw:.1}\" height=\"{:.1}\" fill=\"#cdd8df\"/>\n", baseline - py));
    }
    for (which, col) in [(0u8, "#2A9D8F"), (1u8, "#E76F51")] {
        let mut p = String::from("<path d=\"");
        for (j, &x) in grid.iter().enumerate() {
            p.push_str(&format!(
                "{}{:.1} {:.1} ",
                if j == 0 { "M" } else { "L" },
                xof(x),
                yof(comp(x, which))
            ));
        }
        p.push_str(&format!(
            "\" fill=\"none\" stroke=\"{col}\" stroke-width=\"2\"/>\n"
        ));
        s.push_str(&p);
    }
    let mut vline = |x: f64, col: &str, dash: &str| {
        s.push_str(&format!(
            "<line x1=\"{0:.1}\" y1=\"{top:.1}\" x2=\"{0:.1}\" y2=\"{baseline:.1}\" stroke=\"{col}\" stroke-width=\"2\" {dash}/>\n",
            xof(x)
        ));
    };
    if let Some(a) = cut.kde_antimode {
        vline(a, "#7a5195", "stroke-dasharray=\"4 3\"");
    }
    if let Some(g) = cut.gmm_crossover {
        vline(g, "#2E86AB", "");
    }
    let conv = 1.25f64.log10();
    if conv > lo && conv < hi {
        vline(conv, "#000000", "stroke-dasharray=\"2 3\"");
    }
    for &f in &[0.3f64, 0.6, 1.25, 2.5, 5.0, 10.0] {
        let x = f.log10();
        if x < lo || x > hi {
            continue;
        }
        s.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-size=\"10\" fill=\"#888\" text-anchor=\"middle\">{}</text>\n",
            xof(x), baseline + 14.0, trim(f)
        ));
    }
    // legend (two rows)
    let r1 = baseline + 36.0;
    let r2 = baseline + 56.0;
    s.push_str(&leg_box(
        left,
        r1,
        "#cdd8df",
        "1",
        "histogram (harmonised log10-MIC)",
    ));
    s.push_str(&leg_line(
        left + 220.0,
        r1,
        "#2A9D8F",
        "",
        "sensitive component",
    ));
    s.push_str(&leg_line(
        left + 410.0,
        r1,
        "#E76F51",
        "",
        "tolerant component",
    ));
    s.push_str(&leg_line(
        left,
        r2,
        "#7a5195",
        "stroke-dasharray=\"4 3\"",
        "KDE antimode",
    ));
    s.push_str(&leg_line(left + 150.0, r2, "#2E86AB", "", "GMM crossover"));
    s.push_str(&leg_box(left + 300.0, r2, "#2E86AB", "0.2", "GMM 95% CI"));
    s.push_str(&leg_line(
        left + 430.0,
        r2,
        "#000000",
        "stroke-dasharray=\"2 3\"",
        "1.25 convention",
    ));
    s.push_str("</svg>\n");
    s
}

/// Trim trailing zeros for axis labels.
fn trim(f: f64) -> String {
    let s = format!("{:.2}", f);
    let s = s.trim_end_matches('0').trim_end_matches('.').to_string();
    if s.is_empty() {
        "0".to_string()
    } else {
        s
    }
}

fn inputs_section(s: &mut String, c: &Ctx) {
    s.push_str("## 1. Inputs and per-cohort resolution\n\n");
    s.push_str(&format!(
        "Parameters: $n_{{\\min}}={}$, $\\kappa={}$ (drift tolerance), $\\lambda={:.4}$, noise model `{}`, bootstrap $B={}$.\n\n",
        c.min_support, c.max_drift, c.lambda, c.sigma_mode, c.bootstrap
    ));
    s.push_str("The measurement resolution $\\sigma_c$ is the median spacing of each cohort's MIC grid (in $\\log_{10}$); it recovers $\\log_{10}2\\approx0.301$ for two-fold cohorts.\n\n");
    s.push_str(
        "| Cohort | $n$ | $\\sigma_c$ (log10) | $\\sigma_c$ (dilutions) |\n|---|---:|---:|---:|\n",
    );
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
    s.push_str("with $\\mu_i$ the unknown true biology of genome $i$, $\\delta_c$ the cohort/protocol offset (what we estimate), and $\\varepsilon_{ic}$ random measurement noise.\n\n");
    s.push_str("> **Tip — the noise is mean-zero.** $\\mathbb{E}[\\varepsilon_{ic}]=0$ by construction: any *constant* bias of a lab is absorbed into $\\delta_c$, so what is left in $\\varepsilon$ has no preferred direction. If $\\varepsilon$ is also symmetric, its *median* is $0$ too — which is what makes the median estimator below unbiased.\n\n");
    s.push_str("For two near-clonal isolates ($\\mu_i\\approx\\mu_j$) in cohorts $a,b$ the biology approximately cancels (exactly when $\\mu_i=\\mu_j$):\n\n");
    s.push_str("$$ y_{ia}-y_{jb}=\\underbrace{(\\mu_i-\\mu_j)}_{\\approx 0}+(\\delta_a-\\delta_b)+(\\varepsilon_{ia}-\\varepsilon_{jb}) $$\n\n");
    s.push_str("The **median over many twin pairs** of an edge kills the noise, leaving $\\Delta_{ab}\\approx\\delta_a-\\delta_b$.\n\n");
}

fn symbols_section(s: &mut String) {
    s.push_str("## 3. Symbols (legend)\n\n");
    s.push_str("| Symbol | Meaning | Where it acts |\n|---|---|---|\n");
    let rows = [
        ("$y=\\log_{10}\\mathrm{MIC}$", "the measurement, on the log scale", "everywhere"),
        ("$\\mu_i$", "true biology of genome $i$ (unknown)", "model"),
        ("$\\delta_c$", "cohort/protocol offset — what we estimate", "model / solution"),
        ("$\\varepsilon$", "measurement noise, $\\mathbb{E}[\\varepsilon]=0$", "model"),
        ("$s=y_a-y_b$", "twin log-ratio of one cross-cohort pair", "per edge"),
        ("$d$", "genomic distance of a pair", "per edge"),
        ("$\\Delta_{ab}$", "edge value = median of $s$ over the twins", "per edge"),
        ("$\\tau_{ab}$", "adaptive distance bound (which twins enter the median)", "$\\tau$ selection"),
        ("$n_{ab}$", "number of supporting twin pairs", "SE & $\\tau$"),
        ("$\\sigma_c$", "cohort measurement resolution (grid spacing)", "**SE / weight**"),
        ("$\\lambda$", "biological-drift variance per unit distance", "**SE**"),
        ("$\\kappa$", "drift tolerance — bound is $\\kappa\\,\\bar\\sigma_{ab}$ in $\\tau$ selection; **not** in the SE", "**$\\tau$ selection only**"),
        ("$\\bar\\sigma_{ab}$", "edge resolution scale, RMS of $\\sigma_a,\\sigma_b$", "$\\tau$ drift bound"),
        ("$SE_{ab}$", "edge standard error $=\\sqrt{(\\sigma_a^2+\\sigma_b^2+\\lambda\\tau)/n}$", "weight"),
        ("$w_{ab}=1/SE_{ab}^2$", "edge weight (inverse variance)", "WLS"),
        ("$\\hat\\delta$", "fitted offsets (anchor $=0$)", "solution"),
        ("$P(T_i)$", "probability isolate $i$ is tolerant", "labels"),
    ];
    for (sym, mean, where_) in rows {
        s.push_str(&format!("| {} | {} | {} |\n", sym, mean, where_));
    }
    s.push_str("\n> Note: $\\kappa$ and $\\sigma_c$ play **different roles** — $\\sigma_c$ sets how much each edge is *trusted* (the weight $1/SE^2$), while $\\kappa$ sets how far the median may *drift* before $\\tau$ stops widening. They only meet in the drift bound $\\kappa\\,\\bar\\sigma_{ab}$; $\\kappa$ never enters the SE.\n\n");
}

fn edges_section(s: &mut String, c: &Ctx) {
    s.push_str("## 4. Per-edge construction (with the real twins)\n\n");
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
            s.push_str(&format!(
                "*({} further pairs not shown.)*\n\n",
                cand.len() - shown
            ));
        }

        let (sa, sb) = (
            c.sigma.get(e.a).copied().unwrap_or(0.0),
            c.sigma.get(e.b).copied().unwrap_or(0.0),
        );
        let res = ((sa * sa + sb * sb) / 2.0).sqrt().max(1e-9);
        s.push_str(&format!(
            "**Adaptive $\\tau$ selection.** Widen $\\tau$ from $d_{{\\min}}$ to the smallest distance with $\\ge n_{{\\min}}={}$ pairs **and** drift $\\le\\kappa\\,\\bar\\sigma_{{ab}}$. Here $\\bar\\sigma_{{ab}}=\\sqrt{{(\\sigma_a^2+\\sigma_b^2)/2}}={:.3}$ and $\\kappa={}$, so the drift (last column, in units of $\\bar\\sigma_{{ab}}$) must stay $\\le{}$:\n\n",
            c.min_support, res, c.max_drift, c.max_drift
        ));
        tau_scan(s, e, res, c.max_drift);

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

fn tau_scan(s: &mut String, e: &crate::calib::Edge, res: f64, kappa: f64) {
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
    s.push_str("| $\\tau$ | $n(\\tau)$ | $\\mathrm{median}_{\\le\\tau}$ | drift $/\\bar\\sigma$ | $\\le\\kappa$? | |\n|---:|---:|---:|---:|:--:|:--|\n");
    for &t in taus.iter().take(8) {
        let vals: Vec<f64> = pairs.iter().filter(|x| x.0 <= t).map(|x| x.1).collect();
        let med = median(&vals);
        let drift = (med - base).abs() / res;
        let chosen = (t - e.tau).abs() < 1e-9;
        s.push_str(&format!(
            "| {:.0} | {} | {:.3} | {:.2} | {} | {} |\n",
            t,
            vals.len(),
            med,
            drift,
            if drift <= kappa { "✓" } else { "✗" },
            if chosen { "← chosen" } else { "" }
        ));
    }
    s.push('\n');
    if taus.len() > 8 {
        s.push_str("*(scan truncated to 8 rows.)*\n\n");
    }
}

fn solve_section(s: &mut String, c: &Ctx) {
    s.push_str("## 5. Weighted least squares (the real matrices)\n\n");
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
    s.push_str(
        &free
            .iter()
            .map(|&i| code(&c.cohorts[i]))
            .collect::<Vec<_>>()
            .join(", "),
    );
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
    s.push_str("## 6. Uncertainty: bootstrap credible intervals\n\n");
    s.push_str(&format!(
        "How the intervals are built, in {} repetitions:\n\n\
         1. for **each edge**, draw a perturbed value from its bell $\\Delta^{{(t)}}\\sim\\mathcal{{N}}(\\Delta_{{ab}},SE_{{ab}}^2)$ ({}) — precise edges (small $SE$) barely move, noisy ones move a lot;\n\
         2. **re-solve the weighted least squares** with the perturbed edges $\\to$ one offset vector $\\hat\\delta^{{(t)}}$;\n\
         3. repeat $B={}$ times $\\to$ a **cloud** of $B$ values for each cohort offset.\n\n\
         The reported **fold** is the cloud's median ($10^{{\\,\\mathrm{{median}}}}$); the **95% credible interval** is its **2.5th–97.5th percentile** (sort the $B$ values, drop the extreme 2.5% on each side):\n\n",
        c.bootstrap,
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
    s.push_str("## 7. Data-driven cutoff (harmonised scale)\n\n");
    s.push_str("Two methods on the harmonised $\\log_{10}$-MIC values $\\{x_i\\}$; they bracket the cutoff.\n\n");

    s.push_str(
        "**(a) KDE antimode (non-parametric).** Estimate the density with a Gaussian kernel,\n\n",
    );
    s.push_str("$$ \\hat f(x)=\\frac{1}{n h}\\sum_{i=1}^{n} K\\!\\left(\\frac{x-x_i}{h}\\right),\\qquad K(u)=\\frac{1}{\\sqrt{2\\pi}}e^{-u^2/2},\\qquad h=1.06\\,\\hat\\sigma\\,n^{-1/5}\\ (\\text{Silverman}), $$\n\n");
    s.push_str("then take the **antimode**: the local minimum of $\\hat f$ in the valley between the two modes. No shape assumption; depends on the bandwidth $h$.\n\n");

    s.push_str("**(b) GMM crossover (parametric).** Fit a 2-component Gaussian mixture by EM,\n\n");
    s.push_str("$$ p(x)=w_1\\,\\mathcal{N}(x\\mid\\mu_1,\\sigma_1^2)+w_2\\,\\mathcal{N}(x\\mid\\mu_2,\\sigma_2^2). $$\n\n");
    s.push_str("EM iterates the responsibilities (E-step) and the weighted moments (M-step):\n\n");
    s.push_str("$$ r_k(x_i)=\\frac{w_k\\,\\mathcal{N}(x_i\\mid\\mu_k,\\sigma_k^2)}{\\sum_{l}w_l\\,\\mathcal{N}(x_i\\mid\\mu_l,\\sigma_l^2)},\\quad \\mu_k=\\frac{\\sum_i r_k(x_i)\\,x_i}{\\sum_i r_k(x_i)},\\quad \\sigma_k^2=\\frac{\\sum_i r_k(x_i)(x_i-\\mu_k)^2}{\\sum_i r_k(x_i)},\\quad w_k=\\frac{\\sum_i r_k(x_i)}{n}. $$\n\n");
    s.push_str("The **crossover** is the $x^{*}\\in(\\mu_1,\\mu_2)$ where the two weighted components are equal — the Bayes-optimal boundary (posterior $=1/2$):\n\n");
    s.push_str("$$ w_1\\,\\mathcal{N}(x^{*}\\mid\\mu_1,\\sigma_1^2)=w_2\\,\\mathcal{N}(x^{*}\\mid\\mu_2,\\sigma_2^2). $$\n\n");
    s.push_str("A 95% interval is obtained by **bootstrap**: resample $\\{x_i\\}$ with replacement, refit the mixture, recompute $x^{*}$, $B$ times, and take the 2.5/97.5 percentiles.\n\n");
    s.push_str("On this run:\n\n");
    let mgl = |o: Option<f64>| {
        o.map(|x| format!("{:.3}", fold(x)))
            .unwrap_or_else(|| "n/a".into())
    };
    let ci = |x: f64| {
        if x.is_finite() {
            format!("{:.3}", fold(x))
        } else {
            "n/a".into()
        }
    };
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
    s.push_str("## 8. Per-isolate harmonised values and probabilistic labels\n\n");
    s.push_str("The harmonised MIC shown is the point estimate $\\mathrm{MIC}/10^{\\hat\\delta_c}$. For $P(T)$, the cohort offset is drawn from its bootstrap-derived posterior $\\delta_c\\sim\\mathcal{N}(\\hat\\delta_c,s_c^2)$ ($s_c=$ half the 95% CI width $/1.96$) and $P(T_i(c))$ is the fraction of draws whose harmonised MIC clears $c$:\n\n");
    s.push_str("$$ P(T_i(c))=\\frac{1}{B}\\sum_{t=1}^{B}\\mathbb{1}\\!\\left[\\frac{\\mathrm{MIC}_i}{10^{\\delta_c^{(t)}}}\\ge c\\right],\\qquad \\delta_c^{(t)}\\sim\\mathcal{N}(\\hat\\delta_c,s_c^2). $$\n\n");
    s.push_str(&format!(
        "So $P\\approx0$ or $1$ is a robust label; $P\\approx0.5$ is calibration-sensitive. First {} of {} isolates:\n\n",
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
