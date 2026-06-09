//! Self-contained interactive HTML report for a phenocal run.
//!
//! Emits a single .html file (no external dependencies, no network): the cohort
//! comparison graph drawn in SVG, every twin pair and its log-ratio per edge
//! with the median highlighted, the fitted offsets with bootstrap credible
//! intervals (forest plot), before/after harmonised distributions, and the
//! per-isolate probabilistic labels with a live filter. Open it in any browser.

use crate::calib::{Edge, Solution};
use crate::linalg::{median, quantile};
use std::collections::HashMap;
use std::f64::consts::{LOG10_2, PI};

/// Per-isolate row passed in from main.
pub struct Iso {
    pub sample: String,
    pub cohort: usize,
    pub raw: f64,
    pub harm: f64,
    pub probs: Vec<f64>,
}

/// Everything the report needs.
pub struct Ctx<'a> {
    pub cohorts: &'a [String],
    pub counts: &'a HashMap<String, usize>,
    pub anchor: usize,
    pub anchor_name: &'a str,
    pub edges: &'a [Edge],
    pub sol: &'a Solution,
    pub sigma: &'a [f64],
    pub sigma_mode: &'a str,
    pub lambda: f64,
    pub robust: bool,
    pub nu: f64,
    pub min_support: usize,
    pub max_drift: f64,
    pub bootstrap: usize,
    pub n_cross: usize,
    pub thresholds: &'a [f64],
    pub isolates: &'a [Iso],
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn fold(delta: f64) -> f64 {
    10f64.powf(delta)
}

pub fn render(c: &Ctx) -> String {
    let mut h = String::with_capacity(1 << 18);
    h.push_str("<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">\n");
    h.push_str("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n");
    h.push_str("<title>phenocal report</title>\n");
    h.push_str(STYLE);
    h.push_str("</head>\n<body>\n");

    // ---- header ----
    h.push_str("<header><h1>phenocal &mdash; calibration report</h1>");
    h.push_str("<p class=\"sub\">Genome-anchored phenotype calibration. Everything the run found, inspectable.</p></header>\n");

    // ---- summary cards ----
    h.push_str("<section class=\"cards\">\n");
    card(&mut h, "Cohorts", &format!("{}", c.cohorts.len()));
    card(&mut h, "Anchor (&delta;=0)", &esc(c.anchor_name));
    card(&mut h, "Cross-cohort pairs", &format!("{}", c.n_cross));
    card(&mut h, "Edges", &format!("{}", c.edges.len()));
    card(&mut h, "Residual RMSE", &format!("{:.3} log<sub>10</sub><br><span class=\"dim\">{:.2} dilutions</span>", c.sol.rmse, c.sol.rmse / LOG10_2));
    let noise = if c.robust { format!("Student-t(&nu;={})", c.nu) } else { "Normal".to_string() };
    card(&mut h, "Bootstrap / noise", &format!("B={}<br><span class=\"dim\">{}</span>", c.bootstrap, noise));
    h.push_str("</section>\n");

    h.push_str(&format!(
        "<p class=\"params\">Parameters: n<sub>min</sub>={}, max drift={} dilutions, &sigma; mode=<b>{}</b>, &lambda;={:.4}. \
         Edge SE = &radic;((&sigma;<sub>a</sub><sup>2</sup>+&sigma;<sub>b</sub><sup>2</sup>+&lambda;&middot;&tau;)/n).</p>\n",
        c.min_support, c.max_drift, esc(c.sigma_mode), c.lambda
    ));

    // ---- comparison graph ----
    h.push_str("<h2>Comparison graph</h2>\n");
    h.push_str("<p class=\"note\">Nodes = cohorts (anchor in green). Each edge = a fitted cohort pair; thickness &prop; weight (precision), label = fold &amp; number of supporting twins.</p>\n");
    graph_svg(&mut h, c);

    // ---- forest plot of offsets ----
    h.push_str("<h2>Cohort offsets (fold vs anchor, 95% CrI)</h2>\n");
    forest_svg(&mut h, c);

    // offsets table
    h.push_str("<table class=\"tbl\"><thead><tr><th>cohort</th><th>n</th><th>&sigma;<sub>c</sub></th><th>&delta; (log<sub>10</sub>)</th><th>fold</th><th>95% CrI</th></tr></thead><tbody>\n");
    for (i, name) in c.cohorts.iter().enumerate() {
        let anchor = i == c.anchor;
        h.push_str(&format!(
            "<tr{}><td>{}{}</td><td>{}</td><td>{:.3}</td><td>{:.3}</td><td><b>{:.2}&times;</b></td><td>[{:.2}, {:.2}]</td></tr>\n",
            if anchor { " class=\"anchor\"" } else { "" },
            esc(name),
            if anchor { " &#9875;" } else { "" },
            c.counts.get(name).copied().unwrap_or(0),
            c.sigma.get(i).copied().unwrap_or(0.0),
            c.sol.delta[i],
            fold(c.sol.delta[i]),
            fold(c.sol.lo95[i]),
            fold(c.sol.hi95[i]),
        ));
    }
    h.push_str("</tbody></table>\n");

    // ---- before / after ----
    h.push_str("<h2>Before &rarr; after harmonisation</h2>\n");
    h.push_str("<p class=\"note\">Per cohort: median and inter-quartile range on the log<sub>10</sub> scale. Raw (hollow) vs harmonised (solid). Harmonisation should align the cohorts.</p>\n");
    beforeafter_svg(&mut h, c);

    // ---- per-edge detail (the miraculous debug) ----
    h.push_str("<h2>Per-edge detail &mdash; every twin and its log-ratio</h2>\n");
    h.push_str("<p class=\"note\">Click an edge to expand. Rows that entered the median (within the adaptive &tau;) are highlighted; greyed rows are candidates beyond &tau;. The edge &Delta; is the median of the highlighted log-ratios.</p>\n");
    for e in c.edges {
        edge_block(&mut h, c, e);
    }

    // ---- per-isolate labels ----
    h.push_str("<h2>Per-isolate harmonised values &amp; probabilities</h2>\n");
    h.push_str("<p class=\"note\">P(tolerant) = fraction of bootstrap calibrations in which the harmonised value clears each threshold. Type to filter by sample or cohort.</p>\n");
    h.push_str("<input id=\"flt\" placeholder=\"filter sample or cohort\u{2026}\" oninput=\"filt()\">\n");
    h.push_str("<div class=\"scroll\"><table class=\"tbl\" id=\"isotbl\"><thead><tr><th>sample</th><th>cohort</th><th>raw</th><th>harmonised</th>");
    for t in c.thresholds {
        h.push_str(&format!("<th>P&ge;{}</th>", trim(*t)));
    }
    h.push_str("</tr></thead><tbody>\n");
    for iso in c.isolates {
        let cohort = &c.cohorts[iso.cohort];
        h.push_str(&format!(
            "<tr><td>{}</td><td>{}</td><td>{:.3}</td><td>{:.3}</td>",
            esc(&iso.sample), esc(cohort), iso.raw, iso.harm
        ));
        for p in &iso.probs {
            let cls = if *p >= 0.5 { " class=\"hot\"" } else if *p <= 0.05 { " class=\"cold\"" } else { " class=\"warm\"" };
            h.push_str(&format!("<td{}>{:.2}</td>", cls, p));
        }
        h.push_str("</tr>\n");
    }
    h.push_str("</tbody></table></div>\n");

    h.push_str(&format!(
        "<footer>Generated by phenocal. {} isolates, {} cohorts. Method: de Ruvo, Castelli, Di Pasquale, Radomski (in preparation).</footer>\n",
        c.isolates.len(), c.cohorts.len()
    ));
    h.push_str(SCRIPT);
    h.push_str("</body></html>\n");
    h
}

fn card(h: &mut String, label: &str, value: &str) {
    h.push_str(&format!(
        "<div class=\"card\"><div class=\"lbl\">{}</div><div class=\"val\">{}</div></div>\n",
        label, value
    ));
}

/// Circular layout comparison graph.
fn graph_svg(h: &mut String, c: &Ctx) {
    let n = c.cohorts.len();
    let (w, ht) = (780.0, 480.0);
    let (cx, cy, r) = (w / 2.0, ht / 2.0, 165.0);
    let pos: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let ang = -PI / 2.0 + 2.0 * PI * (i as f64) / (n as f64);
            (cx + r * ang.cos(), cy + r * ang.sin())
        })
        .collect();
    let wmax = c.edges.iter().map(|e| e.weight).fold(1e-9, f64::max);

    h.push_str(&format!(
        "<svg viewBox=\"0 0 {w} {ht}\" class=\"fig\" role=\"img\">\n"
    ));
    // edges first
    for e in c.edges {
        let (x1, y1) = pos[e.a];
        let (x2, y2) = pos[e.b];
        let sw = 1.0 + 7.0 * (e.weight / wmax);
        h.push_str(&format!(
            "<line x1=\"{x1:.1}\" y1=\"{y1:.1}\" x2=\"{x2:.1}\" y2=\"{y2:.1}\" stroke=\"#9bb8c9\" stroke-width=\"{sw:.2}\" stroke-opacity=\"0.8\"/>\n"
        ));
        let (mx, my) = ((x1 + x2) / 2.0, (y1 + y2) / 2.0);
        h.push_str(&format!(
            "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"86\" height=\"30\" rx=\"5\" fill=\"#ffffff\" stroke=\"#cdd8df\"/>\
             <text x=\"{:.1}\" y=\"{:.1}\" class=\"el\">{:.2}&#215; (a/b)</text>\
             <text x=\"{:.1}\" y=\"{:.1}\" class=\"es\">n={}, &#964;={:.0}</text>\n",
            mx - 43.0, my - 15.0,
            mx, my - 3.0, fold(e.delta),
            mx, my + 10.0, e.n, e.tau
        ));
    }
    // nodes
    for (i, name) in c.cohorts.iter().enumerate() {
        let (x, y) = pos[i];
        let anchor = i == c.anchor;
        let fillc = if anchor { "#2A9D8F" } else { "#2E86AB" };
        let cnt = c.counts.get(name).copied().unwrap_or(0);
        h.push_str(&format!(
            "<circle cx=\"{x:.1}\" cy=\"{y:.1}\" r=\"30\" fill=\"{fillc}\"/>\
             <text x=\"{x:.1}\" y=\"{:.1}\" class=\"nn\">{}</text>\
             <text x=\"{x:.1}\" y=\"{:.1}\" class=\"nc\">n={}</text>\n",
            y - 2.0, esc(name), y + 12.0, cnt
        ));
    }
    h.push_str("</svg>\n");
}

/// Forest plot of offsets on a log fold axis.
fn forest_svg(h: &mut String, c: &Ctx) {
    let n = c.cohorts.len();
    let row = 34.0;
    let (w, top, left, right) = (780.0, 20.0, 150.0, 30.0);
    let ht = top + row * (n as f64) + 40.0;
    // log10-fold range
    let mut lo = 0.0f64;
    let mut hi = 0.0f64;
    for i in 0..n {
        lo = lo.min(c.sol.lo95[i]);
        hi = hi.max(c.sol.hi95[i]);
    }
    lo -= 0.1;
    hi += 0.1;
    if (hi - lo).abs() < 1e-6 {
        hi = lo + 1.0;
    }
    let xspan = w - left - right;
    let xof = |d: f64| left + (d - lo) / (hi - lo) * xspan;

    h.push_str(&format!("<svg viewBox=\"0 0 {w} {ht}\" class=\"fig\">\n"));
    // gridlines at nice fold values
    for &f in &[0.1f64, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0] {
        let d = f.log10();
        if d < lo || d > hi {
            continue;
        }
        let x = xof(d);
        let one = (f - 1.0).abs() < 1e-9;
        h.push_str(&format!(
            "<line x1=\"{x:.1}\" y1=\"{top:.1}\" x2=\"{x:.1}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"{}\"/>\
             <text x=\"{x:.1}\" y=\"{:.1}\" class=\"ax\">{}&#215;</text>\n",
            top + row * (n as f64),
            if one { "#E76F51" } else { "#e4e9ec" },
            if one { 1.5 } else { 1.0 },
            top + row * (n as f64) + 16.0,
            trim(f)
        ));
    }
    for i in 0..n {
        let y = top + row * (i as f64) + row / 2.0;
        let anchor = i == c.anchor;
        let col = if anchor { "#2A9D8F" } else { "#2E86AB" };
        let (xl, xm, xh) = (xof(c.sol.lo95[i]), xof(c.sol.delta[i]), xof(c.sol.hi95[i]));
        h.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" class=\"fl\">{}</text>\
             <line x1=\"{xl:.1}\" y1=\"{y:.1}\" x2=\"{xh:.1}\" y2=\"{y:.1}\" stroke=\"{col}\" stroke-width=\"2\"/>\
             <line x1=\"{xl:.1}\" y1=\"{:.1}\" x2=\"{xl:.1}\" y2=\"{:.1}\" stroke=\"{col}\" stroke-width=\"2\"/>\
             <line x1=\"{xh:.1}\" y1=\"{:.1}\" x2=\"{xh:.1}\" y2=\"{:.1}\" stroke=\"{col}\" stroke-width=\"2\"/>\
             <circle cx=\"{xm:.1}\" cy=\"{y:.1}\" r=\"5\" fill=\"{col}\"/>\
             <text x=\"{:.1}\" y=\"{:.1}\" class=\"fv\">{:.2}&#215;</text>\n",
            left - 8.0, y + 4.0, esc(&c.cohorts[i]),
            y - 6.0, y + 6.0,
            y - 6.0, y + 6.0,
            xh + 8.0, y + 4.0, fold(c.sol.delta[i])
        ));
    }
    h.push_str("</svg>\n");
}

/// Before/after median + IQR per cohort on the log10 axis.
fn beforeafter_svg(h: &mut String, c: &Ctx) {
    let n = c.cohorts.len();
    // collect raw-log and harm-log per cohort
    let mut raw: Vec<Vec<f64>> = vec![Vec::new(); n];
    let mut harm: Vec<Vec<f64>> = vec![Vec::new(); n];
    for iso in c.isolates {
        if iso.raw > 0.0 {
            raw[iso.cohort].push(iso.raw.log10());
        }
        if iso.harm > 0.0 {
            harm[iso.cohort].push(iso.harm.log10());
        }
    }
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for v in raw.iter().chain(harm.iter()) {
        for &x in v {
            lo = lo.min(x);
            hi = hi.max(x);
        }
    }
    if !lo.is_finite() {
        lo = -0.5;
        hi = 1.5;
    }
    lo -= 0.1;
    hi += 0.1;
    let row = 40.0;
    let (w, top, left, right) = (780.0, 24.0, 150.0, 60.0);
    let ht = top + row * (n as f64) + 30.0;
    let xspan = w - left - right;
    let xof = |d: f64| left + (d - lo) / (hi - lo) * xspan;

    h.push_str(&format!("<svg viewBox=\"0 0 {w} {ht}\" class=\"fig\">\n"));
    // x ticks (fold of MIC)
    for &f in &[0.1f64, 0.3, 1.0, 3.0, 10.0, 30.0] {
        let d = f.log10();
        if d < lo || d > hi { continue; }
        let x = xof(d);
        h.push_str(&format!(
            "<line x1=\"{x:.1}\" y1=\"{top:.1}\" x2=\"{x:.1}\" y2=\"{:.1}\" stroke=\"#eef2f4\"/>\
             <text x=\"{x:.1}\" y=\"{:.1}\" class=\"ax\">{}</text>\n",
            top + row * (n as f64), top + row * (n as f64) + 14.0, trim(f)
        ));
    }
    for i in 0..n {
        let y = top + row * (i as f64) + row / 2.0;
        h.push_str(&format!("<text x=\"{:.1}\" y=\"{:.1}\" class=\"fl\">{}</text>\n", left - 8.0, y + 4.0, esc(&c.cohorts[i])));
        // raw (hollow) above, harmonised (solid) below
        whisker(h, &mut raw[i], xof, y - 7.0, "#9bb8c9", false);
        whisker(h, &mut harm[i], xof, y + 7.0, "#2E86AB", true);
    }
    // legend
    h.push_str(&format!(
        "<circle cx=\"{:.1}\" cy=\"14\" r=\"4\" fill=\"#fff\" stroke=\"#9bb8c9\" stroke-width=\"2\"/><text x=\"{:.1}\" y=\"18\" class=\"ax\" text-anchor=\"start\">raw</text>\
         <circle cx=\"{:.1}\" cy=\"14\" r=\"4\" fill=\"#2E86AB\"/><text x=\"{:.1}\" y=\"18\" class=\"ax\" text-anchor=\"start\">harmonised</text>\n",
        left + 0.0, left + 8.0, left + 70.0, left + 78.0
    ));
    h.push_str("</svg>\n");
}

fn whisker<F: Fn(f64) -> f64>(h: &mut String, v: &mut [f64], xof: F, y: f64, col: &str, solid: bool) {
    if v.is_empty() {
        return;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = median(v);
    let q1 = quantile(v, 0.25);
    let q3 = quantile(v, 0.75);
    let (x1, xm, x3) = (xof(q1), xof(med), xof(q3));
    h.push_str(&format!(
        "<line x1=\"{x1:.1}\" y1=\"{y:.1}\" x2=\"{x3:.1}\" y2=\"{y:.1}\" stroke=\"{col}\" stroke-width=\"2\"/>",
    ));
    if solid {
        h.push_str(&format!("<circle cx=\"{xm:.1}\" cy=\"{y:.1}\" r=\"4\" fill=\"{col}\"/>"));
    } else {
        h.push_str(&format!("<circle cx=\"{xm:.1}\" cy=\"{y:.1}\" r=\"4\" fill=\"#fff\" stroke=\"{col}\" stroke-width=\"2\"/>"));
    }
    h.push('\n');
}

/// One collapsible edge with the full pair table.
fn edge_block(h: &mut String, c: &Ctx, e: &Edge) {
    let a = esc(&c.cohorts[e.a]);
    let b = esc(&c.cohorts[e.b]);
    let n_used = e.cand.iter().filter(|p| p.used).count();
    h.push_str(&format!(
        "<details class=\"edge\"><summary><b>{a} &harr; {b}</b> &mdash; {:.2}&times; (a/b), \
         n<sub>used</sub>={n_used}/{} &nbsp; d<sub>min</sub>={:.0} &nbsp; &tau;={:.0} &nbsp; SE={:.3} &nbsp; weight={:.1}</summary>\n",
        fold(e.delta), e.cand.len(), e.d_min, e.tau, e.se, e.weight
    ));
    h.push_str(&format!(
        "<table class=\"tbl mini\"><thead><tr><th>sample (a: {a})</th><th>sample (b: {b})</th>\
         <th>distance</th><th>MIC<sub>a</sub></th><th>MIC<sub>b</sub></th><th>log<sub>10</sub>(a/b)</th><th>fold</th><th>in median</th></tr></thead><tbody>\n"
    ));
    for p in &e.cand {
        h.push_str(&format!(
            "<tr class=\"{}\"><td>{}</td><td>{}</td><td>{:.0}</td><td>{:.3}</td><td>{:.3}</td><td>{:.3}</td><td>{:.2}&times;</td><td>{}</td></tr>\n",
            if p.used { "used" } else { "skip" },
            esc(&p.si), esc(&p.sj), p.dist,
            10f64.powf(p.va), 10f64.powf(p.vb),
            p.signed, fold(p.signed),
            if p.used { "&#10003;" } else { "&middot;" }
        ));
    }
    h.push_str(&format!(
        "</tbody><tfoot><tr><td colspan=\"5\">median of highlighted log-ratios &rarr; edge &Delta;</td>\
         <td><b>{:.3}</b></td><td><b>{:.2}&times;</b></td><td></td></tr></tfoot></table></details>\n",
        e.delta, fold(e.delta)
    ));
}

/// Trim trailing zeros from a small float for axis labels.
fn trim(f: f64) -> String {
    let s = format!("{:.2}", f);
    let s = s.trim_end_matches('0').trim_end_matches('.').to_string();
    if s.is_empty() { "0".to_string() } else { s }
}

const STYLE: &str = r#"<style>
:root{--ac:#2E86AB;--gd:#2A9D8F;--bd:#E76F51;}
*{box-sizing:border-box}
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:0;color:#1f2d36;background:#f7fafb;line-height:1.45}
header{background:var(--ac);color:#fff;padding:22px 28px}
header h1{margin:0;font-size:22px}
header .sub{margin:4px 0 0;opacity:.9;font-size:13px}
h2{margin:30px 28px 8px;font-size:17px;border-bottom:2px solid #e4e9ec;padding-bottom:4px}
.note,.params{margin:4px 28px 12px;font-size:13px;color:#5a6b75}
.cards{display:flex;flex-wrap:wrap;gap:12px;padding:18px 28px 0}
.card{background:#fff;border:1px solid #e4e9ec;border-radius:8px;padding:10px 14px;min-width:120px;flex:1}
.card .lbl{font-size:11px;text-transform:uppercase;letter-spacing:.04em;color:#7f9099}
.card .val{font-size:20px;font-weight:600;margin-top:2px}
.dim{font-size:12px;font-weight:400;color:#7f9099}
.fig{display:block;width:calc(100% - 56px);max-width:820px;margin:6px 28px;background:#fff;border:1px solid #e4e9ec;border-radius:8px}
.el{font-size:11px;text-anchor:middle;fill:#1f2d36;font-weight:600}
.es{font-size:9px;text-anchor:middle;fill:#5a6b75}
.nn{font-size:12px;text-anchor:middle;fill:#fff;font-weight:700}
.nc{font-size:9px;text-anchor:middle;fill:#eaf3f5}
.ax{font-size:10px;text-anchor:middle;fill:#7f9099}
.fl{font-size:12px;text-anchor:end;fill:#1f2d36}
.fv{font-size:11px;text-anchor:start;fill:#1f2d36;font-weight:600}
table.tbl{border-collapse:collapse;margin:6px 28px;font-size:13px;background:#fff;width:calc(100% - 56px);max-width:820px}
table.tbl th,table.tbl td{border:1px solid #e4e9ec;padding:5px 9px;text-align:right}
table.tbl th:first-child,table.tbl td:first-child,table.tbl th:nth-child(2),table.tbl td:nth-child(2){text-align:left}
table.tbl thead th{background:#eef3f5;color:#33454f;font-weight:600}
tr.anchor{background:#eafaf6}
.scroll{max-height:460px;overflow:auto;margin:0 28px;width:calc(100% - 56px);max-width:820px;border:1px solid #e4e9ec;border-radius:6px}
.scroll table.tbl{margin:0;width:100%;max-width:none;border:0}
.scroll thead th{position:sticky;top:0}
td.hot{background:#fde4dc;font-weight:600}
td.warm{background:#fff6e8}
td.cold{color:#9aa7ae}
details.edge{margin:6px 28px;max-width:820px;background:#fff;border:1px solid #e4e9ec;border-radius:8px}
details.edge summary{cursor:pointer;padding:9px 12px;font-size:13px}
details.edge[open] summary{border-bottom:1px solid #e4e9ec}
table.mini{font-size:12px;margin:0;width:100%;max-width:none}
tr.used{background:#eef7fb}
tr.skip{color:#9aa7ae;background:#fbfcfd}
table.tbl tfoot td{background:#f2f6f8}
#flt{margin:4px 28px 8px;padding:7px 10px;width:calc(100% - 56px);max-width:360px;border:1px solid #cdd8df;border-radius:6px;font-size:13px}
footer{margin:36px 0 0;padding:16px 28px;background:#eef3f5;color:#5a6b75;font-size:12px}
</style>
"#;

const SCRIPT: &str = r#"<script>
function filt(){
 var q=document.getElementById('flt').value.toLowerCase();
 var rows=document.querySelectorAll('#isotbl tbody tr');
 for(var i=0;i<rows.length;i++){
  var t=rows[i].cells[0].textContent.toLowerCase()+' '+rows[i].cells[1].textContent.toLowerCase();
  rows[i].style.display = t.indexOf(q)>=0 ? '' : 'none';
 }
}
</script>
"#;
