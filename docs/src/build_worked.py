"""Build docs/worked.html: the worked example, every calculation by hand.

    python3 docs/src/build_worked.py docs/worked.html target/release/phenocal

The example is examples/worked_{profiles,phenotypes}.tsv. This script redoes
every step of the algorithm in plain Python (distances, pairs, tau, edge
values, weighted least squares, uncertainty, labels), runs phenocal on the
same files, and stops if any number shown differs from phenocal's output.
"""
import math
import os
import re
import subprocess
import sys
import tempfile

OUT, PHENOCAL = sys.argv[1], sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PROF = os.path.join(ROOT, "examples", "worked_profiles.tsv")
PHEN = os.path.join(ROOT, "examples", "worked_phenotypes.tsv")
MIN_SUPPORT = 3
KAPPA = 1.0
L = math.log10
LOG2 = L(2)
SIGMA = LOG2 / math.sqrt(2)  # fixed per-cohort resolution
LAMBDA = LOG2 ** 2
THRESHOLDS = [0.75, 1.0, 1.25, 1.5, 2.0]
MISSING = {"", "-", "?", "N", "NA", "LNF", "PLOT3", "PLOT5", "LOTSC", "NIPH",
           "NIPHEM", "ALM", "ASM", "PAMA", "EXC", "INV"}


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;")


def f(x, d=3):
    """Number with a real minus sign."""
    s = f"{x:.{d}f}"
    if s.startswith("-") and float(s) == 0:
        s = s[1:]
    return s.replace("-", "\u2212")


def tf(x, d=3):
    """Number for TeX."""
    s = f"{x:.{d}f}"
    if s.startswith("-") and float(s) == 0:
        s = s[1:]
    return s


def median(v):
    s = sorted(v)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def phi(z):
    return 0.5 * (1 + math.erf(z / math.sqrt(2)))


# ------------------------------------------------------------------ inputs
ph_rows = [l.rstrip("\n").split("\t") for l in open(PHEN)][1:]
value = {r[0]: float(r[1]) for r in ph_rows}
cohort = {r[0]: r[2] for r in ph_rows}
cohorts = []
for r in ph_rows:
    if r[2] not in cohorts:
        cohorts.append(r[2])
cidx = {c: i for i, c in enumerate(cohorts)}
y = {s: L(v) for s, v in value.items()}

pr_rows = [l.rstrip("\n").split("\t") for l in open(PROF)]
loci = pr_rows[0][1:]
calls_raw = {r[0]: r[1:] for r in pr_rows[1:]}
order = [r[0] for r in pr_rows[1:]]


def allele(cell):
    s = cell.strip()
    if s.startswith("INF-"):
        s = s[4:]
    if s.isdigit():
        s = s.lstrip("0")
        return ("int", s) if s else None
    return None if s.upper() in MISSING else ("text", s)


calls = {s: [allele(c) for c in v] for s, v in calls_raw.items()}


def dist(a, b):
    return sum(1 for x, z in zip(calls[a], calls[b]) if x and z and x != z)


# cross-cohort pairs, oriented so that the first cohort comes first
pairs = []
for i, a in enumerate(order):
    for b in order[i + 1:]:
        if cohort[a] == cohort[b]:
            continue
        if cidx[cohort[a]] > cidx[cohort[b]]:
            a2, b2 = b, a
        else:
            a2, b2 = a, b
        pairs.append((a2, b2, dist(a2, b2), y[a2] - y[b2]))
n_pairs = len(pairs)

# ------------------------------------------------------------------ edges
edges = {}
for ca in range(len(cohorts)):
    for cb in range(ca + 1, len(cohorts)):
        ps = sorted([p for p in pairs if cidx[cohort[p[0]]] == ca and cidx[cohort[p[1]]] == cb],
                    key=lambda p: p[2])
        if not ps:
            continue
        dmin = ps[0][2]
        base = median([p[3] for p in ps if p[2] == dmin])
        sbar = math.sqrt((SIGMA ** 2 + SIGMA ** 2) / 2)
        steps = []
        chosen = None
        for tau in sorted({p[2] for p in ps}):
            vals = [p[3] for p in ps if p[2] <= tau]
            m = median(vals)
            drift = abs(m - base) / sbar
            ok_n, ok_d = len(vals) >= MIN_SUPPORT, drift <= KAPPA
            steps.append((tau, len(vals), m, drift, ok_n, ok_d))
            if ok_n and ok_d:
                chosen = tau
                break
        used = [p for p in ps if p[2] <= chosen]
        n = len(used)
        delta = median([p[3] for p in used])
        se = math.sqrt((SIGMA ** 2 + SIGMA ** 2 + LAMBDA * chosen) / n)
        edges[(ca, cb)] = dict(pairs=ps, dmin=dmin, base=base, sbar=sbar, steps=steps,
                               tau=chosen, n=n, delta=delta, se=se, w=1 / se ** 2)

# ------------------------------------------------------------------ WLS (anchor = cohort 0)
free = [1, 2]
M = [[0.0, 0.0], [0.0, 0.0]]
rhs = [0.0, 0.0]
for (a, b), e in edges.items():
    for c, sign in ((a, 1), (b, -1)):
        if c in free:
            i = free.index(c)
            M[i][i] += e["w"]
            rhs[i] += sign * e["w"] * e["delta"]
            other = b if c == a else a
            if other in free:
                M[i][free.index(other)] -= e["w"]
det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
inv = [[M[1][1] / det, -M[0][1] / det], [-M[1][0] / det, M[0][0] / det]]
sol = [inv[0][0] * rhs[0] + inv[0][1] * rhs[1], inv[1][0] * rhs[0] + inv[1][1] * rhs[1]]
delta_hat = [0.0, sol[0], sol[1]]
resid = {k: e["delta"] - (delta_hat[k[0]] - delta_hat[k[1]]) for k, e in edges.items()}
rmse = math.sqrt(sum(edges[k]["w"] * r * r for k, r in resid.items()) / sum(e["w"] for e in edges.values()))
sd_an = [0.0, math.sqrt(inv[0][0]), math.sqrt(inv[1][1])]

# ------------------------------------------------------------------ run phenocal and compare
tmp = tempfile.mkdtemp()
pre = os.path.join(tmp, "w")
subprocess.run([PHENOCAL, "--phenotypes", PHEN, "--profiles", PROF, "--anchor", cohorts[0],
                "--min-support", str(MIN_SUPPORT), "--out", pre, "--trace", pre + ".md"],
               check=True, capture_output=True)


def table(name):
    rows = [l.rstrip("\n").split("\t") for l in open(f"{pre}.{name}.tsv")]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def check(what, ours, theirs, tol):
    if abs(ours - theirs) > tol:
        sys.exit(f"MISMATCH {what}: computed {ours} vs phenocal {theirs}")


for r in table("edges"):
    e = edges[(cidx[r["cohort_a"]], cidx[r["cohort_b"]])]
    check("d_min", e["dmin"], float(r["d_min"]), 1e-9)
    check("tau", e["tau"], float(r["tau"]), 1e-9)
    check("n", e["n"], int(r["n_pairs"]), 0)
    check("fold", 10 ** e["delta"], float(r["fold_a_over_b"]), 5e-5)
    check("se", e["se"], float(r["se_log10"]), 5e-5)
    check("weight", e["w"], float(r["weight"]), 5e-4)
trace = open(pre + ".md").read()
wls_tab = trace.split("Solving gives the offsets")[1].split("Weighted residual")[0]
for c in cohorts[1:]:
    m = re.search(rf"`{c}` \| (-?[0-9.]+) \|", wls_tab)
    check("wls " + c, delta_hat[cidx[c]], float(m.group(1)), 5e-5)
check("rmse", rmse, float(re.search(r"RMSE \$=([0-9.]+)\$", trace).group(1)), 5e-5)
n_trace = int(re.search(r"\((\d+) cross-cohort pairs", trace).group(1))
check("pairs", n_pairs, n_trace, 0)
off = {r["cohort"]: r for r in table("offsets")}
rep_delta = {c: float(off[c]["delta_log10"]) for c in cohorts}
rep_lo = {c: L(float(off[c]["fold_lo95"])) for c in cohorts}
rep_hi = {c: L(float(off[c]["fold_hi95"])) for c in cohorts}
for c in cohorts[1:]:  # Monte Carlo agrees with the analytic interval
    i = cidx[c]
    check("mc point " + c, delta_hat[i], rep_delta[c], 0.01)
    check("mc lo " + c, delta_hat[i] - 1.959964 * sd_an[i], rep_lo[c], 0.03)
    check("mc hi " + c, delta_hat[i] + 1.959964 * sd_an[i], rep_hi[c], 0.03)
# Monte Carlo spread of the labC upper bound over 30 seeds
his = []
for seed in range(1, 31):
    subprocess.run([PHENOCAL, "--phenotypes", PHEN, "--profiles", PROF, "--anchor", cohorts[0],
                    "--min-support", str(MIN_SUPPORT), "--quick", "--seed", str(seed),
                    "--out", pre + "_s"], check=True, capture_output=True)
    rows = [l.rstrip("\n").split("\t") for l in open(pre + "_s.offsets.tsv")]
    by = {r[0]: dict(zip(rows[0], r)) for r in rows[1:]}
    his.append(L(float(by[cohorts[2]]["fold_hi95"])))
mc_mean = sum(his) / len(his)
mc_sd = math.sqrt(sum((h - mc_mean) ** 2 for h in his) / (len(his) - 1))
check("mc mean hi", delta_hat[2] + 1.959964 * sd_an[2], mc_mean, 0.005)
harm_tab = {r["sample"]: float(r["value_harmonised"]) for r in table("harmonised")}
for s in value:
    check("harm " + s, 10 ** (y[s] - rep_delta[cohort[s]]), harm_tab[s], 5e-5)
lab_tab = {r["sample"]: r for r in table("labels")}
post_sd = {c: (rep_hi[c] - rep_lo[c]) / (2 * 1.959964) for c in cohorts}


def p_tol(s, t):
    c = cohort[s]
    m = y[s] - rep_delta[c] - L(t)
    if post_sd[c] == 0:
        return 1.0 if m >= 0 else 0.0
    return phi(m / post_sd[c])


for s in value:
    for t in THRESHOLDS:
        key = "P_tol_" + (f"{t:g}".replace(".", "p"))
        check(f"P {s} {t}", p_tol(s, t), float(lab_tab[s][key]), 0.02)

# ------------------------------------------------------------------ html pieces
A, B, C = cohorts
EK = {"AB": (0, 1), "AC": (0, 2), "BC": (1, 2)}
eAB, eAC, eBC = (edges[EK[k]] for k in ("AB", "AC", "BC"))


TABLES = [0]


def tcap(text):
    TABLES[0] += 1
    return f'<p class="tcap"><b>Table S{TABLES[0]}.</b> {text}</p>'


def code(s):
    return f"<code>{esc(s)}</code>"


# phenotypes table: one row per strain
strains = []
for s in order:
    st = s.split("_")[0]
    if st not in strains:
        strains.append(st)
ph_html = []
for st in strains:
    cells = []
    for c in cohorts:
        s = f"{st}_{c[-1]}"
        cells.append(f'<td class="num">{value[s]:g}</td>' if s in value else '<td class="num muted">·</td>')
    ph_html.append(f"<tr><td>{st}</td>{''.join(cells)}</tr>")
ph_html = "".join(ph_html)

# profiles table
prof_html = ["<tr><th>sample</th>" + "".join(f"<th>{i + 1}</th>" for i in range(len(loci))) + "</tr>"]
for s in order:
    cells = []
    for raw, al in zip(calls_raw[s], calls[s]):
        cls = "miss" if al is None else ("inf" if raw.startswith("INF-") else ("mut" if int(raw) >= 100 else ""))
        cells.append(f'<td class="{cls}">{esc(raw)}</td>')
    prof_html.append(f"<tr><td>{s}</td>{''.join(cells)}</tr>")
prof_html = "".join(prof_html)


def locus_table(a, b):
    rows = []
    for j, (ra, rb) in enumerate(zip(calls_raw[a], calls_raw[b])):
        xa, xb = calls[a][j], calls[b][j]
        if xa is None or xb is None:
            verdict, cls = "not compared (missing)", "miss"
        elif xa != xb:
            verdict, cls = "different: +1", "mut"
        else:
            verdict, cls = "same", ""
        rows.append(f'<tr class="{cls}"><td>{j + 1}</td><td>{esc(ra)}</td><td>{esc(rb)}</td><td>{verdict}</td></tr>')
    return "".join(rows)


def pairs_table(e, limit):
    rows = []
    for a, b, d, s in e["pairs"]:
        if d > limit:
            continue
        used = d <= e["tau"]
        rows.append(f'<tr class="{"used" if used else ""}"><td>{a}</td><td>{b}</td><td class="num">{d}</td>'
                    f'<td class="num">{value[a]:g}</td><td class="num">{value[b]:g}</td>'
                    f'<td class="num">{f(s)}</td><td>{"yes" if used else ""}</td></tr>')
    rest = sum(1 for p in e["pairs"] if p[2] > limit)
    dmin_rest = min((p[2] for p in e["pairs"] if p[2] > limit), default=None)
    return "".join(rows), rest, dmin_rest


def tau_table(e):
    rows = []
    for tau, n, m, drift, okn, okd in e["steps"]:
        chosen = tau == e["tau"]
        why = "chosen" if chosen else ("too few pairs" if not okn else "median drifts")
        rows.append(f'<tr class="{"used" if chosen else ""}"><td class="num">{tau}</td><td class="num">{n}</td>'
                    f'<td class="num">{f(m)}</td><td class="num">{f(drift, 2)}</td>'
                    f'<td>{"yes" if okn else "no"}</td><td>{"yes" if okd else "no"}</td><td>{why}</td></tr>')
    return "".join(rows)


def edge_block(key, e, intro):
    a, b = cohorts[EK[key][0]], cohorts[EK[key][1]]
    ptab, rest, drest = pairs_table(e, 3)
    vals = sorted(p[3] for p in e["pairs"] if p[2] <= e["tau"])
    med_expl = ", ".join(tf(v) for v in vals)
    n = e["n"]
    if n % 2:
        mtex = rf"\operatorname{{median}}\{{{med_expl}\}} = {tf(e['delta'])}"
    else:
        mtex = (rf"\operatorname{{median}}\{{{med_expl}\}} \\ &= \frac{{({tf(vals[n // 2 - 1])}) + ({tf(vals[n // 2])})}}{{2}}"
                rf" = {tf(e['delta'])}")
    return f"""
<h3 id="edge-{key}">{code(a)} – {code(b)}</h3>
{intro}
<div class="tablewrap">{tcap(f"Cross-cohort pairs of {esc(a)} and {esc(b)} up to distance 3, with their log-ratio s; highlighted rows enter the median.")}<table>
<thead><tr><th>isolate in {esc(a)}</th><th>isolate in {esc(b)}</th><th class="num">d</th><th class="num">value {esc(a)}</th><th class="num">value {esc(b)}</th><th class="num">s</th><th>used</th></tr></thead>
<tbody>{ptab}</tbody></table></div>
<p class="small">The other {rest} pairs of these two laboratories are different strains, at distance {drest} or more.</p>
<div class="tablewrap">{tcap(f"Candidate bounds τ for {esc(a)}–{esc(b)}, Eq. (4): the first row satisfying both conditions is chosen.")}<table>
<thead><tr><th class="num">τ</th><th class="num">pairs with d ≤ τ</th><th class="num">median of s</th><th class="num">drift / σ̄</th><th>≥ {MIN_SUPPORT} pairs</th><th>drift ≤ κ</th><th></th></tr></thead>
<tbody>{tau_table(e)}</tbody></table></div>
$$ \\begin{{aligned}} \\Delta_{{\\mathrm{{{key[0]}{key[1]}}}}} &= {mtex} \\end{{aligned}} $$
<p>a factor of \\(10^{{{tf(e['delta'])}}} = {10 ** e['delta']:.2f}\\), and</p>
$$ \\mathrm{{SE}} = \\sqrt{{\\frac{{{SIGMA ** 2:.5f} + {SIGMA ** 2:.5f} + {LAMBDA:.5f}\\times {e['tau']}}}{{{n}}}}}
= {e['se']:.4f}, $$
$$ w = \\frac{{1}}{{\\mathrm{{SE}}^2}} = {e['w']:.2f}. $$
"""


# figure: tau selection on the A-C edge
def fig_tau():
    e = eAC
    x0, x1, yt, yb = 150, 660, 40, 250
    dmax = 4
    lo, hi = -1.0, 0.15
    X = lambda d: x0 + d / dmax * (x1 - x0)
    Y = lambda s: yt + (hi - s) / (hi - lo) * (yb - yt)
    o = ['<svg viewBox="0 0 700 320" role="img" aria-labelledby="ft"><title id="ft">Log-ratios of the '
         f'{esc(A)}–{esc(C)} pairs against genomic distance, with the candidate bounds τ</title>']
    # zones
    o.append(f'<rect class="usedzone" x="{X(-0.3):.1f}" y="{yt}" width="{X(e["tau"] + 0.3) - X(-0.3):.1f}" height="{yb - yt}"/>')
    for v in (0.0, -0.301, -0.602, -0.903):
        o.append(f'<line class="track" x1="{x0 - 40}" x2="{x1}" y1="{Y(v):.1f}" y2="{Y(v):.1f}"/>')
        o.append(f'<text class="tick" x="{x0 - 50}" y="{Y(v) + 4:.1f}" text-anchor="end">{f(v)}</text>')
    o.append(f'<line class="axis" x1="{x0 - 40}" x2="{x1}" y1="{yb}" y2="{yb}"/>')
    for d in range(dmax + 1):
        o.append(f'<text class="tick" x="{X(d):.1f}" y="{yb + 18}" text-anchor="middle">{d}</text>')
    o.append(f'<text class="note" x="{(x0 + x1) / 2:.0f}" y="{yb + 40}" text-anchor="middle">genomic distance d (alleles)</text>')
    o.append(f'<text class="note" x="18" y="{(yt + yb) / 2:.0f}" transform="rotate(-90 18 {(yt + yb) / 2:.0f})" text-anchor="middle">log-ratio s</text>')
    # points, jittered when they coincide
    near = [p for p in e["pairs"] if p[2] <= dmax]
    count = {}
    for a, b, d, s in near:
        count[(d, round(s, 6))] = count.get((d, round(s, 6)), 0) + 1
    seen = {}
    for a, b, d, s in near:
        k = (d, round(s, 6))
        j = seen.get(k, 0)
        seen[k] = j + 1
        dx = (j - (count[k] - 1) / 2) * 14
        cls = "pt used" if d <= e["tau"] else "pt"
        o.append(f'<circle class="{cls}" cx="{X(d) + dx:.1f}" cy="{Y(s):.1f}" r="6"/>')
    # medians for each candidate tau
    for tau, n, m, drift, okn, okd in e["steps"]:
        col = "med" if tau == e["tau"] else "tau"
        o.append(f'<line class="{col}" x1="{X(tau) - 26:.1f}" x2="{X(tau) + 26:.1f}" y1="{Y(m):.1f}" y2="{Y(m):.1f}"/>')
        lab = "chosen" if tau == e["tau"] else ("n < 3" if not okn else "drifts")
        o.append(f'<text class="note" x="{X(tau):.1f}" y="{yt - 14}" text-anchor="middle">τ = {tau}: {lab}</text>')
    o.append("</svg>")
    return "".join(o)


# figure: the comparison graph
def fig_graph():
    P = {A: (120, 250), B: (350, 60), C: (580, 250)}
    o = ['<svg viewBox="0 0 700 330" role="img" aria-labelledby="fg"><title id="fg">The comparison graph '
         'with the three edge values and weights</title>']
    labels = {(0, 1): (190, 135, "end"), (0, 2): (350, 292, "middle"), (1, 2): (510, 135, "start")}
    for (a, b), e in edges.items():
        pa, pb = P[cohorts[a]], P[cohorts[b]]
        o.append(f'<line class="edge" x1="{pa[0]}" y1="{pa[1]}" x2="{pb[0]}" y2="{pb[1]}" stroke-width="{1 + e["w"] / 8:.1f}"/>')
        lx, ly, anc = labels[(a, b)]
        o.append(f'<text class="note" x="{lx}" y="{ly}" text-anchor="{anc}">Δ = {f(e["delta"])} ({10 ** e["delta"]:.2f}×)</text>')
        o.append(f'<text class="tick" x="{lx}" y="{ly + 18}" text-anchor="{anc}">τ = {e["tau"]}, n = {e["n"]}, w = {e["w"]:.1f}</text>')
    for c, cls in zip(cohorts, ("nA", "nB", "nC")):
        x, yy = P[c]
        o.append(f'<circle class="node {cls}" cx="{x}" cy="{yy}" r="38"/>')
        o.append(f'<text class="lab" x="{x}" y="{yy + 5}" text-anchor="middle">{esc(c)}</text>')
    o.append("</svg>")
    return "".join(o)


# labels table
lab_rows = []
for s in sorted(value, key=lambda s: (cidx[cohort[s]], int(s.split("_")[0][1:]))):
    c = cohort[s]
    cells = "".join(f'<td class="num">{float(lab_tab[s]["P_tol_" + f"{t:g}".replace(".", "p")]):.2f}</td>'
                    for t in THRESHOLDS)
    lab_rows.append(f'<tr><td>{s}</td><td>{c}</td><td class="num">{value[s]:g}</td>'
                    f'<td class="num">{harm_tab[s]:.3f}</td>{cells}</tr>')
lab_rows = "".join(lab_rows)

ex = "S11_B"  # the isolate worked by hand
ex_c = cohort[ex]
ex_m = y[ex] - rep_delta[ex_c] - L(1.25)
ex_p = p_tol(ex, 1.25)
ex_tab = float(lab_tab[ex]["P_tol_1p25"])

S7 = locus_table("S7_A", "S7_C")
S1 = locus_table("S1_A", "S1_B")

cyc = eAB["delta"] + eBC["delta"] - eAC["delta"]
w1, w2, w3 = eAB["w"], eAC["w"], eBC["w"]

body = f"""
<header class="paper">
  <p class="kicker">Supplementary material</p>
  <h1>Worked example: the calibration computed by hand</h1>
  <p class="affil-note">Supplement to <a href="./">Genome-anchored calibration of a phenotype
  measured across laboratories</a>. Every number below is recomputed by the script that generates
  this page and checked against the output of phenocal.</p>
  <nav class="links">
    <a href="./">Main text</a>
    <a href="https://github.com/genpat-it/phenocal/blob/main/examples/worked_profiles.tsv">Profiles</a>
    <a href="https://github.com/genpat-it/phenocal/blob/main/examples/worked_phenotypes.tsv">Phenotypes</a>
    <a href="https://github.com/genpat-it/phenocal/blob/main/docs/src/build_worked.py">Generating script</a>
  </nav>
</header>
<main>

<nav class="toc supp" aria-label="Contents"><span class="h">Contents</span><ol>
<li><a href="#data">The data</a></li>
<li><a href="#distances">From allelic profiles to distances</a></li>
<li><a href="#pairs">Cross-cohort pairs and log-ratios</a></li>
<li><a href="#edges">One edge per pair of laboratories: τ, Δ, SE</a></li>
<li><a href="#wls">Putting the edges together: weighted least squares</a></li>
<li><a href="#uncertainty">Uncertainty</a></li>
<li><a href="#harmonise">Harmonised values and probabilistic labels</a></li>
<li><a href="#run">Run it yourself</a></li>
</ol></nav>

<p>The main text derives the method; this supplement follows one data set through every
step of it, by hand. The example is deliberately small, but it is constructed so that each
step has an effect: the bound τ must widen once for lack of pairs and is once stopped by
biological drift, and the three laboratories do not agree perfectly, so that the
least-squares solution is a genuine compromise between them. Equations of the main text are
cited as Eq. (1)–(9); those of this supplement are numbered (S1), (S2), ….</p>

<h2 id="data"><span class="n">S1</span>The data</h2>
<p>Three laboratories ({code(A)}, {code(B)}, {code(C)}) measured a phenotype (an MIC, in
mg/L) on 24 isolates of 12 strains. Each strain was measured by two of the laboratories.</p>
<div class="tablewrap">{tcap("Values measured by each laboratory (mg/L); a dot means the strain was not measured there.")}<table class="narrow">
<thead><tr><th>strain</th>{''.join(f'<th class="num">{esc(c)}</th>' for c in cohorts)}</tr></thead>
<tbody>{ph_html}</tbody></table></div>
<p>The file <code>examples/worked_phenotypes.tsv</code> has one row per isolate,
<code>S1_A&nbsp;&nbsp;1&nbsp;&nbsp;labA</code> and so on; the isolate name is the strain followed
by the laboratory.</p>
<p>The genomes are given as cgMLST allelic profiles, 12 loci per isolate, as chewBBACA @@CITE:chewbbaca@@ writes
them (<code>examples/worked_profiles.tsv</code>). Each number is the id of the allele found
at that locus. Highlighted: <span class="k inf">INF-</span> alleles (new alleles inferred by
chewBBACA), <span class="k miss">missing calls</span> and <span class="k mut">alleles that differ from the
rest of the strain</span>.</p>
<div class="tablewrap">{tcap("Allelic profiles, 12 loci (columns) per isolate.")}<table class="prof">{prof_html}</table></div>

<h2 id="distances"><span class="n">S2</span>From allelic profiles to distances</h2>
<p>The distance between two isolates is the number of loci where <b>both</b> have a call
and the two alleles differ. Each call is read as follows:</p>
<ul>
<li><code>INF-1</code> is allele 1: chewBBACA adds the prefix when the allele is new, but it is
an allele like any other.</li>
<li><code>LNF</code>, <code>ASM</code>, <code>NIPH</code>, <code>PLOT3</code>, … and <code>0</code>
or an empty cell are missing calls. A locus missing in either isolate is skipped, so a
missing call never counts as a difference.</li>
</ul>
<p>Isolates <code>S7_A</code> and <code>S7_C</code>, locus by locus:</p>
<div class="tablewrap">{tcap("Isolates S7_A and S7_C compared locus by locus, Eq. (3).")}<table class="narrow">
<thead><tr><th>locus</th><th>S7_A</th><th>S7_C</th><th></th></tr></thead>
<tbody>{S7}</tbody></table></div>
<p>Eleven loci are compared and two differ, so \\(d(\\mathrm{{S7\\_A}}, \\mathrm{{S7\\_C}}) = 2\\).
For <code>S1_A</code> and <code>S1_B</code> the only unusual call is <code>INF-1</code>, which is
allele 1, so \\(d = 0\\): the two isolates have identical profiles. Isolates of different
strains differ at every locus they share, so their distance is 11 or 12.</p>
<p>This is the distance of Eq. (3), computed in the same way by chewBBACA, cgmlst-dists
@@CITE:cgmlstdists@@ and ReporTree @@CITE:reportree@@. phenocal
computes it only for the isolates that have a phenotype, and only between isolates of
different laboratories, since those are the only pairs it uses.</p>

<h2 id="pairs"><span class="n">S3</span>Cross-cohort pairs and log-ratios</h2>
<p>Values are compared on the \\(\\log_{{10}}\\) scale, where a constant factor between
laboratories becomes a constant difference:</p>
$$ y = \\log_{{10}}(\\text{{value}}), \\qquad y_{{ic}} = \\mu_i + \\delta_c + \\varepsilon_{{ic}}. \\tag{{1}} $$
<p>\\(\\mu_i\\) is the true value of strain \\(i\\), \\(\\delta_c\\) the offset of laboratory \\(c\\)
(what we want), \\(\\varepsilon\\) the measurement noise. For a pair of isolates from
laboratories \\(a\\) and \\(b\\) the <b>log-ratio</b> is</p>
$$ s = y_a - y_b = \\log_{{10}}\\frac{{\\text{{value}}_a}}{{\\text{{value}}_b}}. \\tag{{2}} $$
<p>If the two isolates are the same strain, \\(\\mu\\) cancels and \\(s\\) measures
\\(\\delta_a - \\delta_b\\) plus noise. For example <code>S1_A</code> = 1 and <code>S1_B</code> = 2:
\\(s = \\log_{{10}}(1/2) = {tf(L(0.5))}\\), so on this strain {code(B)} reads twice as high
as {code(A)}.</p>
<p>The 9 isolates of {code(A)}, 6 of {code(B)} and 9 of {code(C)} give
\\(9\\times 6 + 9\\times 9 + 6\\times 9 = {n_pairs}\\) cross-cohort pairs, grouped by the two
laboratories they connect. Pairs within one laboratory say nothing about the offsets and
are not used.</p>

<h2 id="edges"><span class="n">S4</span>One edge per pair of laboratories: τ, Δ, SE</h2>
<p>For each pair of laboratories phenocal decides which pairs to trust, summarises them in a
single value \\(\\Delta\\), and attaches a standard error. Three quantities are fixed by the
options of the run (here the defaults, with <code>--min-support 3</code>):</p>
<ul>
<li>\\(\\sigma = \\log_{{10}}(2)/\\sqrt 2 = {SIGMA:.4f}\\), the measurement resolution of a laboratory
working with two-fold dilutions (<code>--sigma fixed</code>);</li>
<li>\\(\\lambda = \\log_{{10}}(2)^2 = {LAMBDA:.4f}\\), how much variance one allele of distance
adds (<code>--lambda</code>);</li>
<li>\\(\\kappa = 1\\), how far the median may drift, in units of
\\(\\bar\\sigma = \\sqrt{{(\\sigma_a^2+\\sigma_b^2)/2}} = {SIGMA:.4f}\\) (<code>--max-drift</code>).</li>
</ul>
<p><b>The rule for τ.</b> Sort the pairs by distance. Let \\(m_0\\) be the median of \\(s\\)
over the closest pairs (distance \\(d_{{\\min}}\\)). Try the distinct distances in increasing
order; τ is the first one for which</p>
$$ n(\\tau) = \\#\\{{d \\le \\tau\\}} \\ge {MIN_SUPPORT} \\qquad\\text{{and}}\\qquad
\\frac{{|\\operatorname{{median}}\\{{s : d \\le \\tau\\}} - m_0|}}{{\\bar\\sigma}} \\le \\kappa. \\tag{{4}} $$
<p>The first condition asks for enough pairs; the second stops τ when more distant pairs
start to pull the median away, the signature of biology rather than protocol. The edge value
and its uncertainty are then</p>
$$ \\Delta = \\operatorname{{median}}\\{{s : d \\le \\tau\\}}, \\qquad
\\mathrm{{SE}} = \\sqrt{{\\frac{{\\sigma_a^2 + \\sigma_b^2 + \\lambda\\,\\tau}}{{n}}}}, \\qquad w = \\frac{{1}}{{\\mathrm{{SE}}^2}}. \\tag{{5, 6}} $$
<p>A median of an even number of values is the mean of the two middle ones.</p>

{edge_block("AB", eAB, f'''<p>One strain (S1) is identical in the two laboratories; two more (S2, S3) differ by one
allele. At τ = 0 there is a single pair, fewer than {MIN_SUPPORT}, so τ widens to 1. The median
does not move, so τ = 1 is accepted. Note that S3 reads the same in both laboratories:
one noisy pair does not move a median of three.</p>''')}

{edge_block("AC", eAC, f'''<p>The interesting edge. S4 and S5 are identical in the two laboratories and both say that
{esc(C)} reads four times higher (\\(s = {tf(L(0.25))}\\)). That is two pairs, not enough. At distance 2
come S6 and S7, which read the <i>same</i> in both laboratories: they are not quite the same
strain, and their biology differs. With them the median moves by a whole dilution (0.301), 1.41 times
σ̄, more than κ = 1 allows, so τ = 2 is rejected. At τ = 3 two more pairs bring the median back, and τ = 3
is accepted with six pairs.</p>''')}

<figure>{fig_tau()}
<figcaption><b>Figure S1.</b> The {esc(A)}–{esc(C)} pairs up to distance 4. Dashes mark the median of the pairs with
d ≤ τ for each candidate τ; τ = 2 is rejected because its median (−0.301) drifts from the
median of the identical pairs (−0.602). The shaded band holds the pairs used.</figcaption></figure>

{edge_block("BC", eBC, f'''<p>Three strains (S10, S11, S12) are identical in the two laboratories, so τ = 0. Two
say that {esc(C)} reads four times higher than {esc(B)}, one says twice as high; the median is
\\(\\log_{{10}}(1/4)\\). With τ = 0 and three pairs this is the most precise edge, with
twice the weight of the others.</p>''')}

<h2 id="wls"><span class="n">S5</span>Putting the edges together: weighted least squares</h2>
<figure>{fig_graph()}
<figcaption><b>Figure S2.</b> The comparison graph: one node per laboratory, one edge per pair of laboratories
that share strains. Line width is proportional to the weight w.</figcaption></figure>
<p>Each edge is an equation in the offsets, \\(\\delta_a - \\delta_b \\approx \\Delta_{{ab}}\\). Fixing the
anchor \\(\\delta_A = 0\\) (any laboratory can be the anchor; it only sets the reference scale):</p>
$$ \\begin{{aligned}}
0 - \\delta_B &\\approx {tf(eAB['delta'])} && (w = {w1:.2f})\\\\
0 - \\delta_C &\\approx {tf(eAC['delta'])} && (w = {w2:.2f})\\\\
\\delta_B - \\delta_C &\\approx {tf(eBC['delta'])} && (w = {w3:.2f})
\\end{{aligned}} $$
<p>Three equations, two unknowns, and they disagree: going round the triangle, the edge values
should add up to zero, but</p>
$$ \\Delta_{{AB}} + \\Delta_{{BC}} - \\Delta_{{AC}} = {tf(eAB['delta'])} + ({tf(eBC['delta'])}) - ({tf(eAC['delta'])}) = {tf(cyc)}. \\tag{{S1}} $$
<p>The {esc(A)}–{esc(B)} edge says {esc(B)} reads twice as high as {esc(A)}; the other two together say it
reads the same. Weighted least squares finds the offsets that minimise the weighted sum of
squared disagreements:</p>
$$ S(\\delta_B, \\delta_C) = w_{{AB}}\\big(\\Delta_{{AB}} + \\delta_B\\big)^2 + w_{{AC}}\\big(\\Delta_{{AC}} + \\delta_C\\big)^2
+ w_{{BC}}\\big(\\Delta_{{BC}} - \\delta_B + \\delta_C\\big)^2. \\tag{{S2}} $$
<p>Setting both derivatives to zero gives the normal equations
\\((A^{{\\top}} W A)\\,\\delta = A^{{\\top}} W \\Delta\\):</p>
$$ \\begin{{bmatrix}} w_{{AB}} + w_{{BC}} & -w_{{BC}} \\\\ -w_{{BC}} & w_{{AC}} + w_{{BC}} \\end{{bmatrix}}
\\begin{{bmatrix}} \\delta_B \\\\ \\delta_C \\end{{bmatrix}} =
\\begin{{bmatrix}} -w_{{AB}}\\Delta_{{AB}} + w_{{BC}}\\Delta_{{BC}} \\\\ -w_{{AC}}\\Delta_{{AC}} - w_{{BC}}\\Delta_{{BC}} \\end{{bmatrix}} $$
$$ \\begin{{bmatrix}} {M[0][0]:.2f} & {M[0][1]:.2f} \\\\ {M[1][0]:.2f} & {M[1][1]:.2f} \\end{{bmatrix}}
\\begin{{bmatrix}} \\delta_B \\\\ \\delta_C \\end{{bmatrix}} =
\\begin{{bmatrix}} {rhs[0]:.3f} \\\\ {rhs[1]:.3f} \\end{{bmatrix}}. \\tag{{S3}} $$
<p>The matrix on the left is the weighted graph Laplacian with the anchor removed. Its
determinant is \\({M[0][0]:.2f}^2 - {abs(M[0][1]):.2f}^2 = {det:.1f}\\), and inverting it,</p>
$$ (A^{{\\top}} W A)^{{-1}} = \\begin{{bmatrix}} {inv[0][0]:.5f} & {inv[0][1]:.5f} \\\\ {inv[1][0]:.5f} & {inv[1][1]:.5f} \\end{{bmatrix}},
\\qquad
\\begin{{bmatrix}} \\hat\\delta_B \\\\ \\hat\\delta_C \\end{{bmatrix}} = \\begin{{bmatrix}} {delta_hat[1]:.4f} \\\\ {delta_hat[2]:.4f} \\end{{bmatrix}}. $$
<p>As factors: {esc(B)} reads \\(10^{{{delta_hat[1]:.4f}}} = {10 ** delta_hat[1]:.2f}\\times\\) and
{esc(C)} \\(10^{{{delta_hat[2]:.4f}}} = {10 ** delta_hat[2]:.2f}\\times\\) higher than {esc(A)}. (Exactly,
\\(\\hat\\delta_B = 0.6\\log_{{10}}2\\) and \\(\\hat\\delta_C = 2.4\\log_{{10}}2\\).) {esc(B)} ends up between
the 2× of its direct edge and the 1× implied by the other two, closer to the latter
because the {esc(B)}–{esc(C)} edge weighs twice as much.</p>
<p>What is left over on each edge:</p>
<div class="tablewrap">{tcap("Edge values, fitted contrasts and residuals of the least-squares solution.")}<table class="narrow">
<thead><tr><th>edge</th><th class="num">Δ</th><th class="num">fitted \\(\\hat\\delta_a - \\hat\\delta_b\\)</th><th class="num">residual</th><th class="num">w</th></tr></thead>
<tbody>{''.join(f'<tr><td>{esc(cohorts[k[0]])}–{esc(cohorts[k[1]])}</td><td class="num">{f(edges[k]["delta"])}</td><td class="num">{f(delta_hat[k[0]] - delta_hat[k[1]])}</td><td class="num">{f(r)}</td><td class="num">{edges[k]["w"]:.2f}</td></tr>' for k, r in resid.items())}</tbody></table></div>
$$ \\text{{weighted RMSE}} = \\sqrt{{\\frac{{\\sum w\\,r^2}}{{\\sum w}}}} = {rmse:.4f}\\ \\log_{{10}}
= {rmse / LOG2:.2f}\\ \\text{{dilutions}}. \\tag{{S4}} $$
<p>A residual of a third of a dilution is large for real data; here it was put in on purpose.
phenocal writes it at the end of every run, as a check on how well the laboratories agree.</p>

<h2 id="uncertainty"><span class="n">S6</span>Uncertainty</h2>
<p>The edge values are uncertain, with standard errors SE. phenocal propagates this by a
parametric bootstrap @@CITE:efron@@: it draws every edge from \\(\\mathcal N(\\Delta, \\mathrm{{SE}}^2)\\), solves the least squares
again, repeats this 10,000 times, and reports the median of the offsets and their
2.5th–97.5th percentiles.</p>
<p>Because the solution is a linear function of the edge values with weights
\\(1/\\mathrm{{SE}}^2\\), the same answer can be written down exactly (Proposition 2 of the main text): the offsets are normal,
with covariance \\((A^{{\\top}} W A)^{{-1}}\\), the matrix inverted above. So</p>
$$ \\mathrm{{sd}}(\\hat\\delta_B) = \\sqrt{{{inv[0][0]:.5f}}} = {sd_an[1]:.4f}, \\qquad
\\mathrm{{sd}}(\\hat\\delta_C) = \\sqrt{{{inv[1][1]:.5f}}} = {sd_an[2]:.4f}, $$
$$ 95\\%:\\ \\hat\\delta \\pm 1.96\\,\\mathrm{{sd}} \\;\\Rightarrow\\;
\\text{{{B}}}: [{10 ** (delta_hat[1] - 1.959964 * sd_an[1]):.2f},\\ {10 ** (delta_hat[1] + 1.959964 * sd_an[1]):.2f}], \\quad
\\text{{{C}}}: [{10 ** (delta_hat[2] - 1.959964 * sd_an[2]):.2f},\\ {10 ** (delta_hat[2] + 1.959964 * sd_an[2]):.2f}]. $$
<p>phenocal's bootstrap gives:</p>
<div class="tablewrap">{tcap("Offsets reported by phenocal (<code>offsets.tsv</code>): bootstrap median and 95% interval, 10,000 draws.")}<table class="narrow">
<thead><tr><th>cohort</th><th class="num">offset (log10)</th><th class="num">fold</th><th class="num">95% interval</th></tr></thead>
<tbody>{''.join(f'<tr><td>{esc(c)}</td><td class="num">{rep_delta[c]:.4f}</td><td class="num">{float(off[c]["fold"]):.2f}×</td><td class="num">[{float(off[c]["fold_lo95"]):.2f}, {float(off[c]["fold_hi95"]):.2f}]</td></tr>' for c in cohorts)}</tbody></table></div>
<p>The differences from the exact values are Monte Carlo error: 10,000 random draws
give each number with some noise, which shrinks with more draws (<code>--bootstrap</code>).
With the default seed the upper bound for {esc(C)} is {float(off[C]["fold_hi95"]):.2f} instead of
{10 ** (delta_hat[2] + 1.959964 * sd_an[2]):.2f}, a rather unlucky draw: over the seeds 1 to 30 the
same bound averages {10 ** mc_mean:.2f}, with a standard deviation of {mc_sd:.4f} in
\\(\\log_{{10}}\\) ({100 * (10 ** mc_sd - 1):.0f}%). The intervals are wide because each edge rests on
three to six pairs; with real collections, edges have tens or hundreds of pairs.</p>
<p>Why a bootstrap at all, if the exact answer is available? Because phenocal also offers a
robust bootstrap with Student-\\(t\\) draws (<code>--robust</code>) and one that resamples clonal
clusters (<code>--clusters</code>), which have no closed form.</p>

<h2 id="harmonise"><span class="n">S7</span>Harmonised values and probabilistic labels</h2>
<p>Each value is brought to the scale of the anchor by removing its laboratory's offset,
using the reported offsets:</p>
$$ \\text{{harmonised}} = \\frac{{\\text{{value}}}}{{10^{{\\hat\\delta_c}}}}, \\qquad\\text{{e.g. }}
\\mathrm{{S11\\_B}}: \\frac{{{value[ex]:g}}}{{10^{{{rep_delta[ex_c]:.4f}}}}} = {harm_tab[ex]:.3f}. $$
<p>Whether an isolate is above a threshold \\(t\\) (say 1.25 mg/L) depends on the offset, which is
uncertain. phenocal treats the offset of laboratory \\(c\\) as normal, centred on
\\(\\hat\\delta_c\\) with standard deviation \\(s_c\\) read off its 95% interval,
\\(s_c = (\\log_{{10}}\\mathrm{{hi}} - \\log_{{10}}\\mathrm{{lo}})/(2 \\times 1.96)\\), and gives the probability</p>
$$ P(\\text{{harmonised}} \\ge t) = \\Phi\\!\\left(\\frac{{y - \\hat\\delta_c - \\log_{{10}} t}}{{s_c}}\\right). $$
<p>For <code>{ex}</code> and \\(t = 1.25\\): \\(s_{{\\mathrm{{B}}}} = ({rep_hi[ex_c]:.4f} - ({rep_lo[ex_c]:.4f}))/3.92 = {post_sd[ex_c]:.4f}\\),</p>
$$ P = \\Phi\\!\\left(\\frac{{{y[ex]:.4f} - {rep_delta[ex_c]:.4f} - {L(1.25):.4f}}}{{{post_sd[ex_c]:.4f}}}\\right)
= \\Phi({ex_m / post_sd[ex_c]:.3f}) = {ex_p:.3f}. $$
<p>phenocal estimates the same probability by drawing the offset 10,000 times and gets
{ex_tab:.3f}. The isolate's harmonised value, {harm_tab[ex]:.2f}, is just above 1.25, but the
offset of {esc(B)} is so uncertain that the label is close to a coin toss. That is the point of
probabilistic labels: they separate isolates whose class is clear from those whose class
depends on the calibration.</p>
<p>All isolates (<code>labels.tsv</code>; isolates of the anchor have a certain offset, so their
probabilities are 0 or 1):</p>
<div class="tablewrap">{tcap("Harmonised values and probabilities P(harmonised ≥ t), Eq. (9), for the default thresholds (<code>labels.tsv</code>).")}<table>
<thead><tr><th>isolate</th><th>cohort</th><th class="num">raw</th><th class="num">harmonised</th>{''.join(f'<th class="num">P(≥ {t:g})</th>' for t in THRESHOLDS)}</tr></thead>
<tbody>{lab_rows}</tbody></table></div>

<h2 id="run"><span class="n">S8</span>Run it yourself</h2>
<p>The two files are in the repository and in the Docker image:</p>
<pre><code>phenocal --phenotypes examples/worked_phenotypes.tsv \\
         --profiles examples/worked_profiles.tsv \\
         --anchor labA --min-support 3 --out worked \\
         --trace worked.md --dashboard worked.html</code></pre>
<p><code>worked.md</code> repeats these calculations for your own data, with the real numbers of
the run; <code>worked.html</code> shows them interactively.</p>
<h2 id="refs">References</h2>
@@REFS@@
</main>
@@FOOTER@@
</body>
</html>
"""

sys.path.insert(0, HERE)
import common  # noqa: E402

extra = """<style>
.crumb { font: .95rem system-ui, sans-serif; margin: 0 0 .4em; }
.small { font-size: .92rem; color: var(--muted); }
table.narrow { width: auto; }
table.prof { font: .85rem ui-monospace, SFMono-Regular, Menlo, monospace; width: auto; }
table.prof td, table.prof th { padding: 3px 7px; text-align: right; }
table.prof td:first-child { text-align: left; }
td.miss, tr.miss td, .k.miss { background: var(--soft); color: var(--muted); }
td.inf, .k.inf { background: var(--accent-soft); }
td.mut, tr.mut td, .k.mut { background: var(--nB); }
.k { padding: .05em .35em; border-radius: 4px; }
tr.used td { background: var(--accent-soft); }
td.muted { color: var(--muted); }
svg .lab { font-size: 15px; font-weight: 600; }
</style>
"""
body = body.replace("@@FOOTER@@", common.FOOTER)
body, order = common.cite(body)
body = body.replace("@@REFS@@", common.references(order))
html = (common.head("Worked example: the calibration computed by hand — phenocal",
                    "Supplementary material: the phenocal calibration worked out by hand on a small "
                    "example, from allelic profiles to probabilistic labels.")
        .replace("</head>", extra + "</head>") + "<body>" + body)
open(OUT, "w").write(html)
print("wrote", OUT, "| all numbers match phenocal")
