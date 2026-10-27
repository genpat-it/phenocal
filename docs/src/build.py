"""Build docs/index.html for the phenocal GitHub Pages site.

    python3 docs/src/build.py docs/index.html docs/src/template.html

Figures are inline SVG whose coordinates are computed here from the toy
example (examples/toy_*.tsv) and the toy run outputs, so every number drawn
is a real one.
"""
import math
import sys

OUT = sys.argv[1]

# ---------------------------------------------------------------- toy data
toy = {  # sample -> (value, cohort, strain)
    "s1_A": (0.5, "labA", 1), "s1_B": (1.0, "labB", 1), "s1_C": (2.5, "labC", 1),
    "s2_A": (1.0, "labA", 2), "s2_B": (2.0, "labB", 2), "s2_C": (5.0, "labC", 2),
    "s3_A": (2.0, "labA", 3), "s3_B": (4.0, "labB", 3), "s3_C": (10.0, "labC", 3),
}
# toy run outputs (phenocal 0.1.2, --anchor labA --min-support 2)
harm = {"s1_A": 0.5, "s1_B": 0.50368, "s1_C": 0.49947,
        "s2_A": 1.0, "s2_B": 1.00736, "s2_C": 0.99895,
        "s3_A": 2.0, "s3_B": 2.01472, "s3_C": 1.99790}
fold_ci = {"labB": (1.0510, 3.7845), "labC": (2.6407, 9.2588)}

L = math.log10


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;")


# ---------------------------------------------------------------- figure 1
def fig_problem():
    x0, x1, lo, hi = 120, 700, L(0.25), L(16)
    X = lambda v: x0 + (L(v) - lo) / (hi - lo) * (x1 - x0)
    rows = [("labA", 50), ("labB", 90), ("labC", 130)]
    rows2 = [("labA", 245), ("labB", 285), ("labC", 325)]
    o = ['<svg viewBox="0 0 720 405" role="img" aria-labelledby="f1t">',
         '<title id="f1t">The same three strains measured by three labs, before and after calibration</title>']
    o.append('<text class="panel" x="10" y="20">Raw values</text>')
    o.append('<text class="panel" x="10" y="215">After calibration (harmonised to labA)</text>')
    for (lab, y), (_, y2) in zip(rows, rows2):
        for yy in (y, y2):
            o.append(f'<line class="track" x1="{x0}" x2="{x1}" y1="{yy}" y2="{yy}"/>')
            o.append(f'<text class="lab" x="{x0 - 14}" y="{yy + 4}" text-anchor="end">{lab}</text>')
    for s, (v, lab, st) in toy.items():
        y = dict(rows)[lab]
        y2 = dict(rows2)[lab]
        o.append(f'<circle class="st{st}" cx="{X(v):.1f}" cy="{y}" r="7"><title>{s}: {v} (raw)</title></circle>')
        o.append(f'<circle class="st{st}" cx="{X(harm[s]):.1f}" cy="{y2}" r="7"><title>{s}: {harm[s]} (harmonised)</title></circle>')
    for yy in (160, 355):
        o.append(f'<line class="axis" x1="{x0}" x2="{x1}" y1="{yy}" y2="{yy}"/>')
        for t in (0.25, 0.5, 1, 2, 4, 8, 16):
            o.append(f'<line class="axis" x1="{X(t):.1f}" x2="{X(t):.1f}" y1="{yy}" y2="{yy + 5}"/>')
            o.append(f'<text class="tick" x="{X(t):.1f}" y="{yy + 19}" text-anchor="middle">{t:g}</text>')
    o.append(f'<text class="tick" x="{(x0 + x1) / 2}" y="400" text-anchor="middle">phenotype value (log scale)</text>')
    lx = 470
    for st in (1, 2, 3):
        o.append(f'<circle class="st{st}" cx="{lx}" cy="15" r="6"/>')
        o.append(f'<text class="tick" x="{lx + 10}" y="19">strain {st}</text>')
        lx += 78
    o.append("</svg>")
    return "\n".join(o)


# ---------------------------------------------------------------- figure 2
def genome(x, y, cls=""):
    # a genome as a bar of allele blocks; identical blocks = identical genome
    cols = [1, 3, 2, 2, 1, 3, 1, 2, 3, 1]
    o = []
    for k, c in enumerate(cols):
        o.append(f'<rect class="al{c} {cls}" x="{x + k * 16}" y="{y}" width="15" height="18" rx="2"/>')
    return "\n".join(o)


def fig_idea():
    o = ['<svg viewBox="0 0 720 230" role="img" aria-labelledby="f2t">',
         '<title id="f2t">Two isolates with the same genome measured in two labs</title>']
    o.append('<text class="panel" x="80" y="24">Lab A</text>')
    o.append('<text class="panel" x="470" y="24">Lab B</text>')
    o.append(genome(80, 45))
    o.append(genome(470, 45))
    o.append('<text class="tick" x="80" y="82">isolate i · allele profile</text>')
    o.append('<text class="tick" x="470" y="82">isolate j · same allele profile (d = 0)</text>')
    o.append('<path class="link" d="M250 54 L460 54" marker-end="url(#arr)" marker-start="url(#arr)"/>')
    o.append('<text class="note" x="355" y="45" text-anchor="middle">genomic distance d = 0</text>')
    o.append('<text class="big" x="80" y="125">y<tspan baseline-shift="sub" font-size="12">i</tspan> = μ + δ<tspan baseline-shift="sub" font-size="12">A</tspan> + ε</text>')
    o.append('<text class="big" x="470" y="125">y<tspan baseline-shift="sub" font-size="12">j</tspan> = μ + δ<tspan baseline-shift="sub" font-size="12">B</tspan> + ε</text>')
    o.append('<text class="tick" x="80" y="146">log of the measured value</text>')
    o.append('<text class="tick" x="470" y="146">same biology μ, other protocol</text>')
    o.append('<rect class="box" x="170" y="170" width="380" height="46" rx="8"/>')
    o.append('<text class="big" x="360" y="199" text-anchor="middle">y<tspan baseline-shift="sub" font-size="12">j</tspan> − y<tspan baseline-shift="sub" font-size="12">i</tspan> ≈ δ<tspan baseline-shift="sub" font-size="12">B</tspan> − δ<tspan baseline-shift="sub" font-size="12">A</tspan>   (μ cancels)</text>')
    o.append('<defs><marker id="arr" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path class="arrowhead" d="M0 0 L10 5 L0 10 z"/></marker></defs>')
    o.append("</svg>")
    return "\n".join(o)


# ---------------------------------------------------------------- figure 3
def fig_tau():
    # illustrative edge: few pairs at d = 0, more as d grows; near pairs scatter
    # around the offset, far pairs drift as biology starts to differ
    pts = [(0, .31), (0, .27), (1, .33), (2, .29), (2, .36), (3, .25), (4, .32), (5, .30),
           (6, .34), (7, .28), (8, .45), (9, .33), (10, .52), (11, .41),
           (12, .58), (13, .36), (14, .66), (15, .49), (16, .71), (17, .55), (18, .80)]
    # the rule of calib.rs build_edge: smallest tau with >= min_support pairs and
    # |median(s | d <= tau) - median(s | d = d_min)| / sigma_bar <= max_drift
    min_support, max_drift, sbar = 5, 1.0, math.log10(2) / math.sqrt(2)
    def median(v):
        v = sorted(v); n = len(v)
        return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2
    dmin = min(d for d, _ in pts)
    base = median([s for d, s in pts if d == dmin])
    tau = None
    for t in sorted({d for d, _ in pts}):
        vals = [s for d, s in pts if d <= t]
        if len(vals) >= min_support and abs(median(vals) - base) / sbar <= max_drift:
            tau = t
            break
    used = [s for d, s in pts if d <= tau]
    med = median(used)
    x0, x1, y0, y1 = 70, 690, 230, 30
    X = lambda d: x0 + d / 19 * (x1 - x0)
    Y = lambda s: y0 - s / 0.9 * (y0 - y1)
    o = ['<svg viewBox="0 0 720 280" role="img" aria-labelledby="f3t">',
         '<title id="f3t">Adaptive relatedness bound tau for one pair of cohorts (illustrative)</title>']
    o.append(f'<rect class="usedzone" x="{x0}" y="{y1}" width="{X(tau + .5) - x0:.1f}" height="{y0 - y1}"/>')
    o.append(f'<line class="axis" x1="{x0}" x2="{x1}" y1="{y0}" y2="{y0}"/>')
    o.append(f'<line class="axis" x1="{x0}" x2="{x0}" y1="{y0}" y2="{y1}"/>')
    for d in range(0, 19, 2):
        o.append(f'<text class="tick" x="{X(d):.1f}" y="{y0 + 17}" text-anchor="middle">{d}</text>')
    for s in (0, .3, .6, .9):
        o.append(f'<text class="tick" x="{x0 - 8}" y="{Y(s) + 4:.1f}" text-anchor="end">{s:g}</text>')
    o.append(f'<text class="tick" x="{(x0 + x1) / 2}" y="{y0 + 38}" text-anchor="middle">genomic distance d (allelic differences)</text>')
    o.append(f'<text class="tick" transform="translate(18 {(y0 + y1) / 2}) rotate(-90)" text-anchor="middle">log-ratio s</text>')
    for d, s in pts:
        cls = "pt used" if d <= tau else "pt"
        o.append(f'<circle class="{cls}" cx="{X(d):.1f}" cy="{Y(s):.1f}" r="5"/>')
    o.append(f'<line class="tau" x1="{X(tau + .5):.1f}" x2="{X(tau + .5):.1f}" y1="{y1}" y2="{y0}"/>')
    o.append(f'<text class="note" x="{X(tau + .5) + 6:.1f}" y="{y1 + 14}">τ = {tau}</text>')
    o.append(f'<line class="med" x1="{x0}" x2="{X(tau + .5):.1f}" y1="{Y(med):.1f}" y2="{Y(med):.1f}"/>')
    o.append(f'<text class="note" x="{X(tau + .5) + 10:.1f}" y="{Y(.1):.1f}">Δ = median of the {len(used)} shaded pairs = {med:.2f}</text>')
    o.append(f'<text class="tick" x="{X(11):.1f}" y="{Y(.86):.1f}" text-anchor="middle">farther pairs: biology drifts in, not used</text>')
    o.append("</svg>")
    return "\n".join(o)


# ---------------------------------------------------------------- figure 4
def fig_graph():
    A, B, C = (120, 170), (600, 170), (360, 40)
    o = ['<svg viewBox="0 0 720 250" role="img" aria-labelledby="f4t">',
         '<title id="f4t">Comparison graph of the toy example</title>']
    def edge(p, q, label, dx, dy, anchor):
        mx, my = (p[0] + q[0]) / 2 + dx, (p[1] + q[1]) / 2 + dy
        return (f'<line class="edge" x1="{p[0]}" y1="{p[1]}" x2="{q[0]}" y2="{q[1]}"/>'
                f'<text class="note" x="{mx}" y="{my}" text-anchor="{anchor}">{label}</text>')
    o.append(edge(A, B, "labB / labA = 2 · 3 pairs at d = 0", 0, 24, "middle"))
    o.append(edge(A, C, "labC / labA = 5", -14, -6, "end"))
    o.append(edge(B, C, "labC / labB = 2.5", 14, -6, "start"))
    for (x, y), name, d, cls in ((A, "labA", "δ = 0 (anchor)", "nA"),
                                  (B, "labB", "δ = log₁₀ 2 · 1.99× [1.05, 3.78]", "nB"),
                                  (C, "labC", "δ = log₁₀ 5 · 5.01× [2.64, 9.26]", "nC")):
        o.append(f'<circle class="node {cls}" cx="{x}" cy="{y}" r="30"/>')
        o.append(f'<text class="nodet" x="{x}" y="{y + 5}" text-anchor="middle">{name}</text>')
        o.append(f'<text class="tick" x="{x}" y="{y + 50 if y > 100 else y - 40 + 0}" text-anchor="middle">{d}</text>' if y > 100 else
                 f'<text class="tick" x="{x + 40}" y="{y + 5}">{d}</text>')
    o.append("</svg>")
    return "\n".join(o)


# ---------------------------------------------------------------- figure 5
def fig_label():
    # s2_B: harmonised 1.00736; labB offset sd from the 95% CrI half-width
    lo, hi = fold_ci["labB"]
    sd = (L(hi) - L(lo)) / (2 * 1.959964)
    mu = L(1.00736)
    thr = L(1.25)
    z = (thr - mu) / sd
    p = 0.5 * math.erfc(z / math.sqrt(2))
    x0, x1, a, b = 60, 690, mu - 3.6 * sd, mu + 3.6 * sd
    X = lambda v: x0 + (v - a) / (b - a) * (x1 - x0)
    pdf = lambda v: math.exp(-0.5 * ((v - mu) / sd) ** 2)
    yb, ht = 190, 150
    pts = [(X(a + (b - a) * k / 200), yb - ht * pdf(a + (b - a) * k / 200)) for k in range(201)]
    curve = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
    tail = [(x, y) for x, y in pts if x >= X(thr)]
    area = f"M{X(thr):.1f},{yb} " + " ".join(f"L{x:.1f},{y:.1f}" for x, y in tail) + f" L{x1},{yb} Z"
    o = ['<svg viewBox="0 0 720 250" role="img" aria-labelledby="f5t">',
         '<title id="f5t">Probability that one isolate is above a threshold</title>']
    o.append(f'<path class="tail" d="{area}"/>')
    o.append(f'<polyline class="curve" points="{curve}"/>')
    o.append(f'<line class="axis" x1="{x0}" x2="{x1}" y1="{yb}" y2="{yb}"/>')
    for t in (0.5, 0.75, 1, 1.25, 1.5, 2):
        if a <= L(t) <= b:
            o.append(f'<line class="axis" x1="{X(L(t)):.1f}" x2="{X(L(t)):.1f}" y1="{yb}" y2="{yb + 5}"/>')
            o.append(f'<text class="tick" x="{X(L(t)):.1f}" y="{yb + 19}" text-anchor="middle">{t:g}</text>')
    o.append(f'<line class="tau" x1="{X(thr):.1f}" x2="{X(thr):.1f}" y1="30" y2="{yb}"/>')
    o.append(f'<text class="note" x="{X(thr) + 6:.1f}" y="44">threshold c = 1.25</text>')
    o.append(f'<text class="big" x="{X(L(1.9)):.1f}" y="{yb - 95}">P = {p:.2f}</text>')
    o.append(f'<text class="tick" x="{X(mu) - 30:.1f}" y="{yb - ht + 8}" text-anchor="end">s2_B harmonised: 1.007</text>')
    o.append(f'<text class="tick" x="{(x0 + x1) / 2}" y="{yb + 40}" text-anchor="middle">harmonised value of s2_B over the uncertainty of the labB offset (log scale)</text>')
    o.append("</svg>")
    return "\n".join(o), p


FLOW = """flowchart TD
    P[/"Phenotypes<br/>sample · value · cohort"/] --> X
    D[/"Genomic distances<br/>pair table or distance matrix"/] --> X
    X["1 · Cross-cohort pairs<br/>genomic distance d<br/>log-ratio s of the two values"] --> T
    subgraph E["per cohort pair"]
        T["2 · Adaptive bound τ<br/>smallest d with enough pairs<br/>and a stable median"] --> V
        V["3 · Edge value Δ<br/>median of s over d ≤ τ<br/>weight 1 / SE²"]
    end
    V --> G["4 · Weighted least squares on the comparison graph<br/>cohort offsets δ, anchor δ = 0"]
    G --> B["5 · Bootstrap the edges and refit<br/>95% credible intervals"]
    B --> H["6 · Harmonised values<br/>value / 10^δ<br/>P(harmonised ≥ c) per threshold"]
    H --> C["Cutoff on the harmonised scale<br/>KDE antimode · GMM crossover"]"""

f5, p_s2b = fig_label()

rows_data = "".join(
    f"<tr><td>{s}</td><td>{c}</td><td class=\"num\">{v:g}</td><td class=\"num\">{harm[s]:.3f}</td></tr>"
    for s, (v, c, _) in toy.items())

html = open(sys.argv[2]).read()
html = (html.replace("@@FIG1@@", fig_problem())
            .replace("@@FIG2@@", fig_idea())
            .replace("@@FIG3@@", fig_tau())
            .replace("@@FIG4@@", fig_graph())
            .replace("@@FIG5@@", f5)
            .replace("@@FLOW@@", esc(FLOW))
            .replace("@@ROWS@@", rows_data)
            .replace("@@PS2B@@", f"{p_s2b:.2f}"))
open(OUT, "w").write(html)
print("wrote", OUT, "P(s2_B >= 1.25) =", round(p_s2b, 4))
