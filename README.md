# phenocal

**Genome-anchored calibration of a continuous phenotype across cohorts/batches,
with uncertainty-propagated probabilistic labels.**

`phenocal` is the reference implementation of the genome-anchored
phenotype-calibration method of de Ruvo et al. (*PhenoCal*; see Citation below).
It is organism-, phenotype-, and schema-agnostic:
given a continuous phenotype, a cohort/batch label per sample, and any genomic
distance that can flag near-clonal pairs, it estimates one calibration offset
per cohort and propagates calibration uncertainty into per-cohort credible
intervals and isolate-level probabilistic labels.

Zero external dependencies; builds offline.

## Idea in one line

Two near-clonal isolates from different cohorts share (approximately) the same
biological phenotype, so the difference in their measured values estimates the
difference in cohort/protocol offset — the biology approximately cancels (exactly at zero genomic distance, approximately for near-clones). Offsets are fitted on a
comparison graph by weighted least squares; one cohort is the anchor (δ=0).

## Build

```bash
cargo build --release
# binary: target/release/phenocal
```

## Quick start (synthetic toy)

```bash
phenocal --phenotypes examples/toy_phenotypes.tsv \
         --pairs examples/toy_pairs.tsv \
         --anchor labA --min-support 2 --out /tmp/toy
```

The toy has three labs measuring three shared strains; `labB` reads ~2× and
`labC` ~5× higher than `labA`. `phenocal` recovers that:

```
labA   1.00x  [1.00, 1.00]   (anchor)
labB   1.99x  [1.05, 3.78]
labC   5.01x  [2.64, 9.26]
```

## Inputs

Two TSV files (headers are matched case-insensitively; common aliases accepted).

**Phenotypes** (`--phenotypes`): one row per sample.

| column            | aliases                | meaning                          |
|-------------------|------------------------|----------------------------------|
| `sample`          | accession, id          | sample identifier                |
| `value`           | mic, phenotype         | raw phenotype (log10 by default) |
| `cohort`          | author, batch, study   | cohort / batch label             |

**Pairs** (`--pairs`): genomic distances between sample pairs (within- and
cross-cohort allowed; only cross-cohort pairs are used).

| column      | aliases        | meaning                  |
|-------------|----------------|--------------------------|
| `sample_i`  | acc_i, i       | first sample             |
| `sample_j`  | acc_j, j       | second sample            |
| `distance`  | dist, d        | genomic distance         |

## Run

```bash
phenocal \
  --phenotypes pheno.tsv \
  --pairs pairs.tsv \
  --out result \
  --anchor USA              # optional; default = largest cohort
  --dashboard result.html   # optional; interactive HTML report
```

### Options

| flag | meaning | default |
|---|---|---|
| `--anchor <COHORT>` | cohort fixed to offset 0 | largest cohort |
| `--already-log` | `value` is already log10 | off (raw, log10 applied) |
| `--min-support <N>` | min near-clonal pairs per edge (n_min) | 5 |
| `--max-drift <D>` | max median drift (κ) for τ selection | 1.0 |
| `--sigma <fixed\|empirical>` | per-cohort resolution σ_c: fixed (log10 2) or estimated from each cohort's grid | fixed |
| `--drift-scale <doubling\|rms\|mean\|max>` | units of the τ drift tolerance (combine the two cohort σ); not load-bearing | rms |
| `--lambda <L>` | biological-drift variance per unit genetic distance in the edge SE | log10(2)² |
| `--robust` | robust bootstrap (Student-t edge perturbations) | off (Normal) |
| `--nu <DF>` | degrees of freedom for `--robust` | 4 |
| `--bootstrap <B>` | bootstrap draws for credible intervals | 10000 |
| `--thresholds <list>` | comma-separated tolerance cutoffs (raw scale) | 0.75,1.0,1.25,1.5,2.0 |
| `--clusters <FILE>` | TSV `sample<TAB>cluster` of clonal-cluster ids; enables the **cluster-effective** inference (cluster-weighted point + cluster-resampled 95% interval), the preferred summary when near-clonal pairs are nested in clusters | off (pair-level only) |
| `--no-cutoff` | skip the (costly) data-driven cutoff estimation; still writes harmonised + labels | off |
| `--seed <S>` | RNG seed | 20260603 |
| `--dashboard <FILE>` | also write the interactive HTML report (open locally) | off |
| `--trace <FILE>` | also write a step-by-step **mathematical trace** in Markdown (renders on GitHub, no compilation) | off |

Run `phenocal --help` for the same list.

### Cluster-effective inference (`--clusters`)

When near-clonal pairs are nested within clonal clusters (e.g. NCBI SNP clusters),
pairs from one cluster are not independent and the pair-level interval is
anticonservative. Pass a `sample<TAB>cluster` table:

```
phenocal --phenotypes pheno.tsv --pairs pairs.tsv --clusters clusters.tsv \
         --anchor USA --out out
```

`out.offsets.tsv` then gains three columns — `cluster_wt_fold` (median over
per-cluster medians, one vote per clonal unit), `cluster_eff_lo95` and
`cluster_eff_hi95` (bootstrap resampling clusters, not pairs). This is the
preferred inferential summary when cluster ids are available.

## Outputs

- `<out>.offsets.tsv` — per cohort: δ (log10), fold vs anchor, 95% credible interval.
- `<out>.sigma.tsv` — per cohort: estimated measurement resolution σ_c (log10 and in dilutions) and the σ mode. With `--sigma empirical` this is read off each cohort's own MIC grid (recovers log10 2 ≈ 0.30 for two-fold cohorts, smaller for finer grids).
- `<out>.edges.tsv` — per cohort pair: d_min, selected τ, n pairs, fold, SE, weight (transparency on what supports each edge).
- `<out>.harmonised.tsv` — per sample: raw and harmonised phenotype (`value / 10^δ`).
- `<out>.labels.tsv` — per sample: harmonised value + probability of tolerance `P(T)` at each threshold.
- `<out>.cutoff.tsv` — data-driven sensitive/tolerant cutoff on the harmonised scale, by two methods: **KDE antimode** (non-parametric density valley) and **GMM crossover** (Bayes-optimal boundary of a 2-component Gaussian mixture, with a bootstrap 95% interval). They bracket the cutoff; phenocal does not commit to a single value (the cutoff is treated as uncertain).

### Interactive dashboard (`--dashboard FILE`)

A single self-contained HTML file (no dependencies, no network — opens in any
browser) to inspect the whole run:

- the **comparison graph** drawn in SVG (nodes = cohorts, anchor highlighted, edge thickness ∝ weight);
- a **forest plot** of the offsets with 95% credible intervals;
- **before → after** harmonisation (median + IQR per cohort on the log scale);
- **per-edge detail** — every candidate twin pair with sample IDs, genomic distance, both MICs and the log-ratio, the rows that entered the adaptive-τ median highlighted, and the resulting edge Δ;
- the **per-isolate** harmonised values and probabilistic labels, with a live filter.

### Mathematical trace (`--trace FILE`)

A step-by-step **Markdown** document of the actual run, with the real numbers and
the formula behind each step: inputs and per-cohort σ_c, the per-edge twin tables
and adaptive-τ selection, the WLS normal-equation matrices (AᵀWA, AᵀWΔ), the
solution, the bootstrap credible intervals, and the data-driven cutoff. It is the
didactic worked example auto-filled with this run's data. Because it is Markdown
(with `$...$` math and tables), it **renders directly on GitHub** — no
compilation, no Docker, just open it.

## Algorithm (per cohort pair, per ruler)

1. Collect cross-cohort pairs with distance `d` and signed log-phenotype
   difference `s = y_a − y_b` (oriented `a` over `b`); keep smallest-`d` instance
   per isolate pair.
2. Adaptive `τ`: smallest distance bound with ≥ `min-support` pairs **and**
   median drift ≤ `max-drift · σ̄_ab` from the `d_min` value, where `σ̄_ab` is the
   edge resolution set by `--drift-scale` (`= log10(2)` for two-fold grids).
3. Edge value `Δ = median(s | d ≤ τ)`; `SE = sqrt((σ_a² + σ_b² + λ·τ)/n)` with
   per-cohort resolution `σ_c` (`--sigma`); `w = 1/SE²`. (Reduces to
   `log10(2)·sqrt((1+τ)/n)` when all cohorts are two-fold.)
4. Weighted least squares on the comparison graph for cohort offsets, anchor δ=0.
5. Bootstrap (sample `Δ* ~ N(Δ, SE²)`, refit) → pair-level 95% credible intervals.
   With `--clusters`, also resample clonal clusters (not pairs) per edge and refit
   → **cluster-weighted point** + **cluster-effective** interval, the preferred
   inferential summary when near-clonal pairs are nested in clusters.
6. Per threshold `c`: `P(T_i(c)) = mean over draws of [ MIC_i / 10^δ ≥ c ]`.

The offsets are identifiable up to an additive constant (gauge); the anchor is a
coordinate convention and does not affect harmonised labels (Lemma 1 of the
paper). A connected graph is required.

## Validation

The PhenoCal manuscript validates the method on a public *Salmonella enterica* /
ciprofloxacin-MIC dataset (9,854 isolates assembled from NCBI Pathogen Detection;
see the paper). Primary validation is a controlled realistic benchmark with known
protocol offsets and biological mismatch, where `phenocal` recovers the imposed
offset markedly better than naive cohort median-shift; an injected-shift experiment
on the real genomes recovers known 2×/4×/8× offsets, and the propagated
probabilistic labels recover the true above-cutoff status. The public surveillance
cohorts (USA, Mexico, South Korea) serve as a weak diagnostic negative control.

The implementation ships an 11-test suite (offset recovery, anchor invariance,
edge selection, cutoff, cluster-effective recovery) run under CI; `cargo test`
passes.

## Citation

de Ruvo A, Castelli P, Di Pasquale A, Radomski N.
*PhenoCal: Genome-Anchored Probabilistic Calibration of Continuous Phenotypes
Across Microbial Genomics Cohorts* (in preparation).
Corresponding author: A. de Ruvo (`a.deruvo@izs.it`).

## License

MIT.
