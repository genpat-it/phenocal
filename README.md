# phenocal

**Genome-anchored calibration of a continuous phenotype across cohorts/batches,
with uncertainty-propagated probabilistic labels.**

`phenocal` implements the genome-anchored MIC calibration method developed by
de Ruvo, Castelli, Di Pasquale and Radomski (*Genome-Anchored Probabilistic MIC
Labels*; see Citation below). It is organism-, phenotype-, and schema-agnostic:
given a continuous phenotype, a cohort/batch label per sample, and any genomic
distance that can flag near-clonal pairs, it estimates one calibration offset
per cohort and propagates calibration uncertainty into per-cohort credible
intervals and isolate-level probabilistic labels.

Zero external dependencies; builds offline.

## Idea in one line

Two near-clonal isolates from different cohorts share (approximately) the same
biological phenotype, so the difference in their measured values estimates the
difference in cohort/protocol offset — biology cancels. Offsets are fitted on a
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
`labC` ~5× higher than `labA`. `phenocal` recovers exactly that:

```
labA   1.00x  [1.00, 1.00]   (anchor)
labB   2.00x  [1.07, 3.71]
labC   4.99x  [2.81, 9.15]
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
  --anchor Ivanova          # optional; default = largest cohort
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
| `--seed <S>` | RNG seed | 20260603 |
| `--dashboard <FILE>` | also write the interactive HTML report | off |

Run `phenocal --help` for the same list.

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

## Algorithm (per cohort pair, per ruler)

1. Collect cross-cohort pairs with distance `d` and signed log-phenotype
   difference `s = y_a − y_b` (oriented `a` over `b`); keep smallest-`d` instance
   per isolate pair.
2. Adaptive `τ`: smallest distance bound with ≥ `min-support` pairs **and**
   median drift ≤ `max-drift` dilutions from the `d_min` value.
3. Edge value `Δ = median(s | d ≤ τ)`; `SE = log10(2)·sqrt((1+τ)/n)`; `w = 1/SE²`.
4. Weighted least squares on the comparison graph for cohort offsets, anchor δ=0.
5. Bootstrap (sample `Δ* ~ N(Δ, SE²)`, refit) → 95% credible intervals.
6. Per threshold `c`: `P(T_i(c)) = mean over draws of [ MIC_i / 10^δ ≥ c ]`.

The offsets are identifiable up to an additive constant (gauge); the anchor is a
coordinate convention and does not affect harmonised labels (Lemma 1 of the
paper). A connected graph is required.

## Validation

On the *Listeria monocytogenes* BC dataset (2099 isolates, 5 cohorts, cgMLST
pairs), `phenocal` reproduces the paper's Table 2: Kragh 0.92×, Palma ~1.17×,
Cooper ~5.0×, He ~11.4× (anchor Ivanova), weighted residual RMSE ≈ 0.10 log10.

## Citation

de Ruvo A†, Castelli P†, Di Pasquale A, Radomski N.
*Genome-Anchored Probabilistic MIC Labels: Uncertainty-Aware Phenotype
Construction for Cross-Cohort Machine Learning* (in preparation).
† equal contribution (co-first authors); either may list their name first on their CV.

## License

MIT.
