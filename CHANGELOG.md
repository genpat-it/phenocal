# Changelog

All notable changes to phenocal. Every release gives byte-identical results
to the previous one on the same input unless stated otherwise.

## [0.1.2] - 2026-10-27

Distance matrices as input, and `cargo install phenocal`.

### Added
- **Distance-matrix input.** `--pairs` now also accepts a square distance
  matrix, so the output of the usual allele-distance tools can be passed as it
  is, for cgMLST and wgMLST alike:
  - chewBBACA: `distance_matrix.tsv` (upper-triangular) and
    `distance_matrix_symmetric.tsv`;
  - ReporTree (INSA): `<prefix>_dist_hamming.tsv`;
  - cgmlst-dists: full matrix, and `--matrix-format lower-tri`/`upper-tri`;
  - cgDist and any tool writing the same TSV layout (first header cell is a
    label, then the sample names in the header and at the start of each row,
    in the same order).

  The layout is detected from the header: a file with `sample_i`, `sample_j`
  and `distance` columns is read as a pair table, anything else as a matrix.
  `--pairs-format long|matrix` forces one or the other.
- **Safe reading of triangular matrices.** cgmlst-dists fills the unused half
  of a triangular matrix with zeros, which a naive reader would take as
  identical genomes. phenocal finds the half that holds the distances (upper,
  lower, or both when they mirror each other), rejects a matrix whose halves
  disagree, and reads each pair once.
- **Sample-name matching.** Matrix names are matched to the phenotype `sample`
  ids exactly, without the FASTA extension (`S1.fasta` → `S1`), or by unique id
  prefix (`S1_shovill.fasta` → `S1`, the assembler suffix chewBBACA keeps from
  the file name). Ambiguous names are an error; the run reports how many names
  matched each way.
- The matrix is streamed in two passes with memory linear in the number of
  samples, so large wgMLST matrices are fine.
- `examples/toy_distance_matrix.tsv`: the toy as a chewBBACA-style matrix; a
  test checks it gives exactly the outputs of `examples/toy_pairs.tsv`.
- Published on crates.io: `cargo install phenocal`.

### Verified
On 2,099 *Listeria monocytogenes* genomes (cgMLST 1,748 loci and wgMLST 6,785
loci; 2.2 million pairs) every matrix layout above gives outputs byte-identical
to the equivalent pair table, and pair tables give outputs byte-identical to
0.1.1.

## [0.1.1] - 2026-10-13

Packaging: phenocal can be installed without a Rust toolchain. Results are
byte-identical to 0.1.0 and to the `ais2026` tag.

### Added
- **One-line installer**:
  `curl -fsSL https://raw.githubusercontent.com/genpat-it/phenocal/main/install.sh | sh`.
  It picks the binary for the platform, checks its SHA-256 and installs it in
  `~/.local/bin` (`PHENOCAL_INSTALL_DIR` and `PHENOCAL_VERSION` change the
  directory and the release).
- **Prebuilt static binaries** for Linux (x86_64, arm64; musl), macOS (Intel,
  Apple Silicon) and Windows (x86_64), attached to each GitHub Release with
  their SHA-256 checksums.
- **Docker image** `ghcr.io/genpat-it/phenocal` for linux/amd64 and
  linux/arm64, usable from Nextflow and other workflow managers; it ships the
  toy example in `/opt/phenocal/examples`.
- `phenocal --version` / `-V`.
- CI on Linux, macOS and Windows, a toy run inside the Docker image, and
  shellcheck of the shell scripts.
- `CITATION.cff` (GitHub "Cite this repository"), this changelog, and README
  badges and installation section.

## `ais2026` tag - 2026-07-29

The code exactly as used in the AIS 2026 paper (*PhenoCal: Natural-Control
Label Calibration for Cross-Cohort Machine Learning*, Springer LNCS). This tag
is frozen and will not move. The code is that of 0.1.0; only documentation
changed:
- README: citation updated to the accepted AIS 2026 paper, the Validation
  section describes the AIS 2026 case study (benzalkonium-chloride MICs of
  *L. monocytogenes* across five laboratory cohorts), corresponding e-mail fixed.
- LICENSE: single-line MIT copyright holder, so GitHub detects the licence.

## [0.1.0] - 2026-06-10

First public version: genome-anchored calibration of a continuous phenotype
measured in several cohorts, with no dependencies.

### Method
- Near-clonal cross-cohort pairs are used as natural controls: the phenotype
  difference between two near-identical genomes measured in different cohorts
  estimates the difference of the cohort offsets.
- Per cohort pair, an adaptive relatedness bound τ (the smallest distance with
  `--min-support` pairs and median drift within `--max-drift`), an edge value
  (median log-ratio) and its standard error.
- Cohort offsets by weighted least squares on the comparison graph, with one
  anchor cohort (`--anchor`, default the largest).
- Bootstrap credible intervals on the offsets (`--bootstrap`, `--robust` for
  Student-t perturbations), propagated to per-isolate probabilities of
  exceeding each threshold (`--thresholds`).
- Per-cohort measurement resolution, fixed or estimated from each cohort's
  grid (`--sigma`), drift-tolerance units (`--drift-scale`) and biological
  drift variance (`--lambda`).
- Data-driven sensitive/tolerant cutoff on the harmonised scale: KDE antimode
  and GMM crossover with a bootstrap interval (`--no-cutoff` skips it).
- Cluster-effective inference (`--clusters`): cluster-weighted point estimate
  and an interval that resamples clonal clusters globally over the graph.

### Outputs
- `offsets`, `sigma`, `edges`, `harmonised`, `labels` and `cutoff` TSV tables.
- `--dashboard`: a self-contained interactive HTML report (comparison graph,
  per-edge twin pairs, forest plot, before/after harmonisation, per-isolate
  labels, cutoff panel).
- `--trace`: a step-by-step Markdown trace of the run with the real numbers
  and formulas, which renders on GitHub.
- `--quick`: offsets, sigma and edges only, for large benchmarks.

### Quality
- Test suite (offset recovery, anchor invariance, edge selection, cutoff,
  cluster-effective recovery) and CI with rustfmt and clippy.

[0.1.2]: https://github.com/genpat-it/phenocal/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/genpat-it/phenocal/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/genpat-it/phenocal/releases/tag/v0.1.0
