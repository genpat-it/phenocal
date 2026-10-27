# Changelog

## [0.1.2] - 2026-10-27

### Added
- `--pairs` also accepts a square distance matrix: the chewBBACA allelic
  distance matrix (cgMLST or wgMLST; `distance_matrix.tsv` upper-triangular or
  `distance_matrix_symmetric.tsv`) and the output of cgmlst-dists and cgDist,
  including triangular matrices with a zero-filled half. The layout is
  detected from the header (`--pairs-format` forces it). Matrix names are
  matched to phenotype ids exactly, without the FASTA extension or by unique id
  prefix (`S1_shovill.fasta` → `S1`).
- `examples/toy_distance_matrix.tsv`.
- Published on crates.io: `cargo install phenocal`.

On 2,099 *L. monocytogenes* genomes (cgMLST and wgMLST), every matrix layout
gives outputs byte-identical to the equivalent pair table, and pair tables give
outputs byte-identical to 0.1.1.

## [0.1.1] - 2026-10-13

Packaging release: results are byte-identical to 0.1.0 and to the `ais2026` tag.

### Added
- Prebuilt static binaries for Linux (x86_64, arm64; musl), macOS (Intel,
  Apple Silicon) and Windows (x86_64), attached to each GitHub Release with
  SHA-256 checksums.
- One-line installer: `curl -fsSL https://raw.githubusercontent.com/genpat-it/phenocal/main/install.sh | sh`.
- Docker image `ghcr.io/genpat-it/phenocal` (linux/amd64, linux/arm64),
  usable from Nextflow.
- `phenocal --version` / `-V`.
- CI on Linux, macOS and Windows, plus Docker and installer checks.
- `CITATION.cff`.

## [0.1.0] - 2026-06-10

First public version. The code used in the AIS 2026 paper is frozen at the
`ais2026` tag.
