# Changelog

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
