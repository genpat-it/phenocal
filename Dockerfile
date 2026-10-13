# phenocal has no dependencies: build one static (musl) binary and copy it
# into a slim runtime image.
FROM rust:1-alpine AS builder
WORKDIR /app
COPY Cargo.toml Cargo.lock ./
COPY src/ ./src/
RUN cargo build --release --locked

FROM debian:bookworm-slim
# procps: `ps`, used by Nextflow to collect task metrics
RUN apt-get update && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/* \
    && useradd -r -s /bin/false phenocal \
    && mkdir -p /data && chown phenocal /data
COPY --from=builder /app/target/release/phenocal /usr/local/bin/phenocal
COPY examples/ /opt/phenocal/examples/
COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
USER phenocal
WORKDIR /data
# Options run phenocal (`docker run IMAGE --phenotypes ...`); a command runs
# that command (`/bin/bash` from Nextflow and other workflow managers).
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

LABEL org.opencontainers.image.title="phenocal" \
      org.opencontainers.image.description="Genome-anchored calibration of a continuous phenotype across cohorts/batches" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.source="https://github.com/genpat-it/phenocal"
