#!/bin/sh
# Install the prebuilt phenocal binary from the GitHub releases.
#
#   curl -fsSL https://raw.githubusercontent.com/genpat-it/phenocal/main/install.sh | sh
#
# Environment variables:
#   PHENOCAL_VERSION       release tag to install (default: latest), e.g. v0.1.1
#   PHENOCAL_INSTALL_DIR   destination directory (default: $HOME/.local/bin)
set -eu

REPO="genpat-it/phenocal"
VERSION="${PHENOCAL_VERSION:-latest}"
INSTALL_DIR="${PHENOCAL_INSTALL_DIR:-$HOME/.local/bin}"

err() { echo "phenocal-install: $*" >&2; exit 1; }

case "$(uname -s)" in
    Linux) os="unknown-linux-musl" ;;
    Darwin) os="apple-darwin" ;;
    *) err "unsupported OS $(uname -s); see https://github.com/$REPO#installation" ;;
esac
case "$(uname -m)" in
    x86_64 | amd64) arch="x86_64" ;;
    aarch64 | arm64) arch="aarch64" ;;
    *) err "unsupported architecture $(uname -m); build from source with cargo" ;;
esac
asset="phenocal-${arch}-${os}.tar.gz"

if [ "$VERSION" = "latest" ]; then
    base="https://github.com/$REPO/releases/latest/download"
else
    base="https://github.com/$REPO/releases/download/$VERSION"
fi

if command -v curl >/dev/null 2>&1; then
    fetch() { curl -fsSL "$1" -o "$2"; }
elif command -v wget >/dev/null 2>&1; then
    fetch() { wget -q "$1" -O "$2"; }
else
    err "curl or wget is required"
fi

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT INT TERM

echo "Downloading $asset ($VERSION)..."
fetch "$base/$asset" "$tmp/$asset" || err "download failed: $base/$asset"
fetch "$base/$asset.sha256" "$tmp/$asset.sha256" || err "download failed: $base/$asset.sha256"

expected="$(cut -d ' ' -f 1 < "$tmp/$asset.sha256")"
if command -v sha256sum >/dev/null 2>&1; then
    actual="$(sha256sum "$tmp/$asset" | cut -d ' ' -f 1)"
elif command -v shasum >/dev/null 2>&1; then
    actual="$(shasum -a 256 "$tmp/$asset" | cut -d ' ' -f 1)"
else
    err "sha256sum or shasum is required to verify the download"
fi
[ "$expected" = "$actual" ] || err "checksum mismatch for $asset"

tar -xzf "$tmp/$asset" -C "$tmp"
mkdir -p "$INSTALL_DIR"
install -m 755 "$tmp/phenocal" "$INSTALL_DIR/phenocal" 2>/dev/null \
    || { cp "$tmp/phenocal" "$INSTALL_DIR/phenocal" && chmod 755 "$INSTALL_DIR/phenocal"; }

echo "Installed $("$INSTALL_DIR/phenocal" --version) to $INSTALL_DIR/phenocal"
case ":$PATH:" in
    *":$INSTALL_DIR:"*) ;;
    *) echo "Note: $INSTALL_DIR is not in your PATH; add it with:"
       echo "  export PATH=\"$INSTALL_DIR:\$PATH\"" ;;
esac
