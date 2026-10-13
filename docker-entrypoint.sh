#!/bin/sh
# Entry point of the phenocal image.
#   docker run IMAGE --phenotypes ... --pairs ...   runs phenocal with these arguments
#   docker run IMAGE /bin/bash -c '...'             runs a command (Nextflow, etc.)
if [ "$#" -eq 0 ]; then
    exec phenocal --help
fi
case "$1" in
    -*) exec phenocal "$@" ;;
esac
if command -v "$1" >/dev/null 2>&1; then
    exec "$@"
fi
exec phenocal "$@"
