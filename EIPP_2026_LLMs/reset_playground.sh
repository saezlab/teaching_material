#!/usr/bin/env bash
#
# Reset the playground/ folder to its starting state (empty directory).
# Students run this before re-doing exercises from scratch.
#
# Usage:
#   ./reset_playground.sh        (asks for confirmation)
#   ./reset_playground.sh -y     (skip confirmation)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PLAYGROUND="$SCRIPT_DIR/playground"

if [ ! -d "$PLAYGROUND" ]; then
    echo "Nothing to reset — $PLAYGROUND does not exist."
    exit 0
fi

# --- Confirmation prompt (unless -y flag is passed) ---
if [ "${1:-}" != "-y" ]; then
    echo "This will DELETE everything inside:"
    echo "  $PLAYGROUND"
    echo ""
    echo "Including .venv, .env, pyproject.toml, uv.lock, and all exercise files."
    echo ""
    read -rp "Are you sure? [y/N] " confirm
    if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
        echo "Aborted."
        exit 0
    fi
fi

# --- Nuke contents, keep the directory ---
echo "Resetting playground..."
rm -rf "$PLAYGROUND"
mkdir -p "$PLAYGROUND"

echo "Done. playground/ is now empty and ready for Exercise 0."
