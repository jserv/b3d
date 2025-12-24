#!/bin/bash
# Generate compressed PNG snapshots from B3D examples
#
# Usage: ./scripts/gen-snapshots.sh <example1> [example2] ...
#
# Runs each SDL2 example with --snapshot= to generate PNG files
# in examples/ directory, then compresses them using ffmpeg.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
EXAMPLES_DIR="$PROJECT_DIR/examples"

if [ $# -eq 0 ]; then
    echo "Usage: $0 <example1> [example2] ..." >&2
    exit 1
fi

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

info() { printf "${GREEN}[INFO]${NC} %s\n" "$*"; }
warn() { printf "${YELLOW}[WARN]${NC} %s\n" "$*"; }
error() { printf "${RED}[ERROR]${NC} %s\n" "$*" >&2; }

# Compress PNG using ffmpeg
# @path: PNG file path
# Returns: compression ratio (0-100), or exits non-zero on failure
compress_png() {
    local path="$1"
    local orig_size new_size tmp_file ratio

    orig_size=$(stat -f%z "$path" 2>/dev/null || stat -c%s "$path" 2>/dev/null)

    # Create temp file with .png extension for ffmpeg format detection
    tmp_file=$(mktemp "${path%.png}.XXXXXX.png")
    trap 'rm -f "$tmp_file"' RETURN

    if ! ffmpeg -y -loglevel error -i "$path" -compression_level 9 "$tmp_file"; then
        return 1
    fi
    mv -f "$tmp_file" "$path"

    new_size=$(stat -f%z "$path" 2>/dev/null || stat -c%s "$path" 2>/dev/null)

    if [ "$orig_size" -gt 0 ]; then
        ratio=$((100 - (new_size * 100 / orig_size)))
    else
        ratio=0
    fi
    echo "$ratio"
}

# Main
main() {
    info "B3D Snapshot Generator"
    info "Output directory: $EXAMPLES_DIR"

    # Check for ffmpeg
    if ! command -v ffmpeg &>/dev/null; then
        error "ffmpeg not found. Install with your package manager:"
        error "  macOS: brew install ffmpeg"
        error "  Debian/Ubuntu: sudo apt install ffmpeg"
        error "  Fedora: sudo dnf install ffmpeg"
        exit 1
    fi

    # Normalize arguments: strip paths, use basename only
    local examples=()
    for arg in "$@"; do
        examples+=("$(basename "$arg" .c)")
    done

    # Build requested examples (show output on failure)
    info "Building examples..."
    local build_log
    build_log=$(mktemp)
    trap 'rm -f "$build_log"' RETURN
    if ! make -C "$PROJECT_DIR" "${examples[@]}" >"$build_log" 2>&1; then
        error "Failed to build examples:"
        cat "$build_log" >&2
        return 1
    fi

    # Generate snapshots
    local success=0
    local failed=0
    local total_orig=0
    local total_new=0

    for example in "${examples[@]}"; do
        local exe="$PROJECT_DIR/$example"
        local png="$EXAMPLES_DIR/${example}.png"

        if [ ! -x "$exe" ]; then
            warn "Skipping $example (not built)"
            ((failed++)) || true
            continue
        fi

        printf "  %-12s " "$example"

        # Run example with snapshot
        if ! "$exe" --snapshot="$png" >/dev/null 2>&1; then
            echo -e "${RED}FAILED${NC}"
            ((failed++)) || true
            continue
        fi

        if [ ! -f "$png" ]; then
            echo -e "${RED}NO OUTPUT${NC}"
            ((failed++)) || true
            continue
        fi

        local orig_size
        orig_size=$(stat -f%z "$png" 2>/dev/null || stat -c%s "$png" 2>/dev/null)
        total_orig=$((total_orig + orig_size))

        # Compress with ffmpeg (continue on failure)
        local ratio
        if ! ratio=$(compress_png "$png"); then
            echo -e "${YELLOW}COMPRESS FAILED${NC} (%d KB, kept original)" "$((orig_size / 1024))"
            total_new=$((total_new + orig_size))
            ((failed++)) || true
            continue
        fi

        local new_size
        new_size=$(stat -f%z "$png" 2>/dev/null || stat -c%s "$png" 2>/dev/null)
        total_new=$((total_new + new_size))
        printf "${GREEN}OK${NC} (%d KB → %d KB, -%d%%)\n" \
               "$((orig_size / 1024))" "$((new_size / 1024))" "$ratio"

        ((success++)) || true
    done

    echo ""
    info "Summary: $success succeeded, $failed failed"
    if [ $total_orig -gt 0 ]; then
        local total_ratio=$((100 - (total_new * 100 / total_orig)))
        info "Total size: $((total_orig / 1024)) KB → $((total_new / 1024)) KB (-${total_ratio}%)"
    fi
    info "Output: $EXAMPLES_DIR/"

    [ $failed -eq 0 ]
}

main "$@"
