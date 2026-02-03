#!/usr/bin/env bash
# Generate review_outputs/ with all option combinations for manual quality check.
# Run from project root: bash scripts/generate_review_outputs.sh

set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
X="python -m xmoltoppm.cli.main"
TRAJ="traj.xyz"
OUT="review_outputs"
mkdir -p "$OUT"
mkdir -p "$OUT/movie_ft0"
mkdir -p "$OUT/movie_ft1"
mkdir -p "$OUT/movie_ft2"

# Circle style: -cc 0 and -cc 1
$X -i "$TRAJ" -s 500 -cc 0 -il 3 -o "$OUT/single_cc0_il3.ppm"
$X -i "$TRAJ" -s 500 -cc 1 -il 3 -o "$OUT/single_cc1_il3.ppm"

# Line style: -il 0, 1, 2, 3, 4 (single frame, same cc)
for il in 0 1 2 3 4; do
  $X -i "$TRAJ" -s 500 -cc 1 -il $il -o "$OUT/single_cc1_il${il}.ppm"
done

# Output format: PPM vs PNG (same options)
$X -i "$TRAJ" -s 500 -cc 1 -il 3 -o "$OUT/single_ppm.ppm"
$X -i "$TRAJ" -s 500 -cc 1 -il 3 -o "$OUT/single_png.png" --png

# All-frames: -af 0 on traj.xyz (all frames), -ft 0 -> 0001.ppm, 0002.ppm, ...
$X -i "$TRAJ" -s 500 -cc 1 -il 3 -af 0 -ft 0 -o "$OUT/movie_ft0/frame.ppm"

# Filename type -ft 1 -> 1.ppm, 2.ppm, ...
$X -i "$TRAJ" -s 500 -cc 1 -il 3 -af 0 -ft 1 -o "$OUT/movie_ft1/frame.ppm"

# Filename type -ft 2 -> molname.ppm (overwrites per frame; last frame wins)
$X -i "$TRAJ" -s 500 -cc 1 -il 3 -af 0 -ft 2 -o "$OUT/movie_ft2/frame.ppm"

# Shadow: -sh 0.6 -il 4 -sf 20 -sd yz
$X -i "$TRAJ" -s 500 -cc 1 -il 4 -sh 0.6 -sf 20 -sd yz -o "$OUT/single_shadow.ppm"

# Boundary lines: -bl 0 and -bl 1
$X -i "$TRAJ" -s 500 -cc 1 -il 3 -bl 0 -o "$OUT/single_bl0.ppm"
$X -i "$TRAJ" -s 500 -cc 1 -il 3 -bl 1 -o "$OUT/single_bl1.ppm"

# Size: small canvas
$X -i "$TRAJ" -s 200 -cc 1 -il 3 -o "$OUT/single_size200.ppm"

# Border: larger border
$X -i "$TRAJ" -s 500 -cc 1 -il 3 -bs 1.0 -o "$OUT/single_bs1.ppm"

echo "Generated all review outputs in $OUT/"
