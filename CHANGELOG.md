# CHANGELOG.md

See also: [README.md](README.md)

All notable changes to this project will be documented in this file.


## [Unreleased]

## [0.1.1] - 2025-10-27
- fix: arguments/flags cleaned up a bit. -o / --output flag added for output directory; omitting means output results to stdout. Replaced -l with -k for consistency (k-mer length)

## [0.1.0] - 2025-10-21
- Initial features:
  - Parse sequences from FASTA via Biopython
  - Output results to CSV or stdout
  - Brute-force repeat finder (nested loops, Hamming distance ≤ m)
  - Command-line interface using argparse
  - Basic pytest

- Known limitations:
  - $O(n^2)$ time complexity; slow for large inputs
  - Inverted repeats not detected
  - Outputs k-mer strings but not positions