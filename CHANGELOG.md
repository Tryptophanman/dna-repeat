# CHANGELOG.md

See also: [README.md](README.md)

All notable changes to this project will be documented in this file.


## [Unreleased]
- A number of fields (columns) have been added to the results:
  - record_id: name of the sequence
  - query_start: start position of the first repeat
  - subject_start: start position of the second repeat
  - query_end: end position of the first repeat
  - subject_end: end position of the second repeat
  - query_seq: first repeat sequence
  - subject_seq: second repeat sequence
  - mismatches: number of mismatches between the repeats
  - kmer_length: length of the repeat in bp
- Added detection of inverted repeats
- Added option to choose to look for direct-only or inverted-only repeats
- Added faster repeat search algorithm "2-bit" by encoding k-mers as 2*k bits and doing bitwise comparisons (AI-assisted!)
- Added sequence validation/checking function and some error handling of exceptions (Empty sequence, invalid characters in sequence, k-mer length longer than sequence itself)
- Added progress bar

## [0.2.0] - 2025-10-27
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