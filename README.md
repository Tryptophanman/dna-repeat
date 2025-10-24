# dna-repeat v0.1.0

## Description
This simple CLI identifies repetitive elements in DNA sequences from a FASTA file. The presence of such repeats poses significant challenges in traditional gene synthesis, as it can cause L and m oligonucleotides to misanneal, resulting in incorrect assembly products. Identifying these elements beforehand can greatly aid in optimizing oligo design and assembly steps, thereby increasing the likelihood of successful synthesis.

In its current form (v0.1.0), the tool reports the record ID, the repeats as a pair of k-mer strings, and the number of mismatches between them. Results are written in CSV format to either stdout (i.e. the terminal) or to an output.csv file.

## Dependencies
- Linux/WSL2
- Python (≥ v3.12)
- Biopython (≥ v1.85)
- pytest
- uv Python package installer/resolver
  
## Installation
```bash
uv venv && uv sync
uv pip install -e .
```

## Usage
```
dna-repeat [-l/--length INT] [-m/--mismatches INT] [-f/--file] INPUT.fasta [OUTPUT_DIR]
```

- INPUT.fasta – input FASTA / Multi-FASTA filepath
- OUTPUT_DIR (optional) – directory to write output.csv into (default: current directory)
- `-l`/`--length` – repeat (k-mer) length (default: 20; allowed range 4–30)
- `-m`/`--mismatches` – allowed mismatches (default: 0; must be ≤ length/2)
- `-f`/`--file` – write CSV to OUTPUT_DIR/output.csv (otherwise, write to terminal (stdout))

Examples:
- Write results to output.csv file in ./data:
```bash
$ dna-repeat -l 23 -m 3 -f data/example.fas data
```

- remove `-f` to write the results to the terminal (stdout):
```bash
$ dna-repeat -l 23 -m 3 data/example.fas
```

- Output to terminal (detailed example)
```bash
$ dna-repeat -l 23 -m 5 fasfiles/seqs.fasta
record_id,seq1,seq2,mismatches
gi|304112|gb|L10209.1|ATHTAGIA,GCAAAAGTCGTTTGAGAAAGAAG,GCAACAGTCGTCTCAGCAAGCAG,5
gi|166595|gb|M55553.1|ATHAGL5A,GAAGATAGGGAGAGGGAAGATAG,GAAGATAGAGATAAAGAGGATAG,5
gi|166595|gb|M55553.1|ATHAGL5A,AAGATAGGGAGAGGGAAGATAGA,AAGATAGAGATAAAGAGGATAGA,5
gi|166595|gb|M55553.1|ATHAGL5A,AGATAGGGAGAGGGAAGATAGAG,AGATAGAGATAAAGAGGATAGAG,5
gi|166595|gb|M55553.1|ATHAGL5A,GATAGGGAGAGGGAAGATAGAGA,GATAGAGATAAAGAGGATAGAGA,5
```

## Future Development
- Switch from `O(n^2)` brute force searching algorithm to something more efficient (such as **block-index filtering (pigeonhole)** and/or **k-mer hashing**) to improve time complexity.
- Add **inverted repeats** detection
- Report **repeat positions**
- Implement **repeat block merging** (chain overlapping k-mer hits into longer repeat blocks)
- **Stretch goal:** visualizations of results (repeat map, dot plot, etc.)
## Release notes
See [CHANGELOG](CHANGELOG.md) for release notes.
## License
MIT - See [LICENSE](LICENSE)