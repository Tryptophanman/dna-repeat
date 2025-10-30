# dna-repeat

## Description
This simple CLI identifies repetitive elements in DNA sequences from a FASTA file. The presence of such repeats poses significant challenges in traditional gene synthesis, as it can cause L and m oligonucleotides to misanneal, resulting in incorrect assembly products. Identifying these elements beforehand can greatly aid in optimizing oligo design and assembly steps, thereby increasing the likelihood of successful synthesis.

In its current form, the tool reports the record ID, the repeats as a pair of k-mer strings, their positions in the sequence, and the number of mismatches between them. Results are written in CSV format to either stdout (i.e. the terminal) or to an output.csv file.

## Dependencies
- Linux/WSL2
- Python (≥ v3.12)
- Biopython (≥ v1.85)
- Pandas (≥ 2.3.3)
- pytest
- UV Python package installer/resolver
  
## Installation
```bash
uv venv && uv sync
uv pip install -e .
```

## Usage
```
dna-repeat [-k/--length INT] [-m/--mismatches INT] [-o/--output DIR] INPUT.fasta
```

- INPUT.fasta – input FASTA / Multi-FASTA filepath
- `-o`/`--output` – write results to output.csv file, into the specifed directory (otherwise, write to terminal (stdout))
- `-k`/`--length` – repeat (k-mer) length (default: 20; allowed range 4–30)
- `-m`/`--mismatches` – allowed mismatches (default: 0; must be ≤ length/2)

Examples:
- Write results to output.csv file in ./data:
```bash
$ dna-repeat -k 23 -m 3 -o data data/example.fas
```

- remove `-o data` to write the results to stdout:
```bash
$ dna-repeat -k 23 -m 3 data/example.fas
```

- Output to terminal (detailed example)
```bash
$ dna-repeat -k 14 -m 2 data/example.fas
record_id,query_start,subject_start,query_end,subject_end,query_seq,subject_seq,mismatches,kmer_length
gi|304112|gb|L10209.1|ATHTAGIA,79,1273,92,1286,CAATGGAGATGATG,CAAGGCAGATGATG,2,14
gi|304112|gb|L10209.1|ATHTAGIA,88,755,101,768,TGATGAGCTCTTCT,TGATGTGTTCTTCT,2,14
gi|304112|gb|L10209.1|ATHTAGIA,89,756,102,769,GATGAGCTCTTCTT,GATGTGTTCTTCTT,2,14
gi|304112|gb|L10209.1|ATHTAGIA,95,98,108,111,CTCTTCTTCTTCTA,TTCTTCTTCTACTA,2,14
gi|304112|gb|L10209.1|ATHTAGIA,96,99,109,112,TCTTCTTCTTCTAC,TCTTCTTCTACTAC,1,14
gi|304112|gb|L10209.1|ATHTAGIA,97,100,110,113,CTTCTTCTTCTACT,CTTCTTCTACTACT,1,14
gi|304112|gb|L10209.1|ATHTAGIA,98,101,111,114,TTCTTCTTCTACTA,TTCTTCTACTACTC,2,14
gi|304112|gb|L10209.1|ATHTAGIA,507,556,520,569,CAGCAGGGATTATG,CAGCTGGGAATATG,2,14
```

## Future Development
- Switch from $O(n^2)$ brute force searching algorithm to something more efficient (such as **block-index filtering (pigeonhole)** and/or **k-mer hashing**) to improve time complexity.
- Add **inverted repeats** detection
- Report **repeat positions**
- Implement **repeat block merging** (chain overlapping k-mer hits into longer repeat blocks)
- **Stretch goal:** visualizations of results (repeat map, dot plot, etc.)
## Release notes
See [CHANGELOG](CHANGELOG.md) for release notes.
## License
MIT - See [LICENSE](LICENSE)