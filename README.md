# dna-repeat v1.0.0

## Description
This simple CLI identifies repetitive elements in DNA sequences from a FASTA file. The presence of such repeats poses significant challenges in traditional gene synthesis, as it can cause L and m oligonucleotides to misanneal, resulting in incorrect assembly products. Identifying these elements beforehand can greatly aid in optimizing oligo design and assembly steps, thereby increasing the likelihood of successful synthesis.

The tool reports the record ID, the repeats as a pair of k-mer strings, their positions in the sequence, and the number of mismatches between them. Results are written in CSV format to either stdout (i.e. the terminal) or to an output.csv file.

## Dependencies
- Linux/WSL2
- Python (≥ v3.12)
- Biopython (≥ v1.85)
- Pandas (≥ 2.3.3)
- tqdm (≥ 4.67.1)
- pytest
- UV Python package
  
## Installation
```bash
uv venv && uv sync
uv pip install -e .
```

## Usage
```
uv run dna-repeat INPUT.fasta [-k/--length INT] [-m/--mismatches INT] [-o/--output DIR] [-i/--inverted-only | -d/--direct-only]
```

- INPUT.fasta – input FASTA / Multi-FASTA filepath
- `-o`/`--output` – write results to an output.csv file in the current directory. Optional: path to desired destination directory for the CSV file (default: write results to terminal)
- `-k`/`--length` – repeat (k-mer) length (default: 20; allowed range 4–30)
- `-m`/`--mismatches` – allowed mismatches (default: 0; must be ≤ length/2)
- `-i`/`--inverted-only` - only look for inverted-repeats (default: look for both direct/inverted)
- `-d`/`--direct-only` - only look for direct-repeats (default: look for both direct/inverted)

## Examples:
- Find direct + inverted repeats of 23 bp in data/example.fas, allowing for up to 3 mismatches. Write to output.csv in the data directory:

```bash
$ uv run dna-repeat data/example.fas -k 23 -m 3 -o data
```

- Omit the `-o data` in order to write the results to the terminal instead of a CSV file:
```bash
$ uv run dna-repeat data/example.fas -k 23 -m 3
```

- Example showing output:
```bash
$ uv run dna-repeat data/example.fas -k 15 -m 2
```
```
Searching for repeats: 100%|████████████████████████████████████████████████████████| 4/4 [00:00<00:00,  4.86sequence/s]

Found repeats in data/example.fas:

record_id,query_start,query_end,subject_start,subject_end,query_seq,subject_seq,mismatches,kmer_length,orientation
gi|304112|gb|L10209.1|ATHTAGIA,89,103,756,770,TGATGAGCTCTTCTT,TGATGTGTTCTTCTT,2,15,direct
gi|304112|gb|L10209.1|ATHTAGIA,96,110,99,113,CTCTTCTTCTTCTAC,TTCTTCTTCTACTAC,2,15,direct
gi|304112|gb|L10209.1|ATHTAGIA,97,111,100,114,TCTTCTTCTTCTACT,TCTTCTTCTACTACT,1,15,direct
gi|304112|gb|L10209.1|ATHTAGIA,98,112,101,115,CTTCTTCTTCTACTA,CTTCTTCTACTACTC,2,15,direct
gi|304112|gb|L10209.1|ATHTAGIA,1095,1109,1466,1480,GAGACAACAGACCCT,AGGGTTTGTTGTCAC,2,15,inverted
gi|1019924|gb|M55554.1|ATHAGL6A,35,49,56,70,ATTCTTGAAAAAAAG,ATTCTTGAGAAGAAG,2,15,direct
gi|1019924|gb|M55554.1|ATHAGL6A,59,73,238,252,CTTGAGAAGAAGAAG,CATCATCTTCTCAAG,2,15,inverted
gi|1019924|gb|M55554.1|ATHAGL6A,60,74,237,251,TTGAGAAGAAGAAGA,TCATCATCTTCTCAA,2,15,inverted
gi|1019924|gb|M55554.1|ATHAGL6A,192,206,197,211,AGAAAGCTTATGAGC,GCTTATGAGCTTTCT,2,15,inverted

No repeats were found in the following sequence(s):
 gi|166595|gb|M55553.1|ATHAGL5A
```

## Future Development
- Switch from $O(n^2)$ brute force searching algorithm to something more efficient (such as **block-index filtering (pigeonhole)** and/or **k-mer hashing**) to improve time complexity.
- **Stretch goal:** visualizations of results (repeat map, dot plot, etc.)
## Release notes
See [CHANGELOG](CHANGELOG.md) for release notes.
## License
MIT - See [LICENSE](LICENSE)