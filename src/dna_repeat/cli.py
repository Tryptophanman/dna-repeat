# src/dna_repeat/cli.py
    # Top-level functions for dna-repeat CLI

import sys
import argparse
import csv
from dna_repeat.core import read_fasta, find_repeats
from pathlib import Path

VERSION: str = 'v0.1.0'     # cross-check with pyproject.toml

class MyArgs(argparse.Namespace):
    '''Provides static type hints for the arguments defined in init_argparser()'''
    input_filepath: Path
    output_directory: Path
    rep_length: int
    allowed_mismatches: int
    f_flag: bool

def init_argparser() -> argparse.ArgumentParser:
    '''Initializes argument parser for CLI'''
    parser = argparse.ArgumentParser(prog='dna-repeat', description=f'dna-repeat {VERSION} - Finds repeats in a DNA sequence')
    parser.add_argument("input_filepath", help='filepath to input file (FASTA)')
    parser.add_argument("output_directory", help='(optional) directory for output file [output.csv]. Default: cwd', nargs='?', default='.')
    parser.add_argument("-l", "--length", help='repeat length in bp (min 4, max 30)', type=int, default=20, dest='rep_length')
    parser.add_argument("-m", "--mismatches", help='number of bp mismatches allowed', type=int, default=0, dest='allowed_mismatches')
    parser.add_argument("-f", "--file", help='make output file (default: stdout)', action='store_true', dest='f_flag')
    return parser

def main(argv: str | None = None) -> int:
    parser = init_argparser()
    custom_namespace = MyArgs()
    args: MyArgs = parser.parse_args(argv, namespace=custom_namespace)

    input_filepath = Path(args.input_filepath).resolve()
    output_directory = Path(args.output_directory).resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    output_filepath: Path = output_directory / 'output.csv'
    rep_length: int = args.rep_length
    allowed_mismatches: int = args.allowed_mismatches

    if not input_filepath.is_file():
        print(f'Input file not found: {input_filepath}', file=sys.stderr)
        return 1
    if not 4 <= rep_length <= 30:
        print('repeat length (--l, --length) must be in the range 4-30', file=sys.stderr)
        return 1
    if allowed_mismatches > rep_length / 2:
        print('m must be <= rep_length / 2', file=sys.stderr)
        return 1       

    writer = csv.writer(open(output_filepath, 'w', newline='') if args.f_flag else sys.stdout)
    writer.writerow(['record_id', 'seq1', 'seq2', 'mismatches'])

    for rec_id, seq in read_fasta(input_filepath):
        hits = find_repeats(seq, rep_length, allowed_mismatches)
        if hits:
            for seq_a, seq_b, mis in hits:
                writer.writerow([rec_id, seq_a, seq_b, mis])
    return 0


if __name__ == "__main__":
    sys.exit(main())

