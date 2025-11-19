# src/dna_repeat/cli.py
# Top-level functions for dna-repeat CLI

import sys
import argparse
from dna_repeat import __version__
from dna_repeat.core import iter_fasta, RepeatHit, clean_and_check, count_fasta_seqs
from dna_repeat.ai import find_repeats_2bit, find_invert_repeats_2bit
from dna_repeat.error import (
    InvalidFASTAError,
    EmptySequenceError,
    InvalidSequenceError,
    InvalidKmerError,
)
from pathlib import Path
import pandas as pd
from tqdm import tqdm


def init_argparser() -> argparse.ArgumentParser:
    """Initializes argument parser for CLI"""
    parser = argparse.ArgumentParser(
        prog="dna-repeat",
        description=f"dna-repeat {__version__} - Finds repeats in a DNA sequence",
    )
    parser.add_argument("input_filepath", help="filepath to input file (FASTA)")
    parser.add_argument(
        "-o",
        "--output",
        help="desired directory for output file (output.csv). Default: cwd",
        nargs="?",
        const=".",
        default=None,
        dest="output_directory",
    )
    parser.add_argument(
        "-k",
        "--length",
        help="repeat length in bp (min 4, max 30)",
        type=int,
        default=20,
        dest="kmer_length",
    )
    parser.add_argument(
        "-m",
        "--mismatches",
        help="number of bp mismatches allowed",
        type=int,
        default=0,
        dest="allowed_mismatches",
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-i",
        "--inverted-only",
        help="only look for inverted repeats",
        action="store_true",
        dest="inverted_only",
    )
    group.add_argument(
        "-d",
        "--direct-only",
        help="only look for direct repeats",
        action="store_true",
        dest="direct_only",
    )
    return parser


def main(argv: str | None = None) -> int:
    parser = init_argparser()
    args = parser.parse_args(argv)

    input_filepath = Path(args.input_filepath).resolve()
    if args.output_directory:
        output_directory = Path(args.output_directory).resolve()
        output_directory.mkdir(parents=True, exist_ok=True)
        output_filepath: Path | None = output_directory / "output.csv"
    else:
        output_filepath = None
    kmer_length: int = args.kmer_length
    allowed_mismatches: int = args.allowed_mismatches

    do_direct = not args.inverted_only
    do_inverted = not args.direct_only

    if not input_filepath.is_file():
        print(f"Input file not found: {input_filepath}", file=sys.stderr)
        return 1
    try:
        number_of_seqs = count_fasta_seqs(input_filepath)
    except InvalidFASTAError as e:
        print(f"InvalidFASTAError: {e.message}", file=sys.stderr)
        if e.details:
            print(f"Details: {e.details}", file=sys.stderr)
        return e.exit_code
    if not 4 <= kmer_length <= 30:
        print("repeat length (-k, --length) must be in the range 4-30", file=sys.stderr)
        return 1
    if allowed_mismatches > kmer_length / 2:
        print("m must be <= kmer_length / 2", file=sys.stderr)
        return 1

    results = []
    errors: list[str] = []
    no_hits: list[str] = []

    try:
        for rec_id, seq in tqdm(
            iter_fasta(input_filepath),
            desc="Searching for repeats",
            total=number_of_seqs,
            unit="sequence",
        ):
            try:
                rec_id, seq = clean_and_check(rec_id, seq, kmer_length)
                hits_found = False
                if do_direct:
                    repeats: list[RepeatHit] = find_repeats_2bit(
                        rec_id=rec_id,
                        seq=seq,
                        kmer_length=kmer_length,
                        allowed_mismatches=allowed_mismatches,
                    )
                    if repeats:
                        results.extend(repeats)
                        hits_found = True
                if do_inverted:
                    repeats: list[RepeatHit] = find_invert_repeats_2bit(
                        rec_id=rec_id,
                        seq=seq,
                        kmer_length=kmer_length,
                        allowed_mismatches=allowed_mismatches,
                    )
                    if repeats:
                        results.extend(repeats)
                        hits_found = True
                if not hits_found:
                    no_hits.append(rec_id)
            except EmptySequenceError as e:
                errors.append(f"{e}")
                continue
            except InvalidSequenceError as e:
                errors.append(f"{e}")
                continue
            except InvalidKmerError as e:
                errors.append(f"{e}")
                continue
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        return 255
    if results:
        df = pd.DataFrame(results)
        if output_filepath:
            df.to_csv(output_filepath, index=False)
            print(f"\nResults have been written to {output_filepath}\n")
        else:
            print(f"\nFound repeats in {input_filepath}:\n")
            df.to_csv(sys.stdout, index=False)
    if no_hits:
        print("\nNo repeats were found in the following sequence(s):")
        for item in no_hits:
            print(f" {item}")
    if errors:
        print("\nThe following sequence(s) returned errors:")
        for item in errors:
            print(f" {item}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
