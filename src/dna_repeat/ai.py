# src/dna_repeat/ai.py
    # Core code either written or inspired by AI. (written ones are marked by "[MADE BY AI]")
    ## Faster repeat-finding algorithm; converts k-mers to 2*k bits and does some bit tricks to check for repeats! About 5x faster than string slicing method

from dna_repeat.constants import BASE_TO_INT
from dna_repeat.core import RepeatHit, reverse_complement

def dna_to_int(seq: str) -> int: 
    '''[MADE BY AI] Executes DNA_seq-to-integer conversion; DNA string -> 2*k bits (int)'''

    # I asked AI to generate a faster way to search for dna repeats than just doing substring searches. It came up with this approach
    # that convertes the k-mer dna string to an integer (2-bit encoded) and then after doing some fancy bitshifting stuff for the
    # hamming distance calculation. So you're going from 8-bit string characters to 2-bit, saving 4x memory for the whole string.
    # Also the comparions of the two substrings should be much faster when going bitwise, rather than comparing string chars one by one.
    # Assuming the CPU handles 64 bit, I guess 32 bp would be the max k-mer length (at least to stay fast, without using RAM).
    # The comments here are written by me (not AI) to help me get what's going on in the code so I can understand.

    x = 0
    for base in seq:
        x = (x << 2) | BASE_TO_INT[base]    # (x << 2) = Left bitwise shift of 2 on integer x. This effectively adds 2 zeros to the right side
                                            # makes space for the incoming integer (2 bits).
                                            # Pipe (|) = Bitwise OR operator. Merges the 2 new bits with those shifted zeros:
                                            # -> ORing a bit with a 0 will always keep the original bit unchanged:
                                            # -> 0|0 = 0, 0|1 = 1, 1|0 = 1, 1|1 = 1 (orig bit is 0, so only changes if incoming is 1)
    return x

def encode_kmers(seq: str, kmer_length: int) -> list[int]:
    '''[MADE BY AI] Returns a list of all k-mers, 2-bit encoded. Uses a rolling window'''
    codes = [0] * (len(seq) - kmer_length + 1)
    mask: int = ((1 << (2 * kmer_length)) - 1)              # Mask of all 1's for keeping a k-mer length of 2*k bits 

    codes[0] = dna_to_int(seq[:kmer_length])                # Get "code" of first k-mer, as int
    for i in range(1, len(seq) - kmer_length + 1):          # Roll through the seq and get the rest of the k-mers and put into list
        val = BASE_TO_INT[seq[i + kmer_length - 1]]         # get int value of next base...
        codes[i] = ((codes[i - 1] << 2) | val) & mask       # ...and tack it onto the right side of the previous code ((codes[i - 1] << 2) | val)
                                                            # The two on the left side are trimmed by ANDing on the mask. This works because the mask
                                                            # is 2 bits shorter; there are effectively 2 0's on the left side. Hence when ANDing, the 
                                                            # result is going to be 2 0's on the left so...trimmed!
    return codes

def hamming_distance_2bit(code_a: int, code_b: int, kmer_length: int) -> int:
    '''[MADE BY AI] Calculates the hamming distance between to 2-bit-encoded sequence using some fancy bit tricks'''
    diff = code_a ^ code_b                              # Find bit mismatches using XOR (^). A mismatch will be either 01, 10, or 11
                                                        # ---
                                                        # The next part does some bit tricks to move 1's to the least significant bit position (LSB)
                                                        # and apply a mask of 01 01 01 by AND so in the end all mismatches types are 01
                                                        # and can be counted by bitcount. Wow this is pretty cool.
                                                        # ---
    m = diff | diff >> 1                                # Bit-shift diff to the right once and then OR it on diff. All mismatches get a 1 in LSB (01 or 11)
    lsb_bit_mask = ((1 << (2 * kmer_length)) - 1) // 3    # Create a bit mask of 01 * k-mer length. (e.g. 0101010101, if k-mer length = 5) repeats of '01'
                                                        # Apparently you can do this by getting a base-2 integer, subtracting 1 and int-dividing by 3...tricky
    return (m & lsb_bit_mask).bit_count()               # Apply mask by ANDing on m. Now all mismatches are 01 and can be simply bit-counted to get # mismatches


def find_repeats_2bit(rec_id: str, seq: str, kmer_length: int, allowed_mismatches: int) -> list[RepeatHit]:
    '''Find direct repeats by 2-bit encoding (bitwise comparisons) and nested loop'''
    hits: list[RepeatHit] = []
    codes = encode_kmers(seq, kmer_length)
    
    for i in range(len(codes)):
        for j in range(i + 1, len(codes)):
            mismatches = hamming_distance_2bit(codes[i], codes[j], kmer_length)
            if mismatches <= allowed_mismatches:
                hits.append(RepeatHit(
                    record_id = rec_id,
                    query_start = i + 1,
                    query_end = i + kmer_length,
                    subject_start = j + 1,
                    subject_end = j + kmer_length,
                    query_seq = seq[i:i + kmer_length],
                    subject_seq = seq[j:j + kmer_length],
                    mismatches = mismatches,
                    kmer_length = kmer_length,
                    orientation = 'direct'
                    ))
    return hits

def find_invert_repeats_2bit(rec_id: str, seq: str, kmer_length: int, allowed_mismatches: int) -> list[RepeatHit]:
    '''Find indirect repeats by 2-bit encoding (bitwise comparisons) and nested loop'''
    seq_rc = reverse_complement(seq)
    hits: list[RepeatHit] = []
    codes = encode_kmers(seq, kmer_length)
    codes_rc = encode_kmers(seq_rc, kmer_length)

    for i in range(len(codes)):
        for j in range(len(codes_rc) - i):
            mismatches = hamming_distance_2bit(codes[i], codes_rc[j], kmer_length)
            if mismatches <= allowed_mismatches:
                hits.append(RepeatHit(
                    record_id = rec_id,
                    query_start = i + 1,
                    query_end = i + kmer_length,
                    subject_start = len(seq) - j - kmer_length + 1,
                    subject_end = len(seq) - j,
                    query_seq = seq[i:i + kmer_length],
                    subject_seq = seq[len(seq) - j - kmer_length:len(seq) - j],
                    mismatches = mismatches,
                    kmer_length = kmer_length,
                    orientation = 'inverted'
                    ))
    return hits