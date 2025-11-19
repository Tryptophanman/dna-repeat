from dna_repeat.ai import encode_kmers, find_repeats_2bit
from dna_repeat.core import RepeatHit

seq = 'TGAAAGCCAGGCAAGTTTTCTGCTTCTTTTGCTTCTTAGTCAGGAGATAGATAGATTACGTTTTTAGAGTGCCAGGCAAGTCTTCTGCTT'
#REPEAT GCCAGGCAAGTTTTCTGCTT
k = 30

def test_encode_kmers():
    result = encode_kmers(seq, k)
    # for item in result:
    #     print(f'{item:0{2*k}b}')    # f'{number:fill_character{min. width of result}output_format_type}
    numbers = result[0]
    binary = f'{numbers:0{2*k}b}'
    print(numbers, '->', binary)
    assert (numbers == 1009460048132210175) and (binary == '111000000010010100101001000010111111110111100111110111111111')

def test_find_repeats_2bit():
    result2: list[RepeatHit] = find_repeats_2bit('stevo', seq, 20, 1)
    assert result2[0].query_seq == 'GCCAGGCAAGTTTTCTGCTT'
