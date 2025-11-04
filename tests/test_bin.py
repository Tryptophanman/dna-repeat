from dna_repeat.core import encode_kmers, find_repeats_2bit

seq = 'TGAAAGCCAGGCAAGTTTTCTGCTTCTTTTGCTTCTTAGTCAGGAGATAGATAGATTACGTTTTTAGAGTGCCAGGCAAGTCTTCTGCTT'
#REPEAT GCCAGGCAAGTTTTCTGCTT
k = 30
result = encode_kmers(seq,k)
# for item in result:
#     print(f'{item:0{2*k}b}')    # f'{number:fill_character{min. width of result}output_format_type}
assert result is not None

result2 = find_repeats_2bit('stevo', seq, 20, 1)
print(result2)
assert result2 is not None
