from dna_repeat.core import find_repeats, hamming_distance

def test_find():
    seq = 'TTTATCTATTCTTGAAAAAAACGACTTTTTCTATTCTTGAACAAA'   # has 20-bp repeat with 2 mismatches
    k = 20
    m = 2
    run = find_repeats(seq, k, m)
    print(run)
    assert any(mis <= m for _,_,mis in run)
    assert len(run) > 0
    assert hamming_distance(run[0][0], run[0][1]) == 2