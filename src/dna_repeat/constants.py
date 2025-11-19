# src/dna_repeat/constants.py

COMPLEMENT: dict[int, int] = str.maketrans('ACGT', 'TGCA')       
'''Translation table for nt. base complement'''

BASE_TO_INT: dict[str, int] = {'A':0,'C':1,'G':2,'T':3}
'''Dictionary for substituting an integer for a nt. base'''