# sequence_tools.py
# small helper utilities for sequence manipulation and codon/amino-acid counting

from typing import Dict, List, Tuple
import math

CODON_TABLE: Dict[str, str] = {
    # Standard genetic code
    'ATA':'I','ATC':'I','ATT':'I','ATG':'M',
    'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'AAC':'N','AAT':'N','AAA':'K','AAG':'K',
    'AGC':'S','AGT':'S','AGA':'R','AGG':'R',
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L',
    'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
    'CGA':'R','CGC':'R','CGG':'R','CGT':'R',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
    'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'GAC':'D','GAT':'D','GAA':'E','GAG':'E',
    'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
    'TTC':'F','TTT':'F','TTA':'L','TTG':'L',
    'TAC':'Y','TAT':'Y','TAA':'*','TAG':'*',
    'TGC':'C','TGT':'C','TGA':'*','TGG':'W',
}

BASES = set("ATGC")

def clean_seq(s: str) -> str:
    return "".join([c for c in s.upper() if c in BASES])

def split_codons(seq: str) -> List[str]:
    # trim to multiple of 3
    L = (len(seq) // 3) * 3
    seq = seq[:L]
    return [seq[i:i+3] for i in range(0, L, 3)]

def translate(seq: str) -> str:
    codons = split_codons(clean_seq(seq))
    aa = [CODON_TABLE.get(c,'X') for c in codons]
    return ''.join(aa)

def codon_counts(seq: str) -> Dict[str,int]:
    counts = {}
    for c in split_codons(clean_seq(seq)):
        counts[c] = counts.get(c,0) + 1
    return counts

def aa_counts(seq: str) -> Dict[str,int]:
    tr = translate(seq)
    counts = {}
    for a in tr:
        counts[a] = counts.get(a,0) + 1
    return counts

def gc_content(seq: str) -> float:
    s = clean_seq(seq)
    if not s:
        return 0.0
    return 100.0 * (s.count('G') + s.count('C')) / len(s)

def codon_freq_vector(seq: str, codon_order: List[str]) -> List[float]:
    counts = codon_counts(seq)
    total = sum(counts.values()) or 1
    return [counts.get(c,0)/total for c in codon_order]
