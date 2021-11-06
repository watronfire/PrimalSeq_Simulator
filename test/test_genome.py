import pytest
from Bio import SeqIO

from primalseq_simulator.genome import Genome

SEQ_REFERENCE = "test/test_genome/wnv_short.fasta"
SEQ_PRIMERS = "test/test_genome/wnv_short.bed"
FIRST_AMPLICON = "AGTAGTGTTTGTGAGGATTAACAACAATTAACACAGTGCGAGCTGTTTCTTAGCACGAAGATCTCGATGTCTAAGAAACCAGGAGGGCCCGGCAAGAGCCGGGCTGTCAATATGCTAAAACGCGGAATGCCCCGCGTGTTGTCCTTGATTGGACTGAAGAGGGCTATGTTGAGCCTGATCGACGGCAAGGGGCCAATACGATTTGTGTTGGCTCTCTTGGCGTTCTTCAGGTTCACAGCAATTGCTCCGACCCGAGCAGTGCTGGATCGATGGAGAGGTGTGAACAAACAAACAGCGATGAAACACCTTCTGAGTTTTAAGAAGGAACT"

@pytest.fixture
def reference_genome():
    return Genome( SEQ_REFERENCE, SEQ_PRIMERS, reads=200 )

def test_sequence_type( reference_genome ):
    want = str
    got = type( reference_genome.sequence )
    assert got == want, f"Genome sequences is of type {got} wanted {want}."

def test_error_when_reads_and_coverage_specified():
    with pytest.raises( ValueError ):
        genome = Genome( SEQ_REFERENCE, SEQ_PRIMERS, reads=100, coverage=100 )

def test_correct_number_of_amplicons_generated( reference_genome ):
    want = 2
    got = len( reference_genome.amplicons )
    assert got == want, f"Genome contains {got} amplicons, wanted {want}."

def test_extracts_subsequence_correctly( reference_genome ):
    want = FIRST_AMPLICON
    got = reference_genome.extract_subsequence( 30, 359 )
    assert got == want, f"Genome extractions incorrect.\nGot:\t{got}\nWant:\t{want}"

def test_amplicon_sequence_correctly_generated( reference_genome ):
    want = FIRST_AMPLICON
    got = reference_genome.amplicons[0].seq
    assert got == want, f"Amplicon sequence incorrect.\nGot:\t{got}\nWant:\t{want}"

def test_amplicons_have_reads( reference_genome ):
    assert all( [i.reads > 0 for i in reference_genome.amplicons] ), f"Amplicons not correctly assigned reads"