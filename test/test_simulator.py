import pytest

from primalseq_simulator.amplicon import Amplicon
from primalseq_simulator.genome import Genome
from primalseq_simulator.simulator import Simulator

SEQ_REFERENCE = "test/test_genome/wnv_short.fasta"
SEQ_PRIMERS = "test/test_genome/wnv_short.bed"
FIRST_AMPLICON = "AGTAGTGTTTGTGAGGATTAACAACAATTAACACAGTGCGAGCTGTTTCTTAGCACGAAGATCTCGATGTCTAAGAAACCAGGAGGGCCCGGCAAGAGCCGGGCTGTCAATATGCTAAAACGCGGAATGCCCCGCGTGTTGTCCTTGATTGGACTGAAGAGGGCTATGTTGAGCCTGATCGACGGCAAGGGGCCAATACGATTTGTGTTGGCTCTCTTGGCGTTCTTCAGGTTCACAGCAATTGCTCCGACCCGAGCAGTGCTGGATCGATGGAGAGGTGTGAACAAACAAACAGCGATGAAACACCTTCTGAGTTTTAAGAAGGAACT"

@pytest.fixture
def reference_genome():
    return Genome( sequence=SEQ_REFERENCE, primers=SEQ_PRIMERS, reads=200 )

@pytest.fixture
def reference_amplicon():
    return Amplicon( sequence=FIRST_AMPLICON, name="test_amplicon", start=30, stop=359, reads=200 )

def test_simulator_generates_seed():
    sim = Simulator()
    assert sim.seed is not None, "Simulator seed was not generated"