import pytest
from Bio import SeqIO

from primalseq_simulator.genome import Genome
from primalseq_simulator.utils import parse_primer_scheme

SEQ_REFERENCE = "test/test_genome/wnv_short.fasta"
SEQ_PRIMERS = "test/test_genome/wnv_short.bed"


@pytest.fixture
def reference_genome():
    return Genome( SEQ_REFERENCE, SEQ_PRIMERS, reads=200 )

def test_empty_init_error():
    with pytest.raises( ValueError ):
        genome = Genome()

def test_sequence_type( reference_genome ):
    assert isinstance( reference_genome.seq, SeqIO.SeqRecord ), f"Genome sequences is of type {type( reference_genome.seq )} wanted SeqIO.SeqRecord"

def test_error_when_reads_and_coverage_specified():
    with pytest.raises( ValueError ):
        genome = Genome( SEQ_REFERENCE, SEQ_PRIMERS, reads=100, coverage=100 )