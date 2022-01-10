import os.path
import pytest
from Bio import SeqIO

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
    return Amplicon( sequence=FIRST_AMPLICON, name="test_amplicon", start=30, stop=359, reads=10 )

@pytest.fixture
def simulator():
    return Simulator()

@pytest.fixture
def simulator_amplicon_output( simulator, reference_amplicon ):
    return simulator.simulate_amplicon( amplicon=reference_amplicon, read_length=250 )

@pytest.fixture
def simulator_genome_output( simulator, reference_genome ):
    return simulator.simulate_reads( genome=reference_genome, read_length=250 )

def test_simulator_generates_seed( simulator ):
    assert simulator.seed is not None, "Simulator seed was not generated"

def test_simulator_generates_temorary_directory( simulator ):
    assert simulator.temp_directory.name is not None, "Simulator did not generate temporary directory"

# I'm going to assume that ART is relatively well tested. I.e. the correct read lengths are generated.
def test_simulator_amplicon_outputs_exists( simulator_amplicon_output ):
    results = [os.path.exists( i ) for i in simulator_amplicon_output]
    assert all( results ), "One or more output files does not exists"

def test_sumulator_amplicon_outputs_nonzero_size( simulator_amplicon_output ):
    results = [os.path.getsize(i) > 0 for i in simulator_amplicon_output]
    assert all( results ), f"One of more output files is empty: {simulator_amplicon_output[0]}, {simulator_amplicon_output[1]}"

def test_simulator_amplicon_ouputs_valid_fastq( simulator_amplicon_output ):
    # check to make sure output files are valid fastqs
    try:
        reads1 = [i for i in SeqIO.parse( simulator_amplicon_output[0], "fastq" )]
        reads2 = [i for i in SeqIO.parse( simulator_amplicon_output[1], "fastq" )]
    except ValueError as exc:
        assert False, f"Output files are not valid fastqs: {exc}"

    # Check if output files contain complimentary reads.
    reads1 = [i.name[:-2] for i in reads1]
    reads2 = [i.name[:-2] for i in reads2]
    assert reads1 == reads2, f"Read names are not congruent: {reads1}\n{reads2}"

def test_simulator_genome_ouputs_exist( simulator_genome_output ):
    results = [os.path.exists(i) for i in simulator_genome_output]
    assert all( results ), f"One or more output files does not exist: {simulator_genome_output}"

def test_sumulator_genome_outputs_nonzero_size( simulator_genome_output ):
    results = [os.path.getsize(i) > 0 for i in simulator_genome_output]
    assert all( results ), f"One of more output files is empty: {simulator_genome_output, simulator_genome_output[1]}"

def test_simulator_genome_ouputs_valid_fastq( simulator_genome_output ):
    try:
        reads1 = [i for i in SeqIO.parse( simulator_genome_output[0], "fastq" )]
        reads2 = [i for i in SeqIO.parse( simulator_genome_output[1], "fastq" )]
    except ValueError as exc:
        assert False, f"Output files are not valid fastqs: {exc}"

    # Check if output files contain complimentary reads.
    reads1 = [i.name[:-2] for i in reads1]
    reads2 = [i.name[:-2] for i in reads2]
    assert reads1 == reads2, f"Read names are not congruent: {reads1}\n{reads2}"