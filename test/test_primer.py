import pytest
from primalseq_simulator.primer import Primer

SEQ_PRIMERS = "test/test_genome/wnv_short.bed"

#@pytest.fixture
#def reference_primers():
#    return Primers( SEQ_PRIMERS )

def test_correct_number_primers_parsed():
    primer = Primer.parse_primer_bed( SEQ_PRIMERS )
    wanted = 2
    got = len( primer )
    assert got == wanted, f"Primer scheme has a length of {got}, wanted length of {wanted}."