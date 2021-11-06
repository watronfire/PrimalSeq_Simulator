import pytest
from primalseq_simulator.primers import Primers

SEQ_PRIMERS = "test/test_genome/wnv_short.bed"

#@pytest.fixture
#def reference_primers():
#    return Primers( SEQ_PRIMERS )

def test_correct_number_primers_parsed():
    primer = Primers( SEQ_PRIMERS )
    wanted = 2
    got = len( primer )
    assert got == wanted, f"Primer scheme has a length of {got}, wanted length of {wanted}."