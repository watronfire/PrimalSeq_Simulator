import numpy as np
from Bio import SeqIO
from primalseq_simulator.primers import Primers

class Genome( object ):
    def __init__( self, seq: str = None , primers: list = None, reads: int = None, coverage: int = None ):

        # Make sure everything is in place before initiating.
        if seq is None:
            raise ValueError( "A sequence must be specified when constructing a genome" )
        if primers is None:
            raise ValueError( "Primers must be given to simulate PrimalSeq reads." )
        if ( reads is not None ) and ( coverage is not None ):
            raise ValueError( "Both reads and coverage cannot be specified" )

        self.seq = SeqIO.read( seq, "fasta" )
        self.primers = Primers( primers )

        if reads is None:
            self.coverage = self.get_coverage( primers, reads=coverage, coverage=True )
        else:
            self.coverage = self.get_coverage( primers, reads=reads, coverage=False )

    @staticmethod
    def get_coverage( primer_scheme : list, reads : int, coverage : bool ) -> np.ndarray:
        """
        Calculates per amplicon coverage given a primer scheme and either a per-amplicon requirement or a per-sample requirement
        Parameters
        ----------
        primer_scheme : list
            A list with length equal to number of amplicons in sample. Contents don't matter at the moment.
        reads : int
            Number of reads to generate in total.
        coverage : bool
            Specifies whether reads refer to per-amplicon coverage (True) or total reads (False).

        Returns
        -------
        np.ndarray
            a list of type:int which correspond to the number of reads to be generated for each amplicon.
        """
        if coverage is not None:
            return np.int_( np.random.normal( reads, reads * 0.05, len( primer_scheme ) ) )
        else:
            distribution = np.ones( len( primer_scheme ) ) / len( primer_scheme )
            return np.random.multinomial( reads, distribution )

