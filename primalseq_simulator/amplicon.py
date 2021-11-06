import Bio

class Amplicon( object ):
    def __init__(self, sequence: Bio.SeqRecord, name: str, start: int, stop: int, reads: int ):
        if len( sequence ) != stop - start:
            raise ValueError( f"Sequence length ({len( sequence )} is not equal to specified length ({stop - start})." )
        self.seq = sequence
        self.name = name
        self.start = start
        self.stop = stop

    def __len__( self ):
        return len( self.seq )

