import Bio

class Amplicon( object ):
    def __init__(self, sequence: str, name: str, start: int, stop: int, reads: int ):
        if len( sequence ) != stop - start:
            raise ValueError( f"Sequence length ({len( sequence )} is not equal to specified length ({stop - start})." )
        self.seq = sequence
        self.name = name
        self.start = start
        self.stop = stop
        self.reads = reads

    def __len__( self ):
        return len( self.seq )

