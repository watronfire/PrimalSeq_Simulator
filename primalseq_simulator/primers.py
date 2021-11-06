
class Primers( object ):
    """
    Attributes
    ----------
    scheme : list
        a list containing entries for each primer pair with the structure [Primer_name,[Left_bound, Right_bound]].
    """
    def __init__( self, loc ):
        """
        Loads a bed file representing a primer scheme.
        Parameters
        ----------
        loc : str
            path to bed file
        """
        self.scheme = []

        with open( loc, "r" ) as primer_file:
            for line in primer_file:
                line_split = line.split( "\t" )
                left = line_split[1]
                right = line_split[2]
                name = line_split[3]

                if "_LEFT" in name:
                    _ = [int( right )]
                elif "_RIGHT" in name:
                    _.append( int( left ) )
                    self.scheme.append( [name.replace( "_RIGHT", "" ), _] )

    def __len__( self ):
        return len( self.scheme )