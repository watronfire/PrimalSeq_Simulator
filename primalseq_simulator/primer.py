
class Primer( object ):
    """
    Attributes
    ----------
    scheme : list
        a list containing entries for each primer pair with the structure [Primer_name,[Left_bound, Right_bound]].
    """
    def __init__( self, name, left, right ):
        """
        Loads a bed file representing a primer scheme.
        Parameters
        ----------
        loc : str
            path to bed file
        """
        self.name = name
        self.left = left
        self.right = right

    @staticmethod
    def parse_primer_bed( loc ):
        scheme = list()

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

                    scheme.append( Primer( name.replace( "_RIGHT", "" ), *_ ) )
        return scheme