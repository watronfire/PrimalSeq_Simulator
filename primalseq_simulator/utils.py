
# This should be removed but its the only think in the file at the moment.
def parse_primer_scheme( primer_loc ):
    """ Loads and parses primer bed file into a list of amplicon boundaries.
    Parameters
    ----------
    primer_loc : str
        The location of the primer bed file

    Returns
    -------
    list
        a list containing entries for each primer pair with the structure [Primer_name,[Left_bound, Right_bound]].
    """
    primer_list = list()

    with open( primer_loc, "r" ) as primer_file:
        for line in primer_file:
            line_split = line.split( "\t" )
            left = line_split[1]
            right = line_split[2]
            name = line_split[3]

            if "_LEFT" in name:
                _ = [ int( right ) ]
            elif "_RIGHT" in name:
                _.append( int( left )  )
                primer_list.append( [name.replace( "_RIGHT", "" ), _] )

        return primer_list