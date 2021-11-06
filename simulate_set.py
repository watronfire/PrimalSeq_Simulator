import argparse
import numpy as np
from Bio import Phylo, SeqIO, SeqRecord
from io import StringIO
from simulate_reads import add_substitution

def generate_substitutions( reference, number ):
    """ Generates a given number of substitutions.
    Parameters
    ----------
    reference : Bio.Seq.SeqRecord
        The reference sequence from which to simulate variants on.
    number : int
        The number of substitutions to generate/

    Returns
    -------
    list
        list of substitutions in the form of a tuple: (Base, Substitution, Ancestral)
    """
    mutation_list = list()
    bases_modified = []
    for _ in range( number ):
        bases = ["A", "T", "C", "G"]
        base_to_modify = np.random.randint( 0, len( reference ) )
        while base_to_modify in bases_modified:
            base_to_modify = np.random.randint( 0, len( reference ) )
        ancestral = reference.seq[base_to_modify]
        bases.remove( ancestral )
        substitution = np.random.choice( bases )
        mutation_list.append( (base_to_modify, substitution, ancestral ) )

    return mutation_list


def generate_alignment( reference, tree ):
    """ Will generate leaf sequences from tree with assigned mutations at each node.
    Parameters
    ----------
    reference : Bio.Seq.SeqRecord
        The reference sequence that serves as root for the tree.
    tree : Bio.Phylo.Tree
        The phylogenetic tree.

    Returns
    -------
    list
        List of SeqRecord objects corresponding the leaves of the tree.
    """
    alignment = []
    for term in tree.get_terminals():
        sequence = reference.seq
        for node in tree.get_path( term ):
            for mutation in node.mutations:
                sequence = add_substitution( sequence, *mutation )
        temp = SeqRecord.SeqRecord( sequence, name=term.name )
        alignment.append( SeqRecord.SeqRecord( sequence, id=term.name, name=term.name ) )
    return alignment


def generate_tree( reference, num_samples, p=0.5, output=None ):
    """ Generates a random strictly bifurcating tree with the specified number of leaves.
    Parameters
    ----------
    reference : Bio.Seq.SeqRecord
        The reference sequence from which to simulate variants on.
    num_samples: int
        Number of leaves to place on the tree
    p: float
        Proability of a successive mutation.
    output: Path
        Location to save tree in newick format. If None, then tree isn't written to file.

    Returns
    -------
    Bio.Phylo.Newick.Tree
        Tree containing num_samples leaves and mutations at each node.
    """
    # Generate strictly bifurcating tree represented by a list of lists.
    t = list( range( num_samples ) )
    for i in range( num_samples - 2 ):
        # TODO enable optional polytomies by pulling size from a geometric distribution.
        selection = list( np.random.choice( len( t ), size=2, replace=False ) )
        selection = sorted( selection, reverse=True )
        node = [t.pop( selection[0] ), t.pop( selection[1] )]
        t.append( node )

    # Convert list-based tree into Bio.Phylo
    tree = str( t )
    tree = tree.replace( "[", "(" )
    tree = tree.replace( "]", ")" )
    tree = Phylo.read( StringIO( tree ), "newick" )
    tree.ladderize()
    Phylo.draw_ascii( tree )
    if output is not None:
        Phylo.write( tree, output, "newick" )

    # Add mutations to each node
    for node in tree.find_clades( order="postorder" ):
        # TODO: there should be an option for this to be user specified.
        number_mutations = np.random.geometric( p=p )
        node.mutations = generate_substitutions( reference, number_mutations )

    alignment = generate_alignment( reference, tree )


def main( arguments ):
    if arguments.seed is not None:
        np.random.seed( arguments.seed )

    ref = SeqIO.read( arguments.reference, "fasta" )
    generate_tree( ref, 10 )


if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Generates a set of PrimalSeq data with varying levels of relatedness." )

    # Initialize positional arguments
    parser.add_argument( "reference", help="reference sequence to simulate reads from" )
    parser.add_argument( "primers", help="bed file containing primer sequences", type=str )

    # Initialize optional arguments
    parser.add_argument( "-o", "--output", help="prefix for output files", required=True )
    parser.add_argument( "-a", "--art", help="location of art executable", required=True )
    parser.add_argument( "-s", "--sequencer", help="sequencer to use error profile for." )
    parser.add_argument( "-v", "--variants", action="store_true", help="whether to insert interhost variants in the simulated reads" )
    parser.add_argument( "-vc", "--variants-count", type=int, help="number of variants to simulated" )
    parser.add_argument( "-vs", "--variants-std", type=int, help="specify the expectation of number of variants to simulated from a poisson distribution" )
    parser.add_argument( "-rs", "--seed", type=int, help="random seed to reproduce results" )
    parser.add_argument( "-fs", "--frequency-scale", type=float, help="scale parameter for the exponential distribution used for variant generation" )
    parser.add_argument( "-vt", "--variant-type", type=str, choices=["all", "major", "minor"], default="all", help="type of variant to generate; restricts variant frequency to specific range." )
    parser.add_argument( "-l", "--limit", type=int, default=1000, help="limit of the number of tries the variant sampler gets." )

    # Determine relatedness of samples
    #sample_linkage = parser.add_mutually_exclusive_group( required=True )

    # Determine coverage
    coverage_type = parser.add_mutually_exclusive_group( required=True )
    coverage_type.add_argument( "-r", "--reads", type=int, help="number of reads to generate" )
    coverage_type.add_argument( "-c", "--coverage", type=int, help="mean number of reads to generate per amplicon" )

    args = parser.parse_args()

    main( args )