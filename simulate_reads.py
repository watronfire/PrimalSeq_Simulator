import argparse
import numpy as np
import tempfile
from subprocess import call, STDOUT
import os
from Bio import SeqIO
import pandas as pd
import sys

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

def get_coverage( primer_scheme, reads, coverage ):
    """ Calculates per amplicon coverage given a primer scheme and either a per-amplicon requirement or a per-sample requirement
    Parameters
    ----------
    primer_scheme : Iterable
        A list with length equal to number of amplicons in sample. Contents don't matter at the moment.
    reads : int
        Number of reads to generate in total.
    coverage : bool
        Specifies whether reads refer to per-amplicon coverage (True) or total reads (False).

    Returns
    -------
    Iterable
        a list of type:int which correspond to the number of reads to be generated for each amplicon.
    """
    if coverage is not None:
        return np.int_( np.random.normal( reads, reads * 0.05, len( primer_scheme ) ) )
    else:
        distribution = np.ones( len( primer_scheme ) ) / len( primer_scheme )
        return np.random.multinomial( reads, distribution )

def add_substitution( sequence, base, alternative, ancestral ):
    """ Function to add a substitution to a str sequence
    Parameters
    ----------
    sequence: str
        The sequence to be modified.
    base: int
        The index position of the base to be substituted.
    alternative: char
        The character to be substituted.
    ancestral: char
        The expected character in sequence at base. Confirms variant table and method is accurate.

    Returns
    -------
    str
        The input sequence with the specified substituted performed.
    """
    if sequence[base] != ancestral:
        sys.exit( "Base in sequence isn't equal to ancestral. Was the variant table generated from this reference?" )

    return sequence[:base] + alternative + sequence[base+1:]

# TODO: need a variant filtering step to remove variant combinations which aren't possible, i.e. where sum( ALT_FREQ ) \
#  is greater than 1.
def extract_amplicon( temp_location, reference, primer_pair, variants=None, reads_assigned=None ):
    """ Extracts the region amplified by a primer pair from the reference sequence and saves to a temporary fasta file.
    Parameters
    ----------
    temp_location : tempfile.TemporaryDirectory
        Location of temporary directory.
    reference : Bio.Seq.SeqRecord
        Reference sequence from which to extract amplicon.
    primer_pair : Iterable
        An entry from the primer pair list with the structure [Primer_name,[Left_bound, Right_bound]].
    variants : pandas.DataFrame, optional
        Table in Variant Call Format containing simulated variants.
    reads_assigned : int, optional
        Total reads to generate for the amplicon.

    Returns
    -------
    Iterable
        A list containing an entry for variant of the amplicon containing the location of the fasta file corresponding to the extracted amplicon and a int of the number of reads to generate.
    """
    seq = str( reference.seq[primer_pair[1][0]:primer_pair[1][1]] )

    seqs = list()
    reads_required = 0

    # If variants are required then generate the variant sequences and calculate the required reads.
    if variants is not None:
        amplicon_variants = variants.loc[variants["POS"].between( primer_pair[1][0], primer_pair[1][1])]
        if len( amplicon_variants ) > 0:
            for index, row in amplicon_variants.iterrows():

                var_amplicon_position = row["POS"] - primer_pair[1][0]

                seq_temp = add_substitution( seq, var_amplicon_position, row["ALT"], row["REF"] )
                reads_temp = np.round( reads_assigned * row["ALT_FREQ"] )
                reads_required += reads_temp

                seqs.append( [seq_temp, reads_temp] )

        if reads_required > reads_assigned + 5:
            sys.exit( "Reads requested cannot be generated. Known limitation of purely unlinked variants." )


    seqs.append( [seq, reads_assigned - reads_required] )

    # Actually write the sequences to file.
    return_list = list()
    for i, j in enumerate( seqs ):
        fp = tempfile.NamedTemporaryFile( dir=temp_location, suffix=".fasta", mode="w+", delete=False )
        fp.write( "> {}\n".format( primer_pair[0] + "_v{}".format( i )  ) )
        fp.write( j[0] + "\n" )
        fp.close()
        return_list.append( [fp.name, j[1]] )
    return return_list

def art_wrapper( executable, temp_location, reference, reads, read_length, sequencer="illumina", seed=None, verbose=False ):
    """ Basic wrapper for the ART simulator.
    Parameters
    ----------
    executable : str
        Location of the ART executable if not in path.
    temp_location : str
        Location of temporary directory.
    reference : str
        Location of fasta file containg sequence of extracted amplicon.
    reads : int
        Number of reads to generate.
    read_length : int
        Length of read to generate
    sequencer : str
        ART-specific name of the sequencer to simulate. Specifies the error profile to use.
    seed : int
        Random seed. Can be specified to reproduce results.
    verbose : bool
        Whether to print end of run summary or not.

    Returns
    -------
    list
        A list contain the names of the outputted files.

    """
    # Get random seed
    if seed is None:
        seed = np.random.get_state()[1][1]

    # Output name should be derived from the input reference name rather than just the primer_pair.
    reference_name = os.path.splitext( os.path.split( reference )[1] )[0]
    output = os.path.join( temp_location, reference_name + "_" )

    # Construct command
    command = "{} -ss MSv3 -amp -p -c {} -l {} -i {} -na -rs {} -o {}".format(
        executable,
        reads,
        read_length,
        reference,
        seed,
        output
    )

    if not verbose:
        command += " -q"

    with open(os.devnull, 'w') as FNULL:
        call( command, stdout=FNULL, shell=True )

    # Remove the amplicon fasta. Shouldn't have to use after this.
    os.remove( reference )

    return ["{}1.fq".format( output ), "{}2.fq".format( output )]

# TODO: This function should write files in compressed format. Use the gzip library to open and write output files.
def combine_and_output( amplicon_locations, output_prefix ):
    """ Combines amplicon reads into a single pair of files.
    Parameters
    ----------
    amplicon_locations : list
        Location of temporary directory containing simulated reads for each amplicon.
    output_prefix : str
        Prefix to append to output file.
    """
    # Iterate through simulated files and append to output file.
    # write the first reads to file
    with open( output_prefix + "_1.fastq", "w" ) as first_reads, open( output_prefix + "_2.fastq", "w" ) as second_read:
        for sim_pair in amplicon_locations:
            with open( sim_pair[0], "r" ) as first:
                for line in first: first_reads.write( line )
            with open( sim_pair[1], "r" ) as second:
                for line in second: second_read.write( line )

def generate_variants( reference, count=None, count_exp=None, frequency_scale=None, var_type="all", retry_limit=1000 ):
    """ Generates a simulated vcf file to use for use in simulating variants. Will be expanded to with optional manual \
    functionality. Currently samples variant frequency from an inverse exponential distribution.
    Parameters
    ----------
    reference : Bio.Seq.SeqRecord
        The reference sequence from which to simulate variants on.
    count : int, optional
        Optional. The number of variants to simulated. If not specified then the number of variants is drawn from a poisson distribution.
    count_exp : int, optional
        The expected number of variants for the poisson distribution. Not used if count is specified.
    frequency_scale : int, optional
        The scale parameter to use to parameterize the exponential distribution. Higher value will induce higher frequencies.
    var_type : {'all', 'major', 'minor'}
        Indicates what type of variants to generate. Useable options are: "minor" - alternative frequency < 0.5, "major" - alternative frequency > 0.5, "all" - alternative frequency < 1.
    retry_limit : int, optional
        The number of tries the sampler can used to generate a suitable alternative frequency according to var_type.

    Returns
    -------
    pd.Dataframe
        Table in Variant Call Format containing simulated variants.
    """
    # Calculate the number of variants to generate. If not specified draw from a poisson distribution centered around county_exp.
    if count is not None:
        variants = count
    else:
        variants = np.random.poisson( count_exp ) if count_exp is not None else np.random.poisson( 10 )

    # Intialize variant requirements depending on the variant type requested.
    # TODO determine an adequate frequency scale for generating major/all variants.
    if var_type == "minor":
        frequency_requirement = lambda x: x < 0.5
        if frequency_scale is None:
            frequency_scale = 0.09
    elif var_type == "major":
        frequency_requirement = lambda x : 1 > x > 0.5
        if frequency_scale is None:
            frequency_scale = 0.12
    elif var_type == "all":
        frequency_requirement = lambda x : True
        if frequency_scale is None:
            frequency_scale = 0.1
    else:
        sys.exit( "'{}' is not a known variant type. Valid options are ['minor', 'major', 'all']".format( var_type ) )

    var_table = list()

    # For each variant required call a random base to modify, a random substitution, and a random frequency.
    for i in range( variants ):
        bases = ["A", "T", "C", "G"]

        base_to_modify = np.random.randint( 0, len( reference ) )
        ancestral = reference.seq[base_to_modify]
        bases.remove( ancestral )
        minor_variant = np.random.choice( bases )
        alternative_frequency = np.random.exponential( scale=frequency_scale )

        # The frequency samplier isn't perfect, so we draw from it until we get a reasonable frequency. If this doesn't \
        # happen, it means you're either really unlucky or the frequency scale is off.
        retry_count = 0
        while not frequency_requirement( alternative_frequency ):
            alternative_frequency = np.random.exponential( scale=frequency_scale )
            if retry_count > retry_limit:
                sys.exit( "Frequency sample limit reached. Distribution scale is either too high for minor variants, or too low for major variants" )
            retry_count += 1

        var_table.append( ( reference.name, base_to_modify, ancestral, minor_variant, alternative_frequency ) )

    # Return a pseudo-VCF table.
    return_df = pd.DataFrame( var_table, columns=["CHROM", "POS", "REF", "ALT", "ALT_FREQ"] )



    return pd.DataFrame( var_table, columns=["CHROM", "POS", "REF", "ALT", "ALT_FREQ"] )

def simulate_reads( executable, primers, ref, output, coverage=None, variants=False, variant_count=None,
                    variant_exp=None, variants_table=None, seed=None, reads=None, frequency_scale=None, var_type="all",
                    retry_limit=1000 ):
    """
    Parameters
    ----------
    executable : str
        Location of the ART executable if not in path.
    primers : str
        The location of the primer bed file
    ref : Bio.Seq.SeqRecord
        The reference sequence from which to simulate variants on.
    output : str
        Prefix to append to output file.
    coverage : int
        Mean number of reads to generate per amplicon.
    variants : bool
        Whether to generate variants or not.
    variant_count : int
        Number of variants to generate.
    variant_exp : int
        Expected number of variants to generate. Used when variant_count not specified and the number of variants is \
        drawn from a poisson distribution.
    variants_table : pandas.DataFrame
        Specify the variants to generate. Will override variants, variant_count, and variant_exp.
    seed : int
        The random seed to use for the simulation.
    reads : int
        The total number of reads to generate. If specified, will override coverage value.
    frequency_scale : int, optional
        The scale parameter to use to parameterize the exponential distribution. Higher value will induce higher frequencies.
    var_type : {'all', 'major', 'minor'}
        Indicates what type of variants to generate. Useable options are: "minor" - alternative frequency < 0.5, "major" - alternative frequency > 0.5, "all" - alternative frequency < 1.
    retry_limit : int, optional
        The number of tries the sampler can used to generate a suitable alternative frequency according to var_type.

    Returns
    -------
    None
        paired fastq files will be generated with output prefix.
    """
    # Initialize random seed if request.
    if seed is not None:
        np.random.seed( seed )
    else:
        # My plan is to print this at the end so you can rerun, but I don't know if that'll do what I want.
        # TODO: test seed saving ability.
        seed_used = np.random.get_state()

    primer_pairs = parse_primer_scheme( primers )

    if reads is None:
        amplicon_cov = get_coverage( primer_pairs, reads=coverage, coverage=True )
    else:
        amplicon_cov =  get_coverage( primer_pairs, reads=reads, coverage=False )

    if variants:
        if variants_table is None:
            variants_table = generate_variants( ref, count=variant_count, count_exp=variant_exp,
                                                frequency_scale=frequency_scale, var_type=var_type,
                                                retry_limit=retry_limit )
        variants_table.to_csv( output + ".vcf", sep="\t", index=False )

    return_list = list()
    with tempfile.TemporaryDirectory() as tmpdirname:
        print( "created temporary directory", tmpdirname )
        for i, pair in enumerate( primer_pairs ):
            af = extract_amplicon( temp_location=tmpdirname,
                                   reference=ref,
                                   primer_pair=pair,
                                   variants=variants_table,
                                   reads_assigned=amplicon_cov[i] )
            for j in af:
                return_list.append( art_wrapper( executable=executable,
                                                 temp_location=tmpdirname,
                                                 reference=j[0],
                                                 reads=j[1],
                                                 read_length=250,
                                                 seed=seed) )
        combine_and_output( return_list, output )

# TODO: Move some of the functionality of this into simulate reads. Right now it is difficult to use the script as an
#  import.
def main( arguments ):
    """
    Parameters
    ----------
    arguments : Namespace
        Output from the argparse.parse_args() function.

    Returns
    -------
    None
        paired fastq files will be generated with output prefix.
    """
    ref = SeqIO.read( arguments.reference, "fasta" )

    simulate_reads( executable=arguments.art, primers=arguments.primers, ref=ref, output=arguments.output,
                    coverage=arguments.coverage, reads=arguments.reads, variants=arguments.variants,
                    variant_count=arguments.variants_count, variant_exp=arguments.variants_std,
                    seed=arguments.seed, frequency_scale=arguments.frequency_scale, var_type=arguments.variant_type,
                    retry_limit=arguments.limit )

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="" )

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

    # Determine coverage
    coverage_type = parser.add_mutually_exclusive_group( required=True )
    coverage_type.add_argument( "-r", "--reads", type=int, help="number of reads to generate" )
    coverage_type.add_argument( "-c", "--coverage", type=int, help="mean number of reads to generate per amplicon" )

    args = parser.parse_args()

    main( args )