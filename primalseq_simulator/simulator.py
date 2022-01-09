import os
import random
import tempfile
from subprocess import call
from time import time

from primalseq_simulator.amplicon import Amplicon
from primalseq_simulator.genome import Genome

EXECUTABLE = "res/art_illumina"

class Simulator( object ):
    def __init__(self, seed=None):
        self.temp_directory = tempfile.TemporaryDirectory()

        if seed:
            self.seed = seed
        else:
            self.seed = time()

    def simulate_reads( self, genome: Genome, read_length: int = 250 ):
        """
        Initiates the simulation process. Iterates through the amplicons of the genome and simulates reads. Then
        collates and returns the reads.
        Returns
        -------
        str
            Path to temporary file containing first reads
        str
            Path to temporary file containing second reads
        """
        # Interate through amplicons in Genome.
        temporary_files = list()
        for amplicon in genome.amplicons:
            temp_reads = self._simulate_amplicon( amplicon, read_length )
            temporary_files.append( temp_reads )
        
        # Combines amplicon reads into a single pair of files.
        #combined = self._combine_reads( temporary_files )
        #return combined

    def _simulate_amplicon( self, amplicon: Amplicon, read_length: int ):
        """
        ART wrapper. Takes in an Amplicon object and generates the required number of reads.
        Parameters
        ----------
        amplicon : Amplicon

        Returns
        -------
        str
            Path to temporary files containing simulated reads
        """
        # Extract amplicon sequence
        reference = self._extract_amplicon_to_file( amplicon )

        # Generate output prefix so ART knows how to name output files.
        output = os.path.join( self.temp_directory.name, amplicon.name + "_" )

        # Generate and call the ART command
        command = f"{EXECUTABLE} -ss MSv3 -amp -p -c {amplicon.reads} -l {read_length} -i {reference} -na -rs {self.seed} -o {output}"
        with open( os.devnull, 'w' ) as FNULL:
            call( command, stdout=FNULL, shell=True )

        # Remove the reference because it won't be used after this and we specified delta=False. Potentially unneeded.
        os.remove( reference )

        return f"{output}1.fq", f"{output}2.fq"

    def _extract_amplicon_to_file( self, amplicon ):
        fp = tempfile.NamedTemporaryFile( dir=self.temp_directory.name, suffix=".fasta", mode="w+", delete=False )
        fp.write( "> {}\n".format( amplicon.name ) )
        fp.write( amplicon.seq + "\n" )
        fp.close()
        return fp.name

    def _combine_reads( self, temporary_files ):
        pass