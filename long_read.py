__author__ = 'Kevin Chan'

import logging
import sys
import os
from external_command_interface import launch_write_command


def prepare_fasta_target(marker_build_dict, index_prefix):
    """
    From a list of target reference packages, concatenate the FASTA files for input to minimap2
    :param marker_build_dict: Dictionary mapping denominator/code names to MarkerBuild objects
    :param index_prefix: Path to write the fasta index
    :return: Path to the fasta index
    """
    fasta_file = os.path.join(index_prefix, "ref.fa")
    reference_package_dir = os.path.join(".", "data", "alignment_data")
    with open(fasta_file, "w") as outfile:
        for marker_build_obj in marker_build_dict.values():
            reference_package_path = os.path.join(reference_package_dir, marker_build_obj.cog)
            with open(reference_package_path, "r") as infile:
                for line in infile:
                    outfile.write(line)
    return fasta_file


def run_minimap(minimap_executable, index_file, query_reads, output_prefix, threads=1):
    """
    Wrapper function for running minimap2 to map long reads to reference index files
    :param minimap_executable: Path to the minimap2 executable, as found by utilities.find_executables()
    :param index_file: Path to a reference minimizer index file for a reference package
    :param query_reads: path to a FASTA or FASTQ file containing the long reads
    :param output_prefix: Prefix string for the SAM file generated
    :param threads: Number of threads to give minimap2 (default is 1)
    :return: Name of the SAM file written
    """
    outpath = os.path.join(output_prefix, "minimap2_aligned.paf")
    # TODO: refactor launch_write_command to keep stdout only
    cmd_list = [minimap_executable, "-t", str(threads), "-x", "map-ont", index_file, query_reads, ">", outpath]
    stdout, minimap_ret_code = launch_write_command(cmd_list)
    if minimap_ret_code != 0:
        logging.error("Long read alignment using {0} did not complete successfully! Command used:\n{1}\n".format(
            minimap_executable, cmd_list))
        sys.exit(5)
    # outfile = open(outpath, "w")
    # outfile.write(stdout)
    # outfile.close()
    return outpath


def extract_minimap_alignments(minimap_matches):
    # TODO: extract reference or query? query.
    return
