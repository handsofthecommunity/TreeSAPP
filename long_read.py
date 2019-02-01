__author__ = "Kevin Chan"

import logging
import sys
import os
from external_command_interface import launch_write_command
from utilities import reformat_string
from fasta import generate_fasta_or_fastq


def prepare_fasta_target(args, marker_build_dict, index_prefix):
    """
    From a list of target reference packages, concatenate the FASTA files for input to minimap2
    :param marker_build_dict: Dictionary mapping denominator/code names to MarkerBuild objects
    :param index_prefix: Path to write the fasta index
    :return: Path to the fasta index
    """
    fasta_file = os.path.join(index_prefix, "ref.fa")
    reference_package_dir = os.path.join(args.treesapp, "data", "alignment_data")
    with open(fasta_file, "w") as outfile:
        for marker_build_obj in marker_build_dict.values():
            reference_package_path = os.path.join(reference_package_dir, marker_build_obj.cog + ".fa")
            with open(reference_package_path, "r") as infile:
                for record in generate_fasta_or_fastq(infile):
                    header, seq, _ = record
                    seq = seq.replace("-", "")
                    outfile.write("{0}\n{1}\n".format(header, seq))
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
    # TODO: clean up this hacky way
    cmd_list = [minimap_executable, "-t", str(threads), "-x", "map-ont", index_file, query_reads, ">", outpath]
    stdout, minimap_ret_code = launch_write_command(cmd_list)
    if minimap_ret_code != 0:
        logging.error("Long read alignment using {0} did not complete successfully! Command used:\n{1}\n".format(
            minimap_executable, cmd_list))
        sys.exit(5)
    return outpath


def extract_minimap_alignments(args, minimap_matches, fasta_reference_map):
    """
    Extracts regions of query sequence which mapped to a reference, and writes FASTA files for homologous regions.
    See extract_hmm_matches for a detailed description.
    :param args: command line arguments as extracted from get_options and check_parser_arguments
    :param minimap_matches: a dict mapping read names to alignment information
    :param fasta_reference_map: a dict mapping query headers (read names) to query read
    :return:
    """
    logging.info("Extracting the quality-controlled DNA sequences... ")
    minimap_input_fastas = list()
    marker_gene_dict = dict()
    numeric_read_index = dict()
    trimmed_query_bins = dict()
    bins = dict()

    for marker in minimap_matches:
        if marker not in numeric_read_index.keys():
            numeric_read_index[marker] = dict()
        numeric_decrementor = -1
        if marker not in marker_gene_dict:
            marker_gene_dict[marker] = dict()

        for paf_obj in sorted(minimap_matches[marker], key=lambda x: x.end - x.start):
            if paf_obj.qname not in fasta_reference_map:
                logging.error("Read name {0} was aligned using minimap2 but not found in the input fasta!".format(paf_obj.qname))
                sys.exit(5)
            read_coordinates = str(paf_obj.qstart) + "_" + str(paf_obj.qend)
            numeric_read_index[marker][numeric_decrementor] = paf_obj.readname + "_" + read_coordinates
            query_sequence = fasta_reference_map[reformat_string(">" + paf_obj.qname)]
            binned = False
            for bin_num in sorted(bins):
                bin_rep = bins[bin_num][0]
                overlap = min(paf_obj.qend, bin_rep.qend) - max(paf_obj.qstart, bin_rep.qstart)
                if (100 * overlap) / (bin_rep.qend - bin_rep.qstart) > 80:
                    bins[bin_num].append(paf_obj)
                    trimmed_query_bins[bin_num] += ">" + str(numeric_decrementor) + "_" + "\n" + \
                                                   query_sequence[paf_obj.qstart:paf_obj.qend+1] + "\n"
                binned = True
                break
            if not binned:
                bin_num = len(bins)
                bins[bin_num] = [paf_obj]
                trimmed_query_bins[bin_num] += ">" + str(numeric_decrementor) + "_" + "\n" + \
                                               query_sequence[paf_obj.qstart:paf_obj.qend+1] + "\n"
            
            # bulk FASTA header format
            # >qname|marker_gene|start_end
            bulk_header = ">" + paf_obj.qname + "|" + \
                          paf_obj.tname + "|" + \
                          paf_obj.qstart + "_" + paf_obj.qend

            marker_gene_dict[marker][bulk_header] = query_sequence[paf_obj.qstart:paf_obj.qend+1]
            numeric_decrementor -= 1

        for group in trimmed_query_bins:
            if trimmed_query_bins[group]:
                marker_query_fa = marker + "_mm2_purified_group" + str(group) + ".fa"
                marker_query_fa = os.path.join(args.output_dir_var, marker_query_fa)
                try:
                    homolog_seq_fasta = open(marker_query_fa, "w")
                except IOError:
                    logging.error("Unable to open {} for writing.\n".format(marker_query_fa))
                    sys.exit(5)
                minimap_input_fastas.append(marker_query_fa)
                homolog_seq_fasta.write(trimmed_query_bins[group])
                homolog_seq_fasta.close()
        trimmed_query_bins.clear()
        bins.clear()
    logging.info("done.\n")

    for marker in marker_gene_dict:
        trimmed_hits_fasta = marker + "_mm2_purified.fa"
        trimmed_hits_fasta = os.path.join(args.output_dir_var, trimmed_hits_fasta)
        logging.debug("\tWriting {0} sequences to {1}\n".format(marker, trimmed_hits_fasta))
    return minimap_input_fastas, numeric_read_index
