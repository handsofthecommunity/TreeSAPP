#!/usr/bin/env python3

import sys
import os
import re
import logging
import json
from collections import namedtuple

from .classy import MarkerBuild, Cluster, ReferencePackage
from .fasta import read_fasta_to_dict
from .utilities import get_hmm_length
from . import HMMER_domainTblParser

__author__ = 'Connor Morgan-Lang'


def parse_ref_build_params(base_dir: str, targets=None):
    """
    Returns a dictionary of MarkerBuild objects storing information pertaining to the build parameters of each marker.

    :param base_dir: Path to the treesapp package directory containing 'data/ref_build_parameters.tsv'
    :param targets: List of refpkg codes that are desired or an empty list suggesting all refpkgs should be used
    """
    ref_build_parameters = base_dir + os.sep + "data" + os.sep + 'ref_build_parameters.tsv'
    if type(targets) is str:
        targets = [targets]

    try:
        param_handler = open(ref_build_parameters, 'r')
    except IOError:
        logging.error("\tUnable to open " + ref_build_parameters + " for reading.\n")
        sys.exit(5)

    header_fields = ["name", "code", "molecule", "sub_model", "marker_info", "cluster_identity", "ref_sequences",
                     "tree_tool", "poly-params", "lowest_reliable_rank", "last_updated", "description"]
    header_re = re.compile("\t".join(header_fields))
    if not header_re.match(param_handler.readline().strip()):
        logging.error("Header of '" + ref_build_parameters + "' is unexpected!")
        sys.exit(5)

    logging.debug("Reading build parameters of reference markers... ")
    skipped_lines = []
    missing_info = []
    marker_build_dict = dict()
    for line in param_handler:
        if header_re.match(line):
            continue
        if line[0] == '#':
            skipped_lines.append(line)
            continue
        marker_build = MarkerBuild()
        marker_build.load_build_params(line, len(header_fields))
        if targets and marker_build.denominator not in targets:
            skipped_lines.append(line)
        else:
            if marker_build.denominator in marker_build_dict:
                logging.debug("Multiple '" + marker_build.denominator + "' codes in " + ref_build_parameters +
                              ". Previous entry in marker_build_dict being overwritten...\n")
            marker_build_dict[marker_build.denominator] = marker_build
            if marker_build.load_pfit_params(line):
                missing_info.append(marker_build)
            marker_build.check_rank()
    param_handler.close()

    logging.debug("done.\n")

    if missing_info:
        logging.debug("Rank distance information missing for:\n\t" +
                      "\n\t".join([mb.cog + '-' + mb.denominator for mb in missing_info]) + "\n")
    if skipped_lines:
        logging.debug("Skipped the following lines:\n\t" +
                      "\n\t".join([line.strip() for line in skipped_lines]) + "\n")

    if len(marker_build_dict) == 0:
        logging.error("No reference package information was parsed.\n" +
                      "Is your target '" + ','.join(targets) + "' in " + ref_build_parameters + "?\n")
        sys.exit(3)
    return marker_build_dict


def load_json_build(json_file: str) -> dict:
    try:
        json_handler = open(json_file, 'r')
    except IOError:
        logging.error("Unable to open JSON file '%s' for reading.\n" % json_file)
        sys.exit(3)

    build_params = json.load(json_handler)
    json_handler.close()

    return build_params


def gather_ref_packages(treesapp_dir: str, marker_build_dict: dict, targets=None) -> dict:
    """
    Returns a dictionary of ReferencePackage instances for each MarkerBuild instance in the marker_build_dict.
    Optionally these markers can be subsetted further by a list of targets

    :param treesapp_dir: Path to the treesapp package directory containing 'data/ref_build_parameters.tsv'
    :param marker_build_dict: A dictionary of MarkerBuild instances indexed by their refpkg codes/denominators
    :param targets: List of refpkg codes that are desired or an empty list suggesting all refpkgs should be used
    """
    refpkg_dict = dict()
    refpkg_data_dir = treesapp_dir + os.sep + "data" + os.sep
    logging.debug("Gathering reference package files... ")

    for denominator in marker_build_dict:
        marker_build = marker_build_dict[denominator]  # type: MarkerBuild
        refpkg = ReferencePackage()
        refpkg.refpkg_code = denominator
        refpkg.prefix = marker_build.cog
        if targets:
            if refpkg.prefix not in targets and refpkg.refpkg_code not in targets:
                continue
        refpkg.num_seqs = marker_build.num_reps
        refpkg.gather_package_files(refpkg_data_dir, marker_build.molecule)
        refpkg.profile_length = get_hmm_length(refpkg.profile)
        refpkg_dict[denominator] = refpkg

    logging.debug("done.\n")

    if len(refpkg_dict) == 0:
        logging.error("No reference package data was found.\n" +
                      "Are there reference packages in '" + refpkg_data_dir + "'?\n")
        sys.exit(3)
    return refpkg_dict


def read_graftm_classifications(assignment_file):
    """
    Function for reading the _read_tax.tsv file generated by graftM.
    Sequences that have either multiple genes and/or subunits encoded or have homologous regions separated by
     a divergent sequence are recorded as "_split_" and only the first split is recorded and analyzed.

    :param assignment_file: Path to the _read_tax.tsv file
    :return: Dictionary indexed by taxonomic lineage whose values are headers of classified sequences
    """
    assignments = dict()
    try:
        assignments_handle = open(assignment_file, 'r')
    except IOError:
        logging.error("Unable to open classification file '" + assignment_file + "' for reading.\n")
        sys.exit(21)
    tax_lines = assignments_handle.readlines()
    assignments_handle.close()

    for line in tax_lines:
        fields = line.strip().split('\t')
        try:
            header, classified = fields
            if re.search("_split_", header):
                split = int(header.split('_')[-1])
                if split > 1:
                    continue
                else:
                    header = re.sub("_split_.*", '', header)
            classified = '; '.join([re.sub('e\d+$', '', taxon) for taxon in classified.split('; ')])
            if header and classified:
                if classified not in assignments:
                    assignments[classified] = list()
                assignments[classified].append(header)
        except ValueError:
            logging.error("Unable to parse line:" + str(line))
            sys.exit(21)

    return assignments


def parse_assignments(classified_lines: list):
    """
    Parses the marker_contig_map.tsv lines loaded to retrieve lineage assignment and marker information

     Now also looks for fragments of identical parent sequences that were individually classified.
     If these are found, the longest fragment is selected for classification.
     Removing redundant fragments is performed to ease downstream classifications as
     the number of query sequences is calculated separately and we don't want to exceed 100% classifications!
      Alternatively, the number of query sequences could be calculated from the classification tables
      but we don't think this is the best route as unclassified seqs would wreak havoc.

    :param classified_lines: A list of classification lines returned by read_marker_classification_table
    :return: A dictionary of lineage information for each assignment, indexed by the marker gene it was classified as
    """
    classified = namedtuple("classified", ["refpkg", "taxon", "length"])
    assignments = dict()
    unique_headers = dict()  # Temporary storage for classified sequences prior to filtering
    dups = set()  # For storing the names of all query sequences that were split and classified separately
    for fields in classified_lines:
        _, header, marker, length, raw_tax, rob_class, _, _, _, _, _ = fields
        if marker and rob_class:
            if marker not in assignments:
                assignments[marker] = dict()
            if rob_class not in assignments[marker]:
                assignments[marker][rob_class] = list()
            if header not in unique_headers:
                unique_headers[header] = None
            # If fragments from the same parent query had identical lengths these would be overwritten anyway
            unique_headers[header] = {int(length): classified(refpkg=marker, taxon=rob_class, length=int(length))}
        else:
            logging.error("Bad line in classification table - no robust taxonomic classification:\n" +
                          '\t'.join(fields) + "\n")
            sys.exit(21)
    for header in unique_headers:
        if len(unique_headers[header]) > 1:
            dups.add(header)
        max_len = max(unique_headers[header])
        best_dat = unique_headers[header][max_len]
        assignments[best_dat.refpkg][best_dat.taxon].append(header)

    if dups:
        logging.debug(str(len(dups)) + " fragments from identical parent sequences were removed post-classification.\n")
    return assignments


def read_marker_classification_table(assignment_file):
    """
    Function for reading the tabular assignments file (currently marker_contig_map.tsv)
    Assumes column 2 is the TreeSAPP assignment and column 3 is the sequence header
    (leaving 1 for marker name and 4 for numerical abundance)

    :param assignment_file: Path to the file containing sequence phylogenetic origin and assignment
    :return: A list of lines that have been split by tabs into lists themselves
    """
    classified_lines = list()
    header = "Sample\tQuery\tMarker\tLength\tTaxonomy\tConfident_Taxonomy\tAbundance\tiNode\tLWR\tEvoDist\tDistances\n"

    try:
        assignments_handle = open(assignment_file, 'r')
    except IOError:
        logging.error("Unable to open classification file '" + assignment_file + "' for reading.\n")
        sys.exit(21)
    # This is the header line
    if assignments_handle.readline() != header:
        logging.error("Header of assignments file is unexpected!\n")
        sys.exit(21)

    # First line in the table containing data
    line = assignments_handle.readline()
    n_fields = len(header.split("\t"))
    while line:
        fields = line.strip().split('\t')
        if len(fields) == n_fields:
            classified_lines.append(fields)
        else:
            logging.error("Unable to parse line:\n" + str(line))
            sys.exit(21)
        line = assignments_handle.readline()
    assignments_handle.close()

    return classified_lines


def best_discrete_matches(matches: list) -> list:
    """
    Function for finding the best alignment in a list of HmmMatch() objects
    The best match is based off of the full sequence score
    :param matches: A list of HmmMatch() objects
    :return: List of the best HmmMatch's
    """
    # Code currently only permits multi-domains of the same gene
    dropped_annotations = list()
    len_sorted_matches = sorted(matches, key=lambda x: x.end - x.start)
    i = 0
    orf = len_sorted_matches[0].orf
    while i+1 < len(len_sorted_matches):
        j = i + 1
        a_match = len_sorted_matches[i]  # type HmmMatch
        while j < len(len_sorted_matches):
            b_match = len_sorted_matches[j]  # type HmmMatch
            if a_match.target_hmm != b_match.target_hmm:
                if HMMER_domainTblParser.detect_orientation(a_match.start, a_match.end,
                                                                     b_match.start, b_match.end) != "satellite":
                    if a_match.full_score > b_match.full_score:
                        dropped_annotations.append(len_sorted_matches.pop(j))
                        j -= 1
                    else:
                        dropped_annotations.append(len_sorted_matches.pop(i))
                        j = len(len_sorted_matches)
                        i -= 1
            j += 1
        i += 1

    if len(len_sorted_matches) == 0:
        logging.error("All alignments were discarded while deciding the best discrete HMM-match.\n")
        sys.exit(3)

    logging.debug("HMM search annotations for " + orf +
                  ":\n\tRetained\t" + 
                  ', '.join([match.target_hmm +
                             " (%d-%d)" % (match.start, match.end) for match in len_sorted_matches]) +
                  "\n\tDropped\t\t" +
                  ', '.join([match.target_hmm +
                             " (%d-%d)" % (match.start, match.end) for match in dropped_annotations]) + "\n")
    return len_sorted_matches


def parse_domain_tables(args, hmm_domtbl_files: list) -> dict:
    """

    :param args:
    :param hmm_domtbl_files: A list of domain table files written by hmmsearch
    :return: Dictionary of HmmMatch objects indexed by their reference package and/or HMM name
    """
    # Check if the HMM filtering thresholds have been set
    thresholds = HMMER_domainTblParser.prep_args_for_parsing(args)

    logging.info("Parsing HMMER domain tables for high-quality matches... ")

    search_stats = HMMER_domainTblParser.HmmSearchStats()
    hmm_matches = dict()
    orf_gene_map = dict()
    optional_matches = list()

    # TODO: Capture multimatches across multiple domain table files
    for domtbl_file in hmm_domtbl_files:
        prefix, reference = re.sub("_domtbl.txt", '', os.path.basename(domtbl_file)).split("_to_")
        domain_table = HMMER_domainTblParser.DomainTableParser(domtbl_file)
        domain_table.read_domtbl_lines()
        distinct_hits = HMMER_domainTblParser.format_split_alignments(domain_table, search_stats)
        purified_hits = HMMER_domainTblParser.filter_poor_hits(thresholds, distinct_hits, search_stats)
        complete_hits = HMMER_domainTblParser.filter_incomplete_hits(thresholds, purified_hits, search_stats)
        HMMER_domainTblParser.renumber_multi_matches(complete_hits)

        for match in complete_hits:
            match.genome = reference
            if match.orf not in orf_gene_map:
                orf_gene_map[match.orf] = dict()
            try:
                orf_gene_map[match.orf][match.target_hmm].append(match)
            except KeyError:
                orf_gene_map[match.orf][match.target_hmm] = [match]
            if match.target_hmm not in hmm_matches.keys():
                hmm_matches[match.target_hmm] = list()
    search_stats.num_dropped()
    for orf in orf_gene_map:
        if len(orf_gene_map[orf]) == 1:
            for target_hmm in orf_gene_map[orf]:
                for match in orf_gene_map[orf][target_hmm]:
                    hmm_matches[target_hmm].append(match)
                    search_stats.seqs_identified += 1
        else:
            search_stats.multi_alignments += 1
            # Remove all the overlapping domains - there can only be one highlander
            for target_hmm in orf_gene_map[orf]:
                optional_matches += orf_gene_map[orf][target_hmm]
            retained = 0
            for discrete_match in best_discrete_matches(optional_matches):
                hmm_matches[discrete_match.target_hmm].append(discrete_match)
                retained += 1
            search_stats.dropped += (len(optional_matches) - retained)
            search_stats.seqs_identified += retained
            optional_matches.clear()

    logging.info("done.\n")

    alignment_stat_string = search_stats.summarize()

    if search_stats.seqs_identified == 0 and search_stats.dropped == 0:
        logging.warning("No alignments found! TreeSAPP is exiting now.\n")
        sys.exit(0)
    if search_stats.seqs_identified == 0 and search_stats.dropped > 0:
        logging.warning("No alignments (" + str(search_stats.seqs_identified) + '/' + str(search_stats.dropped) +
                        ") met the quality cut-offs! TreeSAPP is exiting now.\n")
        alignment_stat_string += "\tPoor quality alignments:\t" + str(search_stats.bad) + "\n"
        alignment_stat_string += "\tShort alignments:\t" + str(search_stats.short) + "\n"
        logging.debug(alignment_stat_string)
        sys.exit(0)

    alignment_stat_string += "\n\tNumber of markers identified:\n"
    for marker in sorted(hmm_matches):
        alignment_stat_string += "\t\t" + marker + "\t" + str(len(hmm_matches[marker])) + "\n"
        # For debugging:
        # for match in hmm_matches[marker]:
        #     match.print_info()

    logging.debug(alignment_stat_string)
    return hmm_matches


def read_colours_file(annotation_file: str, refpkg_name: str) -> (dict, bool):
    """
    Read annotation data from 'annotation_file' and store it in marker_subgroups under the appropriate
    marker and data_type.

    :param annotation_file: Path to an iTOL-compatible annotation file (e.g. colours_styles_file.txt)
    :param refpkg_name: Name of the reference package
    :return: A dictionary of lists where each list is populated by tuples with start and end leaves
    """
    try:
        style_handler = open(annotation_file, 'r')
    except IOError:
        logging.error("Unable to open " + annotation_file + " for reading!\n")
        sys.exit(5)

    clusters = dict()
    field_sep = ''
    internal_nodes = False

    line = style_handler.readline()
    # Skip the header
    while line.strip() != "DATA":
        header_fields = line.strip().split(' ')
        if header_fields[0] == "SEPARATOR":
            if header_fields[1] == "SPACE":
                field_sep = ' '
            elif header_fields[1] == "TAB":
                field_sep = '\t'
            else:
                logging.error("Unknown separator used in " + annotation_file + ": " + header_fields[1] + "\n")
                sys.exit(5)
        line = style_handler.readline()
    # For RGB
    range_line_rgb = re.compile(r"^(\d+)_" + re.escape(refpkg_name) + r"\|(\d+)_" + re.escape(refpkg_name) + re.escape(field_sep) +
                                "range" + re.escape(field_sep) +
                                r".*\)" + re.escape(field_sep) +
                                "(.*)$")
    single_node_rgb = re.compile(r"^(\d+)_" + re.escape(refpkg_name) + re.escape(field_sep) +
                                 "range" + re.escape(field_sep) +
                                 r".*\)" + re.escape(field_sep) +
                                 "(.*)$")
    internal_node_rgb = re.compile(r"^(\d+)" + re.escape(field_sep) +
                                   "range" + re.escape(field_sep) +
                                   r".*\)" + re.escape(field_sep) +
                                   "(.*)$")

    # For hexadecimal
    range_line = re.compile(r"^(\d+)_" + re.escape(refpkg_name) + r"\|(\d+)_" + re.escape(refpkg_name) + re.escape(field_sep) +
                            "range" + re.escape(field_sep) +
                            "#[0-9A-Za-z]{6}" + re.escape(field_sep) +
                            "(.*)$")
    single_node = re.compile(r"^(\d+)_" + re.escape(refpkg_name) + re.escape(field_sep) +
                             "range" + re.escape(field_sep) +
                             "#[0-9A-Za-z]{6}" + re.escape(field_sep) +
                             "(.*)$")
    internal_node = re.compile(r"^(\d+)" + re.escape(field_sep) +
                               "range" + re.escape(field_sep) +
                               "#[0-9A-Za-z]{6}" + re.escape(field_sep) +
                               "(.*)$")

    # Begin parsing the data from 4 columns
    line = style_handler.readline().strip()
    while line:
        if range_line.match(line):
            style_data = range_line.match(line)
            start, end, description = style_data.groups()
        elif range_line_rgb.match(line):
            style_data = range_line_rgb.match(line)
            start, end, description = style_data.groups()
        elif single_node.match(line):
            style_data = single_node.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
        elif single_node_rgb.match(line):
            style_data = single_node_rgb.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
        elif internal_node.match(line):
            style_data = internal_node.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
            internal_nodes = True
        elif internal_node_rgb.match(line):
            style_data = internal_node_rgb.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
            internal_nodes = True
        else:
            logging.error("Unrecognized line formatting in " + annotation_file + ":\n" + line + "\n")
            sys.exit(5)

        description = style_data.groups()[-1]
        if description not in clusters.keys():
            clusters[description] = list()
        clusters[description].append((start + "_" + refpkg_name,
                                      end + "_" + refpkg_name))

        line = style_handler.readline().strip()

    style_handler.close()

    logging.debug("\tParsed " + str(len(clusters)) + " clades from " + annotation_file + "\n")

    return clusters, internal_nodes


def xml_parser(xml_record, term):
    """
    Recursive function for parsing individual xml records

    :param xml_record:
    :param term:
    :return:
    """
    # TODO: Finish this off - would be great for consistently extracting data from xml
    value = None
    if isinstance(xml_record, str):
        return value
    if term not in xml_record.keys():
        for record in xml_record:
            value = xml_parser(record, term)
            if value:
                return value
            else:
                continue
    else:
        return xml_record[term]
    return value


def read_phylip_to_dict(phylip_input):
    header_dict = dict()
    tmp_seq_dict = dict()
    seq_dict = dict()
    x = 0

    try:
        phylip = open(phylip_input, 'r')
    except IOError:
        logging.error("Unable to open the Phylip file (" + phylip_input + ") provided for reading!\n")
        sys.exit(5)

    line = phylip.readline()
    try:
        num_sequences, aln_length = line.strip().split(' ')
        num_sequences = int(num_sequences)
        aln_length = int(aln_length)
    except ValueError:
        logging.error("Phylip file is not formatted correctly!\n" +
                      "Header must contain 2 space-separated fields (number of sequences and alignment length).\n")
        sys.exit(5)

    line = phylip.readline()
    while line:
        line = line.strip()
        if len(line.split()) == 2:
            # This is the introduction set: header, sequence
            header, sequence = line.split()
            header_dict[x] = header
            tmp_seq_dict[x] = sequence
            x += 1
        elif 60 >= len(line) >= 1:
            tmp_seq_dict[x] += line
            x += 1
        elif line == "":
            # Reset accumulator on blank lines
            x = 0
        else:
            logging.error("Unexpected line in Phylip file:\n" + line + "\n")
            sys.exit(5)
        line = phylip.readline()

        if x > num_sequences:
            logging.error("Accumulator has exceeded the number of sequences in the file (according to header)!\n")
            sys.exit(5)

    # Check that the alignment length matches that in the header line
    if num_sequences != len(tmp_seq_dict):
        logging.error("Number of lines declared in Phylip header (" + str(num_sequences) +
                      ") does not match number of sequences parsed (" + str(len(tmp_seq_dict)) + ")!\n")
        sys.exit(5)

    x = 0
    while x < num_sequences-1:
        if len(tmp_seq_dict[x]) != aln_length:
            logging.error(header_dict[x] +
                          " sequence length exceeds the stated multiple alignment length (according to header)!\n" +
                          "sequence length = " + str(len(tmp_seq_dict[x])) +
                          ", alignment length = " + str(aln_length) + "\n")
            sys.exit(5)
        else:
            pass
        x += 1

    phylip.close()
    for x in header_dict:
        seq_dict[header_dict[x]] = tmp_seq_dict[x]
    return seq_dict


def read_stockholm_to_dict(sto_file):
    """

    :param sto_file: A Stockholm-formatted multiple alignment file
    :return: A dictionary with sequence headers as keys and sequences as values
    """
    seq_dict = dict()

    try:
        sto_handler = open(sto_file, 'r')
    except IOError:
        logging.error("Unable to open " + sto_file + " for reading!\n")
        sys.exit(3)

    line = sto_handler.readline()
    while line:
        line = line.strip()
        if re.match("^[#|/].*", line):
            # Skip the header (first line) as well as secondary structure lines
            pass
        elif not line:
            pass
        else:
            try:
                seq_name, sequence = line.split()
            except ValueError:
                logging.error("Unexpected line format in " + sto_file + ":\n" + line + "\n")
                sys.exit(3)

            if seq_name not in seq_dict:
                seq_dict[seq_name] = ""
            seq_dict[seq_name] += re.sub('\.', '-', sequence.upper())
        line = sto_handler.readline()

    return seq_dict


def read_uc(uc_file):
    """
    Function to read a USEARCH cluster (.uc) file

    :param uc_file: Path to a .uc file produced by USEARCH
    :return: Dictionary where keys are numerical identifiers and values are Cluster objects
        The Cluster object
    """
    cluster_dict = dict()
    rep_len_map = dict()
    try:
        uc = open(uc_file, 'r')
    except IOError:
        logging.error("Unable to open USEARCH cluster file " + uc_file + " for reading!\n")
        sys.exit(13)

    logging.debug("Reading usearch cluster file... ")

    # Find all clusters with multiple identical sequences
    for line in uc:
        cluster_type, num_id, length, identity, _, _, _, cigar, header, representative = line.strip().split("\t")
        if cluster_type == "S":
            cluster_dict[num_id] = Cluster(header)
            rep_len_map[header] = length
        elif cluster_type == "H":
            cluster_dict[num_id].members.append([header, identity])
        elif cluster_type == "C":
            pass
        else:
            logging.error("Unexpected cluster type '" + str(cluster_type) + "' in " + uc_file + "\n")
            sys.exit(13)

    uc.close()
    logging.debug("done.\n")
    return cluster_dict


def read_rpkm(rpkm_output_file):
    """
    Read the CSV file written by rpkm. A header and line with unmapped reads is expected and are skipped.
    Each line is expected to have 4 elements:
     Sample ID, sequence name, number of reads recruited, RPKM

    :param rpkm_output_file: A file path
    :return: Dictionary mapping contig names to floats
    """
    rpkm_values = dict()

    try:
        rpkm_stats = open(rpkm_output_file)
    except IOError:
        logging.error("Unable to open " + rpkm_output_file + " for reading!\n")
        sys.exit(13)

    # Skip the header
    next(rpkm_stats)
    # Skip the line with unaligned reads
    next(rpkm_stats)
    for line in rpkm_stats:
        # Line format is Sample ID (output file name), sequence name, number of reads recruited, RPKM
        try:
            _, seq_name, _, rpkm = line.strip().split(',')
        except ValueError:
            n_values = str(len(line.split(',')))
            logging.error("Unexpected line format in RPKM file - should contain 4 elements, "
                          "" + n_values + " encountered. Offending line:\n" + line + "\n")
            sys.exit(13)
        rpkm_values[seq_name] = float(rpkm)
    rpkm_stats.close()
    return rpkm_values


def validate_alignment_trimming(msa_files: list, unique_ref_headers: set, queries_mapped=False, min_seq_length=30):
    """
    Parse a list of multiple sequence alignment (MSA) files and determine whether the multiple alignment:
        1. is shorter then the min_seq_length (30 by default)
        2. is missing any reference sequences
    The number of query sequences discarded - these may have been added by hmmalign or PaPaRa - is returned via a string

    NOTE: Initially designed for sequence records with numeric names (e.g. >4889) but accomodates other TreeSAPP formats

    :param msa_files: A list of either Phylip or FASTA formatted MSA files
    :param unique_ref_headers: A set of all headers that were in the untrimmed MSA
    :param queries_mapped: Boolean indicating whether sequences should be present in addition to reference sequences.
           While query sequences _could_ be identified as any that are not in unique_ref_headers,
           queries have names that are negative integers for more rapid and scalable identification
    :param min_seq_length: Optional minimum unaligned (no '-'s) length a sequence must exceed to be retained
    :return: 1. Dictionary indexed by MSA file name mapping to FASTA-dictionaries
             2. A string mapping the number of query sequences removed from each MSA file
             3. A string describing the number of sequences discarded
    """
    discarded_seqs_string = ""
    successful_multiple_alignments = dict()
    failed_multiple_alignments = list()
    n_refs = len(unique_ref_headers)
    for multi_align_file in msa_files:
        filtered_multi_align = dict()
        discarded_seqs = list()
        num_queries_retained = 0
        n_retained_refs = 0
        n_msa_refs = 0
        f_ext = multi_align_file.split('.')[-1]

        # Read the multiple alignment file
        if re.search("phy", f_ext):  # File is in Phylip format
            seq_dict = read_phylip_to_dict(multi_align_file)
        elif re.match("^f", f_ext):  # This is meant to match all fasta extensions
            seq_dict = read_fasta_to_dict(multi_align_file)
        elif f_ext == "mfa":  # This is meant to match a multiple alignment in FASTA format
            seq_dict = read_fasta_to_dict(multi_align_file)
        else:
            logging.error("Unable to detect file format of " + multi_align_file + ".\n")
            sys.exit(13)

        # Parse the MSA dict and ensure headers are integer-compatible
        multi_align = dict()
        for seq_name in seq_dict:
            seq = seq_dict[seq_name]
            try:
                if int(seq_name) > 0:
                    n_msa_refs += 1
            except ValueError:
                if re.match(r"^_\d+", seq_name):
                    leaf_num = re.sub("^_", '-', seq_name)
                # The section of regular expresion after '_' needs to match denominator and refpkg names
                elif re.match(r"^\d+_\w{2,10}$", seq_name):
                    leaf_num = seq_name.split('_')[0]
                else:
                    logging.error("Unexpected sequence name '" + seq_name +
                                  "' detected in " + multi_align_file + ".\n")
                    sys.exit(13)
                if int(leaf_num) > 0:
                    n_msa_refs += 1
            multi_align[seq_name] = seq
        if len(multi_align) == 0:
            logging.warning("No sequences were read from " + multi_align_file + ".\n" +
                            "The untrimmed alignment will be used instead.\n")
            failed_multiple_alignments.append(multi_align_file)
            continue
        # The numeric identifiers make it easy to maintain order in the Phylip file by a numerical sort
        for seq_name in sorted(multi_align, key=lambda x: int(x.split('_')[0])):
            seq_dummy = re.sub('-', '', multi_align[seq_name])
            if len(seq_dummy) < min_seq_length:
                discarded_seqs.append(seq_name)
            else:
                filtered_multi_align[seq_name] = multi_align[seq_name]
                # The negative integers indicate this is a query sequence
                if seq_name[0] == '-':
                    num_queries_retained += 1
                else:
                    n_retained_refs += 1
        discarded_seqs_string += "\n\t\t" + multi_align_file + " = " + str(len(discarded_seqs))
        if len(discarded_seqs) == len(multi_align.keys()):
            # Throw an error if the final trimmed alignment is shorter than min_seq_length, and therefore empty
            logging.warning("Multiple sequence alignment in " + multi_align_file +
                            " is shorter than minimum sequence length threshold (" + str(min_seq_length) +
                            ").\nThe untrimmed alignment will be used instead.\n")
            failed_multiple_alignments.append(multi_align_file)
        elif n_refs > n_msa_refs:
            # Testing whether there were more sequences in the untrimmed alignment than the trimmed one
            logging.warning("Reference sequences in " + multi_align_file + " were removed during alignment trimming " +
                            "suggesting either truncated sequences or the initial reference alignment was terrible.\n" +
                            "The untrimmed alignment will be used instead.\n")
            failed_multiple_alignments.append(multi_align_file)
        elif n_refs > n_retained_refs:
            logging.warning("Reference sequences shorter than the minimum character length (" +
                            str(min_seq_length) + ") in " + multi_align_file +
                            " were removed after alignment trimming.\n" +
                            "The untrimmed alignment will be used instead.\n")
            failed_multiple_alignments.append(multi_align_file)
        # Ensure that there is at least 1 query sequence retained after trimming the multiple alignment
        elif queries_mapped and num_queries_retained == 0:
            logging.warning("No query sequences in " + multi_align_file + " were retained after trimming.\n")
        else:
            successful_multiple_alignments[multi_align_file] = filtered_multi_align

        if multi_align_file in successful_multiple_alignments:
            discarded_seqs_string += " (retained)"
        else:
            discarded_seqs_string += " (removed)"

    return successful_multiple_alignments, failed_multiple_alignments, discarded_seqs_string


def multiple_alignment_dimensions(seq_dict, mfa_file):
    """
    Checks to ensure all sequences are the same length and returns a tuple of (nrow, ncolumn)

    :param seq_dict: A dictionary containing headers as keys and sequences as values
    :param mfa_file: The name of the multiple alignment FASTA file being validated
    :return: tuple = (nrow, ncolumn)
    """
    sequence_length = 0
    for seq_name in seq_dict:
        sequence = seq_dict[seq_name]
        if sequence_length == 0:
            sequence_length = len(sequence)
        elif sequence_length != len(sequence) and sequence_length > 0:
            logging.error("Number of aligned columns is inconsistent in " + mfa_file + "!\n")
            sys.exit(3)
        else:
            pass
            # Sequence is the right length, carrying on
    return len(seq_dict), sequence_length


def read_seq_taxa_table(seq_names_to_taxa: str):
    seq_lineage_map = dict()
    try:
        handler = open(seq_names_to_taxa, 'r')
    except IOError:
        logging.error("Unable to open '" + seq_names_to_taxa + "' for reading!\n")
        sys.exit(3)

    for line in handler:
        try:
            seq_name, lineage = line.strip().split("\t")
        except (ValueError, IndexError):
            logging.error("Bad line encountered in '" + seq_names_to_taxa + "' - expected two tab-separated fields:\n" +
                          line)
            sys.exit(3)
        if seq_name[0] == '>':
            seq_name = seq_name[1:]
        seq_lineage_map[seq_name] = lineage
    handler.close()
    return seq_lineage_map
