#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import argparse
import sys
import os
import inspect
import shutil
import glob
from numpy import sqrt

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
sys.path.insert(0, cmd_folder + os.sep + ".." + os.sep)
from fasta import get_headers
from external_command_interface import launch_write_command
import file_parsers
from classy import prep_logging, ReferencePackage
from entrez_utils import *
from lca_calculations import compute_taxonomic_distance, all_possible_assignments, \
    optimal_taxonomic_assignment, grab_graftm_taxa
from utilities import fish_refpkg_from_build_params


class ClassifiedSequence:
    def __init__(self, header):
        self.name = header
        self.ref = ""
        self.ncbi_tax = ""
        self.assigned_lineage = ""
        self.true_lineage = ""
        self.tax_dist = 0


class ConfusionTest:
    def __init__(self, gene_list):
        self._MAX_TAX_DIST = -1
        self.header_regex = None
        self.ref_packages = {key: ReferencePackage() for key in gene_list}
        self.fn = {key: [] for key in gene_list}
        self.fp = {key: [] for key in gene_list}
        self.tp = {key: [] for key in gene_list}  # This will be a list of ClassifiedSequence instances
        self.tax_lineage_map = dict()
        self.dist_wise_tp = dict()
        self.num_total_queries = 0

    def get_info(self, verbose=False):
        info_string = "Reference packages being tested:\n"
        info_string += ", ".join(list(self.ref_packages.keys())) + "\n"

        self.check_dist()
        info_string += "Stats based on taxonomic distance < " + str(self._MAX_TAX_DIST) + "\n"

        if self.num_total_queries > 0:
            info_string += str(self.num_total_queries) + " query sequences being used for testing.\n"

        if verbose:
            for refpkg in self.ref_packages:
                info_string += self.marker_classification_summary(refpkg)

        return info_string

    def marker_classification_summary(self, refpkg_name):
        """
        Provide a classification summary for a specific marker gene, refpkg_name
        :param refpkg_name:
        :return: A string summarizing the classification performance of a single marker/refpkg_name
        """
        self.check_dist()
        self.check_refpkg_name(refpkg_name)
        num_tp, remainder = self.get_true_positives_at_dist(refpkg_name)

        summary_string = "\nSummary for reference package '" + str(refpkg_name) + "':\n"
        summary_string += "\tTrue positives\t\t" + str(len(self.tp[refpkg_name])) + "\n"
        summary_string += "Stats based on taxonomic distance <" + str(self._MAX_TAX_DIST) + ":\n"
        summary_string += "\tTrue positives\t\t" + str(num_tp) + "\n"
        summary_string += "\tFalse positives\t\t" + str(self.get_false_positives(refpkg_name) + remainder) + "\n"
        summary_string += "\tFalse negatives\t\t" + str(self.get_false_negatives(refpkg_name)) + "\n"
        summary_string += "\tTrue negatives\t\t" + str(self.get_true_negatives(refpkg_name)) + "\n"
        return summary_string

    def retrieve_lineages(self, group=1):
        if not self.header_regex:
            logging.error("Unable to parse taxonomic identifiers from header without a regular expression.\n")
            sys.exit(19)

        # Gather the unique taxonomy IDs and store in EntrezRecord instances
        entrez_records = list()
        for marker in self.tp:
            for tp_seq in self.tp[marker]:  # type: ClassifiedSequence
                header = tp_seq.name
                try:
                    tax_id = self.header_regex.search(header).group(group)
                except (KeyError, TypeError):
                    logging.error("Header '" + str(header) + "' in test FASTA doesn't match the supported format.\n")
                    sys.exit(5)

                # Only make Entrez records for new NCBI taxonomy IDs
                if tax_id not in self.tax_lineage_map:
                    e_record = EntrezRecord(header, "")
                    e_record.ncbi_tax = tax_id
                    e_record.bitflag = 3
                    entrez_records.append(e_record)

        # Query the Entrez database for these unique taxonomy IDs
        fetch_lineages_from_taxids(entrez_records)

        for e_record in entrez_records:  # type: EntrezRecord
            self.tax_lineage_map[e_record.ncbi_tax] = clean_lineage_string(e_record.lineage)
        return

    def bin_true_positives_by_taxdist(self, refpkg_name=None):
        """
        Defines the number of true positives at each taxonomic distance x where 0 <= x <= 7,
        since there are 7 ranks in the NCBI taxonomic hierarchy.
        All sequences correctly classified (at the gene level) are assigned a taxonomic distance,
        so the sum of dist_wise_tp[x] for all x will equal the number of all true positives.

        :return: None
        """
        if not refpkg_name:
            marker_set = self.tp
        else:
            marker_set = [refpkg_name]
        for marker in marker_set:
            self.dist_wise_tp[marker] = dict()
            for tp_inst in self.tp[marker]:  # type: ClassifiedSequence
                # Find the optimal taxonomic assignment
                optimal_taxon = optimal_taxonomic_assignment(self.ref_packages[marker].taxa_trie,
                                                             self.tax_lineage_map[tp_inst.ncbi_tax])
                tp_inst.tax_dist, status = compute_taxonomic_distance(tp_inst.assigned_lineage, optimal_taxon)
                if status > 0:
                    logging.debug("Lineages didn't converge between:\n" +
                                  tp_inst.assigned_lineage + "\n" +
                                  self.tax_lineage_map[tp_inst.ncbi_tax] + "\n")
                try:
                    self.dist_wise_tp[marker][tp_inst.tax_dist].append(tp_inst.name)
                except KeyError:
                    self.dist_wise_tp[marker][tp_inst.tax_dist] = [tp_inst.name]
        return

    def check_dist(self):
        if self._MAX_TAX_DIST < 0:
            logging.error("ConfusionTest's _MAX_TAX_DIST has yet to be set.\n")
            sys.exit(5)
        return

    def check_refpkg_name(self, refpkg_name):
        if refpkg_name not in self.ref_packages:
            logging.error(refpkg_name + " is not found in the names of markers to be tested.\n")
            sys.exit(9)
        return

    def get_true_positives_at_dist(self, refpkg_name=None):
        """
        Calculates the sum of all true positives at a specified maximum taxonomic distance and less.
        Sequences classified at a distance greater than self._MAX_TAX_DIST are counted as false negatives,
        since it is as if the sequences were not classified at all.

        :return: The sum of true positives at taxonomic distance <= max_distance
        """
        self.check_dist()
        all_tp_headers = set()
        remainder_headers = set()
        if refpkg_name:
            marker_set = [refpkg_name]
        else:
            marker_set = self.dist_wise_tp
        for ref_name in marker_set:
            for tax_dist in sorted(self.dist_wise_tp[ref_name], key=int):  # type: int
                if tax_dist <= self._MAX_TAX_DIST:
                    all_tp_headers.update(set(self.dist_wise_tp[ref_name][tax_dist]))
                else:
                    remainder_headers.update(set(self.dist_wise_tp[ref_name][tax_dist]))
        return len(all_tp_headers), len(remainder_headers)

    def get_true_negatives(self, refpkg_name=None):
        acc = 0
        if refpkg_name:
            marker_set = [refpkg_name]
        else:
            marker_set = self.ref_packages
        for marker in marker_set:
            acc += len(self.fn[marker])
            acc += len(self.fp[marker])
            acc += len(self.tp[marker])
        return self.num_total_queries - acc

    def get_false_positives(self, refpkg_name=None):
        if refpkg_name:
            return len(self.fp[refpkg_name])
        else:
            unique_fp = set()
            return len([unique_fp.update(self.fp[marker]) for marker in self.ref_packages])

    def get_false_negatives(self, refpkg_name=None):
        if refpkg_name:
            return len(self.fn[refpkg_name])
        else:
            unique_fn = set()
            return len([unique_fn.update(self.fn[marker]) for marker in self.ref_packages])

    def bin_headers(self, test_seq_names, assignments, annot_map, marker_build_dict):
        """
        Function for sorting/binning the classified sequences at T/F positives/negatives based on the
        :param test_seq_names: List of all headers in the input FASTA file
        :param assignments: Dictionary mapping taxonomic lineages to a list of headers that were classified to this lineage
        :param annot_map: Dictionary mapping reference package (gene) name keys to database names values
        :param marker_build_dict:
        :return: None
        """
        # False positives: those that do not belong to the annotation matching a reference package name
        # False negatives: those whose annotations match a reference package name and were not classified
        # True positives: those with a matching annotation and reference package name and were classified
        # True negatives: those that were not classified and should not have been
        mapping_dict = dict()
        positive_queries = dict()

        for refpkg in annot_map:
            marker = marker_build_dict[refpkg].cog
            orthos = annot_map[refpkg]  # List of all orthologous genes corresponding to a reference package
            for gene in orthos:
                if gene not in mapping_dict:
                    mapping_dict[gene] = []
                mapping_dict[gene].append(marker)

        logging.info("Labelling true test sequences... ")
        for header in test_seq_names:
            try:
                ref_g, tax_id = self.header_regex.match(header).groups()
            except (TypeError, KeyError):
                logging.error("Header '" + header + "' in test FASTA file does not match the supported format.\n")
                sys.exit(5)
            # Keep the name in
            if ref_g in mapping_dict:
                markers = mapping_dict[ref_g]
                ##
                # TODO: This leads to double-counting and therefore needs to be deduplicated later
                ##
                for marker in markers:
                    if marker not in positive_queries:
                        positive_queries[marker] = []
                    positive_queries[marker].append(header)
        logging.info("done.\n")

        logging.info("Assigning test sequences to the four class conditions... ")
        for marker in assignments:
            positives = set(positive_queries[marker])
            true_positives = set()
            refpkg = fish_refpkg_from_build_params(marker, marker_build_dict).denominator
            for tax_lin in assignments[marker]:
                classified_seqs = assignments[marker][tax_lin]
                for seq_name in classified_seqs:
                    try:
                        ref_g, tax_id = self.header_regex.match(seq_name).groups()
                    except (TypeError, KeyError, AttributeError):
                        logging.error(
                            "Classified sequence name '" + seq_name + "' does not match the supported format.\n")
                        sys.exit(5)
                    if ref_g in mapping_dict and marker in mapping_dict[ref_g]:
                        # Populate the relevant information for the classified sequence
                        tp_inst = ClassifiedSequence(seq_name)
                        tp_inst.ncbi_tax = tax_id
                        tp_inst.ref = marker
                        tp_inst.assigned_lineage = tax_lin

                        # Add the True Positive to the relevant collections
                        self.tp[refpkg].append(tp_inst)
                        true_positives.add(seq_name)
                    else:
                        self.fp[refpkg].append(seq_name)

            # Identify the False Negatives using set difference - those that were not classified but should have been
            self.fn[refpkg] = list(positives.difference(true_positives))
        logging.info("done.\n")
        return


def get_arguments():
    parser = argparse.ArgumentParser(add_help=False,
                                     description="A wrapper script for calculating Matthews correlation coefficient for"
                                                 "TreeSAPP or GraftM. Currently only supports testing proteins.")
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("--fasta_input", required=True,
                               help="FASTA-formatted file used for testing the classifiers")
    # required_args.add_argument("--reference_marker", required=True,
    #                            help="Short-form name of the marker gene to be tested (e.g. mcrA, pmoA, nosZ)")
    required_args.add_argument("--annot_map", required=True,
                               help="Path to a tabular file mapping markers being tested to their database annotations."
                                    " First column is the ")

    optopt = parser.add_argument_group("Optional options")
    optopt.add_argument("--tool", default="treesapp", required=False,
                        choices=["treesapp", "graftm", "diamond"],
                        help="Classify using one of the tools: treesapp [DEFAULT], graftm, or diamond.")
    optopt.add_argument("--gpkg_dir", required=False, default=None,
                        help="Path to a directory containing all GraftM reference packages to test."
                             " Files must follow 'name.gpkg' scheme and 'name' is in the first column of --annot_map")
    optopt.add_argument("--output", required=False, default="./MCC_output/",
                        help="Path to a directory for writing output files")
    optopt.add_argument("-p", "--pkg_path", required=False, default=None,
                        help="The path to the TreeSAPP-formatted reference package(s) [ DEFAULT = TreeSAPP/data/ ].")
    # optopt.add_argument('-m', '--molecule',
    #                     help='the type of input sequences (prot = Protein [DEFAULT]; dna = Nucleotide; rrna = rRNA)',
    #                     default='prot',
    #                     choices=['prot', 'dna', 'rrna'])
    # optopt.add_argument("-l", "--length",
    #                     required=False, type=int, default=0,
    #                     help="Arbitrarily slice the input sequences to this length. "
    #                          "Useful for testing classification accuracy for fragments.")

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument("-T", "--num_threads", default=4, required=False,
                                    help="The number of threads to use when running either TreeSAPP or GraftM")
    miscellaneous_opts.add_argument('--overwrite', action='store_true', default=False,
                                    help='overwrites previously processed output folders')
    miscellaneous_opts.add_argument('-v', '--verbose', action='store_true', default=False,
                                    help='Prints a more verbose runtime log')
    miscellaneous_opts.add_argument("-h", "--help",
                                    action="help",
                                    help="Show this help message and exit")

    args = parser.parse_args()

    if args.output[-1] != os.sep:
        args.output += os.sep
    args.treesapp = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep + ".." + os.sep
    if args.tool == "treesapp" and not args.pkg_path:
        args.pkg_path = args.treesapp + "data" + os.sep
    args.targets = ["ALL"]

    if sys.version_info > (2, 9):
        args.py_version = 3
    else:
        args.py_version = 2

    if args.gpkg_dir and args.gpkg_dir[-1] != os.sep:
        args.gpkf_dir += os.sep
    return args


def validate_command(args):
    if args.overwrite:
        if os.path.exists(args.output):
            shutil.rmtree(args.output)

    if args.tool in ["diamond", "graftm"]:
        if not args.gpkg_dir:
            logging.error(args.tool + " specified but not GraftM reference package directory specified.\n")
            sys.exit(17)
        elif not os.path.isdir(args.gpkg_dir):
            logging.error(args.gpkg_dir + " GraftM reference package directory does not exist!\n")
            sys.exit(17)
        elif len(glob.glob(args.gpkg_dir + "*gpkg")) == 0:
            logging.error("Not GraftM reference packages found in " + args.gpkg_dir + ".\n")
            sys.exit(17)
    elif args.tool == "treesapp" and args.gpkg_dir:
        logging.warning("--gpkg specific but tool selected is 'treesapp'... it will be ignored.\n")

    return


def read_annotation_mapping_file(annot_map_file):
    annot_map = dict()
    try:
        annot_map_handler = open(annot_map_file)
    except IOError:
        logging.error("Unable to open " + annot_map_file + " for reading!\n")
        sys.exit(3)

    # Assuming the first column is the reference package name and the second is the database annotation name
    n = 0
    for line in annot_map_handler:
        n += 1
        if line[0] == '#':
            continue
        elif not line:
            continue
        else:
            fields = line.split("\t")
            try:
                annot_map[fields[0]] = fields[1].split(',')
            except KeyError:
                logging.error("Insufficient number of fields on line " + str(n) + " in " + annot_map_file + "!\n" +
                              "File must have the reference package name and the database name in"
                              " the first two columns, respectively. Any number of columns can follow.\n")
                sys.exit(9)

    annot_map_handler.close()
    return annot_map


def calculate_matthews_correlation_coefficient(tp: int, fp: int, fn: int, tn: int):
    numerator = float((tp * tn) - (fp * fn))
    denominator = float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return round(numerator/sqrt(denominator), 3)


def main():
    args = get_arguments()
    log_name = args.output + os.sep + "MCC_log.txt"
    prep_logging(log_name, args.verbose)
    validate_command(args)

    ##
    # Read the file mapping reference package name to the database annotations
    ##
    pkg_name_dict = read_annotation_mapping_file(args.annot_map)
    marker_build_dict = file_parsers.parse_ref_build_params(args)
    test_obj = ConfusionTest(pkg_name_dict.keys())

    ##
    # Load the taxonomic trie for each reference package
    ##
    if args.tool == "treesapp":
        for pkg_name in test_obj.ref_packages:
            refpkg = test_obj.ref_packages[pkg_name]
            marker = marker_build_dict[pkg_name].cog
            refpkg.gather_package_files(marker, args.pkg_path)
            test_obj.ref_packages[pkg_name].taxa_trie = all_possible_assignments(test_obj.ref_packages[pkg_name].lineage_ids)
    else:
        for gpkg in glob.glob(args.gpkg_dir + "*gpkg"):
            pkg_name = str(os.path.basename(gpkg).split('.')[0])
            tax_ids_file = gpkg + os.sep + pkg_name + "_taxonomy.csv"
            test_obj.ref_packages[pkg_name].taxonomic_tree = grab_graftm_taxa(tax_ids_file)

    ##
    # Run the specified taxonomic analysis tool and collect the classifications
    ##
    assignments = {}
    test_fa_prefix = '.'.join(args.fasta_input.split('.')[:-1])
    if args.tool == "treesapp":
        ref_pkgs = ','.join(pkg_name_dict.keys())
        classification_table = os.sep.join([args.output, "TreeSAPP_output", "final_outputs", "marker_contig_map.tsv"])
        if not os.path.isfile(classification_table):
            # TODO: Replace with a call to the treesapp main function
            classify_call = [args.treesapp + "treesapp.py", "-i", args.fasta_input,
                             "-t", ref_pkgs, "-T", str(args.num_threads),
                             "-m", "prot", "--output", args.output + "TreeSAPP_output" + os.sep,
                             "--trim_align", "--overwrite"]
            launch_write_command(classify_call, False)
        classification_lines = file_parsers.read_marker_classification_table(classification_table)
        assignments = file_parsers.parse_assignments(classification_lines)
    else:
        # Since you are only able to analyze a single reference package at a time with GraftM, this is ran iteratively
        for gpkg in glob.glob(args.gpkg_dir + "*gpkg"):
            pkg_name = str(os.path.basename(gpkg).split('.')[0])
            if pkg_name not in pkg_name_dict:
                logging.warning("'" + pkg_name + "' not in " + args.annot_map + " and will be skipped...\n")
                continue
            output_dir = os.sep.join([args.output, "GraftM_output", pkg_name]) + os.sep
            classification_table = output_dir + test_fa_prefix + os.sep + test_fa_prefix + "_read_tax.tsv"
            if not os.path.isfile(classification_table):
                classify_call = ["graftM", "graft",
                                 "--forward", args.fasta_input,
                                 "--graftm_package", gpkg,
                                 "--input_sequence_type", "aminoacid",
                                 "--threads", str(args.num_threads),
                                 "--output_directory", output_dir,
                                 "--force"]
                launch_write_command(classify_call, False)

            assignments[pkg_name] = file_parsers.read_graftm_classifications(classification_table)

    if len(assignments) == 0:
        logging.error("No sequences were classified by " + args.tool + ".\n")
        sys.exit(3)

    logging.info("Reading headers in " + args.fasta_input + "... ")
    test_seq_names = get_headers(args.fasta_input)
    logging.info("done.\n")
    test_obj.num_total_queries = len(test_seq_names)
    eggnog_re = re.compile(r"^>?(COG[A-Z0-9]+|ENOG[A-Z0-9]+)_(\d+)\..*")
    test_obj.header_regex = eggnog_re

    ##
    # Bin the test sequence names into their respective confusion categories (TP, TN, FP, FN)
    ##
    test_obj.bin_headers(test_seq_names, assignments, pkg_name_dict, marker_build_dict)
    test_seq_names.clear()

    ##
    # Parse the taxonomic IDs from EggNOG headers and map taxonomic lineage information to classified sequences
    ##
    _TAXID_GROUP = 2
    test_obj.retrieve_lineages(_TAXID_GROUP)
    test_obj.bin_true_positives_by_taxdist()

    # test_obj._MAX_TAX_DIST = 2
    # print(test_obj.get_info(True))
    ##
    # Report the MCC score across different taxonomic distances - should increase with greater allowed distance
    ##
    d = 0
    mcc_string = "Tax.dist\tMCC\tTrue.Pos\tTrue.Neg\tFalse.Pos\tFalse.Neg\n"
    while d < 7:
        test_obj._MAX_TAX_DIST = d
        num_tp, remainder = test_obj.get_true_positives_at_dist()
        num_fp = test_obj.get_false_positives() + remainder
        num_fn = test_obj.get_false_negatives()
        num_tn = test_obj.get_true_negatives()
        mcc = calculate_matthews_correlation_coefficient(num_tp, num_fp, num_fn, num_tn)
        mcc_string += "\t".join([str(x) for x in [d, mcc, num_tp, num_tn, num_fp, num_fn]]) + "\n"
        d += 1
    logging.info(mcc_string)
    return


main()
