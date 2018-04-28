#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'


import sys
import argparse
import os
import re
import inspect
import glob
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
sys.path.insert(0, cmd_folder + os.sep + ".." + os.sep)
from entish import get_node
from treesapp import parse_ref_build_params, jplace_parser


def get_arguments():
    parser = argparse.ArgumentParser(add_help=False,
                                     description="This script is generally used for layering extra annotations "
                                                 "beyond taxonomy (such as Subgroup or Metabolism) to TreeSAPP outputs."
                                                 " This is accomplished by adding an extra column (to all rows) of an "
                                                 "existing marker_contig_map.tsv and annotating the relevant sequences")
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("-c", "--colours_style",
                               help="The colours_style file exported from iTOL with the annotation information. "
                                    "For the variable name to be automatically inferred (rather than through `names`). "
                                    "Format of the file should be `marker`_`var`.txt. For example: mcrA_Metabolism.txt "
                                    "would create a new column in marker_contig_map.tsv named 'Metabolism'.",
                               required=True,
                               nargs='+')
    required_args.add_argument("-o", "--output",
                               help="The TreeSAPP output directory.",
                               required=True)

    optopt = parser.add_argument_group("Optional options")
    optopt.add_argument('-n', '--names',
                        help='The names corresponding to each of the colours_style files.'
                             ' Provide a comma-separated list if multiple colours_style files.',
                        required=False,
                        default=None,
                        nargs='+')

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument('-v', '--verbose',
                                    action='store_true',
                                    default=False,
                                    help='Prints a more verbose runtime log')
    miscellaneous_opts.add_argument("-h", "--help",
                                    action="help",
                                    help="Show this help message and exit")

    args = parser.parse_args()
    args.treesapp = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep

    # Adding a dummy value for format_read_fasta args namespace
    args.min_seq_length = 10

    if sys.version_info > (2, 9):
        args.py_version = 3
    else:
        args.py_version = 2

    if not os.path.isdir(args.output):
        os.makedirs(args.output)

    return args


def check_arguments(args):
    """
    Check that the required files (e.g. jplace, marker_contig_map, annotation files) exist
    :param args:
    :return:
    """
    various_dir = args.output + os.sep + "various_outputs" + os.sep
    jplace_files = glob.glob(various_dir + "*jplace")
    if len(jplace_files) == 0:
        sys.stderr.write("ERROR: could not find .jplace files in " + various_dir + "\n")
        sys.exit()
    if not os.path.isfile(os.sep.join([args.output, "final_outputs", "marker_contig_map.tsv"])):
        sys.stderr.write("ERROR: could not find a classification file in " + args.output + os.sep + "final_outputs\n")
        sys.exit()
    for annot_f in args.colours_style:
        if not os.path.isfile(annot_f):
            sys.stderr.write("ERROR: " + annot_f + " does not exist!\n")
            sys.exit()
    return jplace_files


def parse_marker_classification_table(args):
    """
    Function to read marker_contig_map.tsv and gather the relevant information for adding extra annotations
    This function is different from Clade_exclusion_analyzer::read_marker_classification_table(assignment_file)
    as we are interested in all fields in this function.
    :param args:
    :return:
    """
    master_dat = dict()
    field_order = dict()
    marker_classification_file = os.sep.join([args.output, "final_outputs", "marker_contig_map.tsv"])
    try:
        classifications = open(marker_classification_file, 'r')
    except IOError:
        sys.stderr.write("ERROR: Unable to open " + marker_classification_file + " for reading!\n")
        sys.exit()

    header = classifications.readline()
    header_fields = header.strip().split("\t")
    x = 0
    for field in header_fields:
        field_order[x] = field
        x += 1
    line = classifications.readline()
    while line:
        fields = line.strip().split("\t")
        if len(fields) != len(header_fields):
            sys.stderr.write("ERROR: Inconsistent number of columns in table! Offending line:\n" + line + "\n")
            sys.exit()
        master_dat[fields[0]] = dict()
        y = 1
        while y < len(fields):
            master_dat[fields[0]][header_fields[y]] = fields[y]
            y += 1
        line = classifications.readline()

    classifications.close()
    # for query in master_dat:
    #     print(master_dat[query])
    return master_dat, field_order


def read_colours_file(args, annotation_file):
    """
    Read annotation data from 'annotation_file' and store it in marker_subgroups under the appropriate
    marker and data_type.
    :param args:
    :param annotation_file:
    :return: A dictionary of lists where each list is populated by tuples with start and end leaves
    """
    try:
        style_handler = open(annotation_file, 'r')
    except IOError:
        sys.stderr.write("ERROR: Unable to open " + annotation_file + " for reading!\n")
        sys.exit()

    clusters = dict()

    range_line = re.compile("^(\d+)\|(\d+)\srange\s.*\)\s(.*)$")
    single_node = re.compile("^(\d+)\srange\s.*\)\s(.*)$")
    line = style_handler.readline()
    # Skip the header
    while line and not range_line.match(line):
        line = style_handler.readline()
    # Begin parsing the data from 4 columns
    while line:
        style_data = range_line.match(line)
        if style_data:
            description = style_data.group(3)
            if description not in clusters.keys():
                clusters[description] = list()
            clusters[description].append((style_data.group(1), style_data.group(2)))
        elif single_node.match(line):
            node, description = single_node.match(line).groups()
            clusters[description].append((node, node))
        else:
            sys.stderr.write("WARNING: Unrecognized line formatting in " + annotation_file + ":\n")
            sys.stderr.write(line + "\n")
            sys.stderr.write("The annotations in this line will be skipped...\n")
        line = style_handler.readline()
    style_handler.close()

    if args.verbose:
        sys.stdout.write("\tParsed " + str(len(clusters)) +
                         " clades from " + annotation_file + "\n")

    return clusters


def create_node_map(jplace_tree_string):
    """
    Loads a mapping between all nodes (internal and leaves) and all leaves
    :return:
    """
    no_length_tree = re.sub(":[0-9.]+{", ":{", jplace_tree_string)
    node_map = dict()
    node_stack = list()
    leaf_stack = list()
    x = 0
    num_buffer = ""
    while x < len(no_length_tree):
        c = no_length_tree[x]
        if re.search(r"[0-9]", c):
            while re.search(r"[0-9]", c):
                num_buffer += c
                x += 1
                c = no_length_tree[x]
            node_stack.append([str(num_buffer)])
            num_buffer = ""
            x -= 1
        elif c == ':':
            # Append the most recent leaf
            current_node, x = get_node(no_length_tree, x + 1)
            node_map[current_node] = node_stack.pop()
            leaf_stack.append(current_node)
        elif c == ')':
            # Set the child leaves to the leaves of the current node's two children
            while c == ')' and x < len(no_length_tree):
                if no_length_tree[x + 1] == ';':
                    break
                current_node, x = get_node(no_length_tree, x + 2)
                node_map[current_node] = node_map[leaf_stack.pop()] + node_map[leaf_stack.pop()]
                leaf_stack.append(current_node)
                x += 1
                c = no_length_tree[x]
        x += 1
    return node_map


def parse_clades_from_tree(args, jplace_tree_string, clusters):
    """

    :param args:
    :param jplace_tree_string:
    :param clusters: Dictionary with the cluster names for keys and a tuple containing leaf boundaries as values
    :return:
    """
    clade_members = dict()
    leaf_annotation_map = dict()
    # all_annotated_leaves_map = dict()
    leaves_in_clusters = set()

    internal_node_map = create_node_map(jplace_tree_string)

    # TODO: Make this code more efficient because it currently isn't but it works. Everything is small so no pressure.
    # Create a dictionary to map the cluster name (e.g. Function, Activity, Class, etc) to the leaf nodes
    for cluster in clusters.keys():
        if cluster not in leaf_annotation_map:
            leaf_annotation_map[cluster] = list()
            clade_members[cluster] = set()
        for frond_tips in clusters[cluster]:
            start, end = frond_tips
            # Find the minimum set that includes both start and end
            warm_front = dict()
            for inode in internal_node_map:
                clade = internal_node_map[inode]
                if start in clade:
                    if start == end:
                        clade_members[cluster].add(inode)
                    warm_front[len(clade)] = clade
            for size in sorted(warm_front, key=int):
                if end in warm_front[size]:
                    leaf_annotation_map[cluster] += warm_front[size]
                    break

    # # Intersect the clusters for sequences that map near the root
    # for primary_cluster in leaf_annotation_map:
    #     all_annotated_leaves_map[primary_cluster] = leaf_annotation_map[primary_cluster]
    #     for secondary_cluster in leaf_annotation_map:
    #         if primary_cluster == secondary_cluster:
    #             continue
    #         hybrid = sorted([primary_cluster, secondary_cluster])
    #         if ','.join(hybrid) not in all_annotated_leaves_map:
    #             clade_members[','.join(hybrid)] = set()
    #             all_annotated_leaves_map[','.join(hybrid)] = leaf_annotation_map[primary_cluster] + leaf_annotation_map[secondary_cluster]

    # Map the internal nodes (from Jplace tree string) to the cluster names
    for cluster in clade_members:
        for inode in internal_node_map:
            contained = True
            for leaf in internal_node_map[inode]:
                if leaf not in leaf_annotation_map[cluster]:
                    contained = False
                    break
            if contained:
                clade_members[cluster].add(str(inode))
                leaves_in_clusters.add(str(inode))

    if args.verbose:
        sys.stdout.write("\tCaptured " + str(len(leaves_in_clusters)) + " nodes in clusters.\n")

    diff = len(internal_node_map) - len(leaves_in_clusters)
    if diff != 0:
        unannotated = list()
        sys.stderr.write("WARNING: the following internal nodes were not mapped to annotation groups:\n")
        for inode in internal_node_map:
            contained = False
            for cluster in clade_members:
                if str(inode) in clade_members[cluster]:
                    contained = True
                    break
            if not contained:
                unannotated.append(str(inode))
        sys.stderr.write("\t" + ', '.join(unannotated) + "\n")
        sys.stderr.flush()

    return clade_members


def map_queries_to_annotations(marker_tree_info, marker_build_dict, jplace_files_to_parse, master_dat):
    num_unclassified = 0
    for jplace in jplace_files_to_parse:
        file_name = os.path.basename(jplace)
        jplace_info = re.match("RAxML_portableTree.([A-Z][0-9]{4})_(.*).jplace", file_name)
        gene_code, contig_name = jplace_info.groups()
        if contig_name not in master_dat.keys():
            num_unclassified += 1
            continue
        marker = marker_build_dict[gene_code].cog
        for data_type in marker_tree_info:
            if marker in marker_tree_info[data_type]:
                metadata_placement = set()
                itol_datum = jplace_parser(jplace)
                for pquery in itol_datum.placements:
                    for placement in pquery['p']:
                        node = str(placement[0])
                        for group in marker_tree_info[data_type][marker]:
                            if node in marker_tree_info[data_type][marker][group]:
                                metadata_placement.add(group)

                if len(metadata_placement) == 0:
                    metadata_placement.add("Unknown")
                master_dat[contig_name][data_type] = ';'.join(sorted(metadata_placement))
    if num_unclassified > 0:
        sys.stdout.write("Number of placed sequences that were unclassified: " + str(num_unclassified) + "\n")
    return master_dat


def write_classification_table(args, field_order, master_dat):
    """
    Writes data in master_dat to a new tabular file with original and extra annotation information
    :param args:
    :param field_order:
    :param master_dat:
    :return:
    """
    output_file = os.sep.join([args.output, "final_outputs", "extra_annotated_marker_contig_map.tsv"])
    try:
        table_handler = open(output_file, 'w')
    except IOError:
        sys.stderr.write("ERROR: Unable to open " + output_file + " for writing!\n")
        sys.exit()

    fields = list()
    # Prepare the new header and write it to the new classification table
    for order in sorted(field_order.keys()):
        fields.append(field_order[order])
    table_handler.write("\t".join(fields) + "\n")

    # Now parse the classification data in master_dat
    for query_name in master_dat:
        query_classification_list = [query_name]
        field_acc = 1
        while field_acc < len(fields):
            query_classification_list.append(str(master_dat[query_name][field_order[field_acc]]))
            field_acc += 1
        table_handler.write("\t".join(query_classification_list) + "\n")

    table_handler.close()

    return


def main():
    ##
    # Worklow:
    #   1. Read data/tree_data/ref_build_parameters.tsv to get marker codes, denominators, and more (oh my!)
    #   2. Read the marker_contig_map.tsv file from the output directory to create the master data structure
    #   3. For each of the colours_styles files provided (potentially multiple for the same marker):
    #       3.1) Add the annotation variable to master_dat for every sequence (instantiate with "NA")
    #       3.2) Read the .jplace file for every sequence classified as marker
    #       3.3) Add the annotation information to every sequence classified as marker in master_dat
    #   4. Write the new classification file called "extra_annotated_marker_contig_map.tsv"
    ##
    args = get_arguments()
    jplace_files = check_arguments(args)
    var_dir = args.output + os.sep + "various_outputs" + os.sep

    marker_subgroups = dict()
    unique_markers_annotated = set()
    marker_tree_info = dict()
    jplace_files_to_parse = list()
    jplace_tree_strings = dict()
    marker_build_dict = parse_ref_build_params(args)
    master_dat, field_order = parse_marker_classification_table(args)
    # structure of master dat:
    # {"Sequence_1": {"Field1": x, "Field2": y, "Extra": n},
    #  "Sequence_2": {"Field1": i, "Field2": j, "Extra": n}}
    for annot_f in args.colours_style:
        # Determine the marker being annotated
        marker = data_type = ""
        for gene_code in marker_build_dict:
            marker = marker_build_dict[gene_code].cog
            if re.search("^" + re.escape(marker) + "_\w+.txt$", os.path.basename(annot_f)):
                data_type = re.search("^" + re.escape(marker) + "_(\w+).txt$", os.path.basename(annot_f)).group(1)
                unique_markers_annotated.add(gene_code)
                break
            else:
                marker = data_type = ""
        if marker and data_type:
            if data_type not in marker_subgroups:
                marker_subgroups[data_type] = dict()
            marker_subgroups[data_type][marker] = read_colours_file(args, annot_f)
        else:
            sys.stderr.write("ERROR: Unable to parse the marker and/or annotation type from " + annot_f + "\n")
            sys.exit()
    # Instantiate every query sequence in marker_contig_map with an empty string for each data_type
    for data_type in marker_subgroups:
        for query_seq in master_dat:
            master_dat[query_seq][data_type] = "NA"
    # Update the field_order dictionary with new fields
    field_acc = len(field_order)
    for new_datum in sorted(marker_subgroups.keys()):
        field_order[field_acc] = new_datum
        field_acc += 1

    # Load the query sequence annotations
    for data_type in marker_subgroups:
        if data_type not in marker_tree_info:
            marker_tree_info[data_type] = dict()
        for gene_code in unique_markers_annotated:
            for jplace in jplace_files:
                if re.match(re.escape(var_dir) + "RAxML_portableTree." + re.escape(gene_code) + "_.*.jplace", jplace):
                    jplace_files_to_parse.append(jplace)
                if gene_code not in jplace_tree_strings:
                    jplace_tree_strings[gene_code] = jplace_parser(jplace).tree

            marker = marker_build_dict[gene_code].cog
            if marker in marker_subgroups[data_type]:
                marker_tree_info[data_type][marker] = parse_clades_from_tree(args,
                                                                             jplace_tree_strings[gene_code],
                                                                             marker_subgroups[data_type][marker])
            else:
                pass
    marker_subgroups.clear()
    master_dat = map_queries_to_annotations(marker_tree_info, marker_build_dict, jplace_files_to_parse, master_dat)
    write_classification_table(args, field_order, master_dat)


main()
