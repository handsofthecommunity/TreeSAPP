__author__ = 'Connor Morgan-Lang'

import sys
import re
import _tree_parser
import os
from .utilities import Autovivify, mean
from ete3 import Tree
from scipy import log2


def get_node(tree, pos):
    node = ""
    pos += 1
    c = tree[pos]
    while c != '}':
        node += c
        pos += 1
        c = tree[pos]
    return int(node), pos


def map_internal_nodes_leaves(tree):
    """
    Loads a mapping between all nodes (internal and leaves) and all leaves
    :return:
    """
    no_length_tree = re.sub(":[0-9.]+{", ":{", tree)
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


def find_mean_pairwise_distances(children):
    pairwise_dists = list()
    for rleaf in children:
        for qleaf in children:
            if rleaf.name != qleaf.name:
                pairwise_dists.append(rleaf.get_distance(qleaf))
    return sum(pairwise_dists) / len(pairwise_dists)


def get_tip_distances(parent_node):
    children = parent_node.get_leaves()
    distances = [parent_node.get_distance(child) for child in children]
    return distances


def find_cluster(lost_node: Tree, intra_distances: list = []):
    """
    Recursively calculates whether a node in a tree is the root of a cluster,
    where a cluster is defined as a sub-tree (clade) whose members satisfy the condition:
     the distance to the parent of the current subtree's root multiplied by the logarithm base 2 of the number of
     cousins is greater than the mean(intra-cluster root-to-tip distances).
    Adding a large number of members to the clade is penalized, as well as large distances from the existing subtree to
    a new parent.

    :param lost_node: A node within a tree, for which we want to orient
    :param intra_distances: A list with the current set of leaf-tip distances
    :return: Tree node, a list of float distances
    """
    parent = lost_node.up
    if lost_node.is_root() or parent.is_root():
        return lost_node, intra_distances

    if not intra_distances:
        # If this is the initial attempt at finding lost_node's cluster, find the intra-cluster leaf distances
        intra_distances = get_tip_distances(lost_node)

    # Penalty for increasing the size of the clade is log-base 2
    cousins = lost_node.get_sisters()[0].get_leaf_names()
    parent_dist = parent.get_distance(lost_node)
    cost = parent_dist * log2(len(cousins) + 1)

    if mean(intra_distances) > cost:
        return lost_node, intra_distances

    # Add the distance from the parent to the
    intra_distances = [dist+parent_dist for dist in intra_distances]
    for cousin in cousins:
        intra_distances.append(parent.get_distance(cousin))
    return find_cluster(parent, intra_distances)


def subtrees_to_dictionary(subtrees_string, tree_info):
    subtree_list = subtrees_string.split(';')
    for subtree in subtree_list:
        node = subtree.split(')')[-1]
        tree_info['subtree_of_node'][node] = subtree
    return tree_info


def create_tree_info_hash():
    tree_info = Autovivify()
    return tree_info


def create_tree_internal_node_map(tree_string):
    """
    Loads a mapping between all nodes (internal and leaves) and all leaves
    :return:
    """

    def get_node(tree, pos):
        node = ""
        pos += 1
        c = tree[pos]
        while not re.match("\d", c):
            pos += 1
            c = tree[pos]
        while re.match("\d", c):
            node += c
            pos += 1
            c = tree[pos]
        return node, pos

    tree_string = re.sub("-\d+", '-', tree_string)
    internal_node_counter = 0
    node_map = dict()
    node_stack = list()  # This stack handles the sibling nodes
    x = 0
    while x < len(tree_string):
        c = tree_string[x]
        if re.match("\d", c):
            num_buffer, x = get_node(tree_string, x - 1)
            node_map[internal_node_counter] = [str(num_buffer)]
            node_stack.append(internal_node_counter)
            internal_node_counter += 1
            x -= 1
        elif c == ')':
            while c == ')' and x < len(tree_string):
                if tree_string[x + 1] == ';':
                    break
                node_map[internal_node_counter] = node_map[node_stack.pop()] + node_map[node_stack.pop()]
                node_stack.append(internal_node_counter)
                internal_node_counter += 1
                x += 1
                c = tree_string[x]
        x += 1
    return node_map


def format_children_assignments(children_assignments, tree_info):
    children_of_nodes = children_assignments.split(';')
    for family_string in children_of_nodes:
        parent, children = family_string.split('=')
        for node in children.split(','):
            tree_info['children_of_node'][parent][node] = 1
    return tree_info


def format_parent_assignments(parent_assignments, tree_info):
    parents_of_nodes = parent_assignments.split(',')
    for pair in parents_of_nodes:
        node, parent = pair.split(':')
        tree_info['parent_of_node'][node] = parent
    return tree_info


def format_subtrees(subtrees):
    terminal_children_of_reference = Autovivify()
    subtree_list = subtrees.split(',')
    for subtree in subtree_list:
        nodes = subtree.split(' ')
        node_ints = [int(x) for x in nodes]
        sorted_node_strings = [str(i) for i in sorted(node_ints)]
        terminal_children_of_reference[' '.join(sorted_node_strings) + ' '] = 1
    return terminal_children_of_reference


def deconvolute_assignments(reference_tree_assignments):
    tree_info = create_tree_info_hash()
    children_assignments, parent_assignments, subtrees = reference_tree_assignments.strip().split('\n')
    tree_info = format_children_assignments(children_assignments, tree_info)
    tree_info = format_parent_assignments(parent_assignments, tree_info)
    terminal_children_of_reference = format_subtrees(subtrees)
    return tree_info, terminal_children_of_reference


def read_and_map_internal_nodes_from_newick_tree(reference_tree_file, denominator):
    # Using the C++ _tree_parser extension:
    reference_tree_elements = _tree_parser._read_the_reference_tree(reference_tree_file)
    internal_node_map = create_tree_internal_node_map(reference_tree_elements)
    return internal_node_map


def read_and_understand_the_reference_tree(reference_tree_file, denominator):
    # Using the C++ _tree_parser extension:
    reference_tree_elements = _tree_parser._read_the_reference_tree(reference_tree_file)
    reference_tree_assignments = _tree_parser._get_parents_and_children(reference_tree_elements)
    if reference_tree_assignments == "$":
        sys.stderr.write("Poison pill received from " + denominator + "\n")
        sys.stderr.flush()
        return denominator, None
    else:
        reference_tree_info, terminal_children_of_reference = deconvolute_assignments(reference_tree_assignments)
        return denominator, terminal_children_of_reference


def annotate_partition_tree(code_name, fasta_replace_dict, bipart_tree):
    try:
        tree_txt = open(bipart_tree, 'r')
    except IOError:
        raise IOError("Unable to open RAxML bipartition tree " + bipart_tree + " for reading.")

    tree = tree_txt.readline()
    tree_txt.close()
    for mltree_id_key in fasta_replace_dict.keys():
        tree = re.sub('\(' + mltree_id_key + "_" + code_name, '(' + fasta_replace_dict[mltree_id_key].organism, tree)
        tree = re.sub(',' + mltree_id_key + "_" + code_name, ',' + fasta_replace_dict[mltree_id_key].organism, tree)

    tree_output_dir = os.path.dirname(bipart_tree)
    annotated_tree_name = tree_output_dir + os.sep + "RAxML_bipartitions_annotated." + code_name
    try:
        annotated_tree = open(annotated_tree_name, 'w')
    except IOError:
        raise IOError("Unable to open the annotated RAxML tree " + annotated_tree_name + " for writing!")

    annotated_tree.write(tree)
    annotated_tree.close()

    return


def tree_leaf_distances(tree: Tree):
    # max_dist_threshold equals the maximum path length from root to tip in its clade
    leaf_distances = []
    d = 0.0
    max_dist = None
    topology_only = False
    for post, n in tree.iter_prepostorder(is_leaf_fn=None):
        if n is tree:
            continue
        if post:
            d -= n.dist
        else:
            if n.is_leaf():
                total_d = d + n.dist if not topology_only else d
                leaf_distances.append(total_d)
                if max_dist is None or total_d > max_dist:
                    max_dist = total_d
            else:
                d += n.dist if not topology_only else 1.0
    return max_dist, leaf_distances


def index_tree_edges(tree: str):
    edge_index = dict()
    dist = ""
    edge = ""
    i = 0
    n = len(tree)
    while i < n:
        if tree[i] in [':', '{']:
            i += 1
            if dist:
                while re.match(r"[0-9]", tree[i]):
                    edge += tree[i]
                    i += 1
                edge_index[edge] = float(dist)
                dist = ""
                edge = ""
            else:
                while re.match(r"[0-9.]", tree[i]):
                    dist += tree[i]
                    i += 1
        else:
            i += 1

    return edge_index