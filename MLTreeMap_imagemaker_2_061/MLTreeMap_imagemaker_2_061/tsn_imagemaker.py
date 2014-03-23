#!/usr/bin/python

# MLTreeMap Imagemaker TSN v. 0.0

# Include the necessary libraries
try:
    import sys
    import math
    import svgwrite
    import argparse
    import os
    import re
    import shutil
except:
    sys.exit('Could not load some user defined module functions')

class Autovivify(dict):
    """
    In cases of Autovivify objects, enable the referencing of variables 
    (and sub-variables) without explicitly declaring those variables 
    beforehand.
    """


    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def os_type():
    """Return the operating system of the user."""
    x = sys.platform
    if x:
        hits = re.search(r'darwin', x, re.I)
        if hits:
            return 'mac'
        hits = re.search(r'win', x, re.I)
        if hits:
            return 'win'
        hits = re.search(r'linux', x, re.I)
        if hits:
            return 'linux'


def pathDelim():
    """Return the path deliminator based on the user's operating system."""
    ostype = os_type()
    if ostype == 'win':
        return '\\'
    if ostype in ['linux', 'mac']:
        return '/'


PATHDELIM = str(pathDelim())


def cmp_to_key(mycmp, num_children):
    """In case of Python before version 2.7"""
    class K(object):
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj, num_children) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj, num_children) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj, num_children) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj, num_children) <= 0
        def __ge__(self, other):
            return mycmp(self.obj, other.obj, num_children) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj, num_children) != 0
    return K


def optimize_node_recursively_sort_helper(a, b, nr_terminal_children_per_child):
    """
    The original optimize_node_recursively uses a complex sort.
    To imitate its function, I chose to make this helper method.
    """
    sort_value_a = int(nr_terminal_children_per_child[a])
    sort_value_b = int(nr_terminal_children_per_child[b])
    if sort_value_a == sort_value_b:
        sort_value_a = int(a)
        sort_value_b = int(b)
    return cmp(sort_value_a,sort_value_b)


def getParser():
    """TK"""
    parser = argparse.ArgumentParser(description='example input:\n' + \
      'python tsn_imagemaker.py -i example_input/\n\n' + \
      'mandatory input option:\n' + \
      '-i: either a) a path to a directory with input files ' + \
      'or b) to a specific input file.\n' + \
      'note that in case a) all input files of the same type ' + \
      '(i.e. p, g, n, h, r, m, d) will be concatenated.\n\n' + \
      'optional input options:\n' + \
      '-o: output directory (default: output/)\n' + \
      '-d: use different colors for different datasets ' + \
      '(0: don\'t use this mode (default),' + \
      '1: use this mode,' + \
      '11: use this mode but concatenate percent values ' + \
      'over datasets). ' + \
      'See the user guide for a more detailed description.\n' + \
      '-t: display trees inclusive text and labels ' + \
      '(1, default) ' + \
      'or without text (0)\n' + \
      '-b: parameter for the size of the placement bubbles ' + \
      '(default: 1)\n' + \
      '-r: display 16s and 18s rRNA hits in different trees ' + \
      '(2, default) ' + \
      'or one tree of life (1)\n' + \
      '-h: high quality placement bubbles (1, default) or ' + \
      'low quality bubbles (0). ' + \
      'High quality bubbles work perfectly with inkscape ' + \
      'but can cause trouble with Adobe Illustrator.')
    parser.add_argument('-i', '--input', required=True, help='either \
a) a path to a directory with input files or \
b) to a specific input file')
    parser.add_argument('-o', '--output', default='output', help='output \
directory (default: output/)')
    parser.add_argument('-d', '--colors', default=0, choices=[0,1,11], \
                        help='use different colors for different datasets \
(0: don\'t use this mode (default),\
1: use this mode,\
11: use this mode but concatenate percent values \
over datasets)')
    parser.add_argument('-t', '--text', default=1, choices=[0,1], \
                        help='display trees inclusive text and labels \
(1, default) or without text(0)')
    parser.add_argument('-b', '--bubble_size', default=1, type=int, \
                        help='parameter for the size of the placement \
bubbles (default: 1)')
    parser.add_argument('-r', '--rRNA', default=2, choices=[1,2], \
                        help='display 16s and 18s rRNA hits in different \
trees (2, default) or one tree of life (1)')
    parser.add_argument('-hq', '--hq_bubbles', default=1, choices=[0,1], \
                        help='high quality bubbles work perfectly with inkscape but can cause trouble with Adobe Illustrator')
    return parser

def read_user_input(parser):
    """TK"""
    args = parser.parse_args()

    # Remove any pre-existing path delim on the output directory
    while re.search(r'/\Z', args.output) or re.search(r'\\\Z', args.output):
        args.output = args.output[:-1]

    # Add the path delim to the output directory
    args.output += PATHDELIM

    # Check to see if output directory already exists;
    # if not, create output directory
    while os.path.isdir(args.output):
        print 'WARNING: Your output directory ' + str(args.output) + \
              'already exists! Overwrite [1], quit [2], or \
change directory [3]?'
        answer = int(raw_input())

        while not answer in [1,2,3]:
            answer = int(raw_input('Invalid input. Please choose \
1, 2, or 3.\n'))

        if answer == 1:
            print 'Do you really want to overwrite the old output directory?'
            print 'All data in it will be lost!'
            answer2 = raw_input('Yes [y] or no [n]?\n')

            while not answer2 in ['y', 'n']:
                answer2 = raw_input('Invalid input. Please choose y or n.\n')

            if answer2 == 'y':
                shutil.rmtree(args.output)
            else:
                sys.exit('Exit MLTreeMap Imagemaker\n')
        elif answer == 2:
            sys.exit('Exit MLTreeMap Imagemaker\n')
        else:
            args.output = raw_input('Please enter the path to the new \
directory.\n')

            while re.search(r'/\Z', args.output) \
                  or re.search(r'\\\Z', args.output):
                args.output = args.output[:-1]

            args.output += PATHDELIM

    os.makedirs(args.output)
    return args


def get_input_files(args):
    input_files = Autovivify()
    input = args.input
    rRNA_display_option = args.rRNA
    input2 = args.input

    # Ensure input file name is precursed by the path delim
    while re.search(r'\A/', input2) \
          or re.search(r'\A\\', input2):
        input2 = input2[:1]

    #input2 = input2

    # If input is an individual file, get file name
    if re.search(r'.*/((.)_.*RAxML_.+.txt)\Z', input2):
        if os.path.isfile(input2):
            sys.exit('ERROR: You entered a file as input (' + str(input2) + '), but it does not exist!!\n')
        temp = re.search(r'.*/((.)_.*RAxML_.+.txt)\Z', input2)
        filename = temp.group(1)
        denominator = temp.group(2)
        if denominator in ['b','a'] and rRNA_display_option == 1:
            denominator = 'A'
        input_files[denominator][input] = filename

    # Else, if input is a directory, get file names of contents
    else:
        try:
            _files = os.listdir(input2)
            files = []
            for _file in _files:
                files.append(_file)
        except OSError:
            sys.exit('ERROR: Your input directory ' + str(input2) + \
                     ' does not exist or contains no relevant files!\n')

        for filename in files:
            if re.search(r'\A(.)_.*RAxML_.+.txt\Z', filename):
                temp = re.search(r'\A(.)_.*RAxML_.+.txt\Z', filename)
                denominator = temp.group(1)
                filename_long = str(input) + str(filename)
                if denominator in ['b','a'] and rRNA_display_option == 1:
                    denominator = 'A'
                input_files[denominator][filename_long] = filename

    return input_files


def concatenate_files(args, input_files):
    """TK"""
    output_dir = args.output
    rRNA_display_option = args.rRNA
    concatenated_input_files = Autovivify()
    text_of_denominator = Autovivify()
    colors = []
    used_colors = Autovivify()
    try:
        input = open('tree_data/available_dataset_colors.txt', 'r')
    except IOError:
        sys.exit('ERROR: tree_data/available_dataset_colors.txt does not exist!\n')
    count = 0
    flag = 0

    for line in input:
        if flag == 0:
            flag = 1
            continue
        line = line.strip()
        red, green, blue = line.split(' ')
        try:
            blue
        except NameError:
            continue
        color_tag = 'rgb(' + str(red) + ', ' + str(green) + ', ' + str(blue) + ')'
        colors.insert(count,color_tag)
        count += 1

    nr_of_colors = len(colors)

    for denominator in sorted(input_files):
        output_filename = str(denominator) + '_concatenated_RAxML_outputs.txt'
        output_filename_long = str(output_dir) + str(output_filename)
        try:
            output = open(output_filename_long, 'w')
        except IOError:
            sys.exit('ERROR: Can\'t create ' + str(output_filename_long) + '\n')
        if args.colors > 0:
            try:
                output2 = open(str(output_dir) + str(denominator) + '_color_legend.txt', 'w')
            except IOError:
                sys.exit('ERROR: Can\'t create the color legend file\n')
        concatenated_input_files[denominator][output_filename_long] = output_filename
        percentages_of_texts = Autovivify()
        nr_of_files = 0

        for filename_long in sorted(input_files[denominator]):
            nr_of_files += 1
            if nr_of_files > nr_of_colors and args.colors > 0:
                sys.exit('ERROR: This directory contains more datasets than can be displayed with \
different colors!\n')
            attachment = 'colorcode_' + str(colors[0])
            used_colors[colors[0]] = 1
            if args.colors > 0:
                attachment = 'colorcode_' + str(colors[nr_of_files-1])
                used_colors[colors[nr_of_files-1]] = 1
                output2.write(str(colors[nr_of_files-1]) + '\t' + str(filename_long) + '\n')
            try:
                file = open(filename_long, 'r')
            except IOError:
                sys.exit('ERROR: Can\'t open ' + str(filename_long))

            for fileLine in file:
                fileLine = fileLine.rstrip()
                weight = ''
                text = ''
                if re.search('\APlacement weight (.+)%:(.+)', fileLine):
                    weight = re.search('\APlacement weight (.+)%:(.+)', fileLine).group(1)
                    text = re.search('\APlacement weight (.+)%:(.+)', fileLine).group(2)
                    text += attachment
                elif re.search('\A# .+ analysis, (.+):', fileLine):
                    text_of_denominator[denominator] = re.search('\A# .+ analysis, (.+):', \
                                                       fileLine).group(1)
                    if re.search('rRNA', re.search('\A# .+ analysis, (.+):', fileLine).group(1)) and \
                       rRNA_display_option == 1:
                        text_of_denominator[denominator] = '16s rRNA & 18s rRNA'
                    continue
                elif re.search('\A# Phylogenetic analysis(.*)', fileLine):
                    text_of_denominator[denominator] = 'MLTreeMap tree of life'
                    if re.search('GEBA', re.search('\A# Phylogenetic analysis(.*)', fileLine)):
                        text_of_denominator[denominator] = 'GEBA tree of life'
                    if re.search('fungi', re.search('\A# Phylogenetic analysis(.*)', fileLine)):
                        text_of_denominator[denominator] = 'Fungi'
                    continue
                else:
                    continue
                if percentages_of_texts[text]:
                    percentages_of_texts[text] += int(weight)
                else:
                    percentages_of_texts[text] = int(weight)
            file.close()

        if args.colors > 0:
            output2.close()
        check = 0

        for text in sorted(percentages_of_texts):
            weight = percentages_of_texts[text]
            relative_weight = weight / nr_of_files
            relative_weight = (int(relative_weight * 10000 + 0.5)) / 10000
            output.write('Placement weight ' + str(relative_weight) + '%:' \
                         + str(text) + '\n')
            check += relative_weight

        print str(denominator) + ' files: sum of percentages = ' + \
              str(check)
        output.close()

    return concatenated_input_files, text_of_denominator, used_colors


def run_the_imagemaker(args, concatenated_input_files, text_of_denominator, \
                      used_colors):
    """TK"""
    param_scale_bubble = args.bubble_size
    output_dir = args.output
    bubble_type = args.hq_bubbles
    color_mode = args.colors
    text_mode = args.text

    for denominator in sorted(concatenated_input_files):

        for input_filename_long in \
          sorted(concatenated_input_files[denominator]):
            input_filename = \
              concatenated_input_files[denominator][input_filename_long]
            # TK Verified process inputs match
            message = run_visualization(input_filename, \
              input_filename_long, output_dir, denominator, \
              param_scale_bubble, text_of_denominator, bubble_type, \
              used_colors, color_mode, text_mode)
            print str(message)


class NEWICK_tree():


    def __init__(self):
        # TK lines 14-18?
        self.node_id_counter = -1
        self.parent_of_node = Autovivify()
        self.children_of_node = Autovivify()
        self.terminal_children_of_node = Autovivify()
        self.branch_length_of_node = Autovivify()
        self.terminal_nodes = Autovivify()
        self.bootstrap_support_of_node = Autovivify()
        self.cumulative_branch_length_of_node = Autovivify()
        self.y_coord_of_node = Autovivify()
        self.y_coord_of_node_min = Autovivify()
        self.y_coord_of_node_max = Autovivify()
        self.x_coord_of_node = Autovivify()
        self.nodes_of_x_coords = Autovivify()
        self.species_display_order = Autovivify()
        self.x_coord_of_node_connector_circular = Autovivify()
        self.y_coord_of_node_connector_circular = Autovivify()
        self.connecting_points = Autovivify()
        self.x_coord_of_node_min_circular = Autovivify()
        self.y_coord_of_node_min_circular = Autovivify()
        self.x_coord_of_node_max_circular = Autovivify()
        self.y_coord_of_node_max_circular = Autovivify()
        self.x_coord_of_node_circular = Autovivify()
        self.y_coord_of_node_circular = Autovivify()


    def read_tree_topology_from_file(self, filename):
        """This reads the topology of the species tree from a flat-file"""
        input_concatenated = ''
        try:
            input = open(filename, 'r')
        except IOError:
            sys.exit('ERROR: Cannot open inputfile ' + str(filename) + '\n')

        for line in input:
            line = line.strip()
            line = line.replace('\s+', '')
            input_concatenated += line

        input.close()
        # TK Verified input_concatenated
        self.parse_tree_topology_from_string(input_concatenated)


    def parse_tree_topology_from_string(self, tree_string):
        """TK"""
        tree_input_symbols = list(str(tree_string))
        # TK Verified tree_input_symbols
        self.species_tree_parse_node(tree_input_symbols)
        # TK Verified self
        self.connect_terminal_children_to_node(-1)
        # TK Verified self
        self.assign_cumulative_branchlength_to_node(-1,0)
        self.compute_average_cumulative_branch_length()
        # TK Verified self


    def species_tree_parse_node(self, inputarray):
        """
        This is a recursive routine to parse the species-tree topology out
        of an array of input symbols. The input symbols are brackets,
        letters, and numbers. The routine was initially written by Daniel,
        and is quite clever.
        """
        this_node = self.node_id_counter
        self.node_id_counter -= 1
        nr = 0
        children_of_this_node = []

        # TK Verified this_node and self.node_id_counter
        # First symbol must be an open bracket
        next_item = inputarray.pop(0)
        if not next_item == '(':
            sys.exit('ERROR: #1 error parsing species tree\n')

        # Loop over all subnodes of this node
        while True:
            # Next symbol must be either an open bracket or species_name
            try:
                next_item = inputarray.pop(0)
            except IndexError:
                sys.exit('ERROR: #2 error parsing species tree\n')

            # If the next item is an open bracket, call yourself recursively
            if next_item == '(':
                inputarray.insert(0, next_item)
                # TK Verified inputarray
                sub_node = self.species_tree_parse_node(inputarray)
                # TK Verified sub_node
                self.parent_of_node[sub_node] = this_node
                children_of_this_node.append(sub_node)

            # If the next item is a word character, this is
            # a terminal node (a leaf is a species)
            elif re.search('\w+', next_item):
                node_summary = next_item

                next_item = inputarray.pop(0)
                while re.search('[\w\.\:]+', next_item):
                    node_summary += next_item
                    next_item = inputarray.pop(0)

                species = None
                branch_length = None

                # If there are branch-lengths, note them
                if re.search('\A(\w+)\:([\d\.]+\Z)', str(node_summary)):
                    temp = re.search('\A(\w+)\:([\d\.]+\Z)', node_summary)
                    species = temp.group(1)
                    branch_length = temp.group(2)
                else:
                    species = node_summary

                try:
                    species
                except NameError:
                    sys.exit('ERROR: #4a error parsing species tree\n')
                children_of_this_node.append(species)
                self.parent_of_node[species] = this_node
                try:
                    branch_length
                except NameError:
                    pass
                else:
                    self.branch_length_of_node[species] = branch_length
                self.terminal_nodes[species] = 1
                try:
                    next_item
                except NameError:
                    sys.exit('ERROR: #3 error parsing species tree\n')
                inputarray.insert(0,next_item)
                # TK Verified species, branch_length, next_item, inputarray, children_of_this_node
            # If anything else happens, it's illegal
            else:
                sys.exit('ERROR: #4 error parsing species tree\n')

            # One subnode has been parsed. Next symbol must be either a
            # comma (in which case we loop to the next subnode) or a
            # closing bracket (in which case we are done).
            try:
                next_item = inputarray.pop(0)
            except IndexError:
                sys.exit('ERROR: #5 error parsing species tree\n')
            # TK Verified next_item
            if next_item == ',':
                continue
            if next_item == ':':
                continue
            if next_item == ')':
                # Get the branch_length
                branch_summary = ''
                next_item = inputarray.pop(0)

                while inputarray and re.search('[\d\.\:]+', next_item):
                    branch_summary += next_item
                    next_item = inputarray.pop(0)

                temp = re.search('\A(\d*)\:([\d\.]+\Z)', branch_summary)
                bootstrap_support = None
                branch_length = None
                if temp:
                    bootstrap_support = temp.group(1)
                    branch_length = temp.group(2)
                try:
                    branch_length
                except NameError:
                    pass
                else:
                    self.branch_length_of_node[this_node] = branch_length
                try:
                    bootstrap_support
                except NameError:
                    pass
                else:
                    self.bootstrap_support_of_node[this_node] = \
                      bootstrap_support
                inputarray.insert(0, next_item)
                # TK Verified inputarray
                break

            sys.exit('ERROR: #6 error parsing species tree\n')

        for child in children_of_this_node:
            self.children_of_node[this_node][child] = 1

        return this_node


    def connect_terminal_children_to_node(self, node):
        """TK"""
        if re.search('\A[\-\d]+\Z', str(node)) and float(node) < 0:

            for child in self.children_of_node[node]:
                if re.search('\A[\-\d]+\Z', str(child)) and float(child) < 0:
                    self.connect_terminal_children_to_node(child)

                    for subchild in self.terminal_children_of_node[child]:
                        self.terminal_children_of_node[node][subchild] = 1

                else:
                    self.terminal_children_of_node[node][child] = 1
                    self.terminal_children_of_node[child][child] = 1

        else:
            self.terminal_children_of_node[node][node] = 1
            self.terminal_children_of_node[node][node] = 1


    def assign_cumulative_branchlength_to_node(self, node, \
      previous_cumulative_length):
        """TK"""
        branch_length_this_node = 0
        if self.branch_length_of_node[node]:
            branch_length_this_node = self.branch_length_of_node[node]
        new_cumulative_length = float(branch_length_this_node) + \
          float(previous_cumulative_length)
        self.cumulative_branch_length_of_node[node] = new_cumulative_length
        if re.search('\A[\-\d]+\Z', str(node)) and float(node) < 0:

            for child in self.children_of_node[node]:
                self.assign_cumulative_branchlength_to_node(child, \
                  new_cumulative_length)


    def compute_average_cumulative_branch_length(self):
        """TK"""
        valid_node_counter = 0
        branch_length_sum = 0

        for node in self.cumulative_branch_length_of_node:
            # Internal nodes aren't considered here
            if re.search('\A[\-\d]+\Z', str(node)) and float(node) < 0:
                continue
            valid_node_counter += 1
            branch_length_sum += self.cumulative_branch_length_of_node[node]

        self.average_cumulative_branch_length = branch_length_sum / \
          valid_node_counter


    def print_as_verbose_text(self):
        """TK"""

        for node in sorted(self.children_of_node):
            branch_length = 'undefined'
            if self.branch_length_of_node[node]:
                branch_length = self.branch_length_of_node[node]
            children = sorted(self.children_of_node[node])
            nr_children = len(children)
            terminal_children = sorted(self.terminal_children_of_node[node])
            nr_terminal_children = len(terminal_children)
            message = 'node ' + str(node) + ' has branch-length ' + \
              str(branch_length) + ', ' + str(nr_children) + ' children ['
            message += '.'.join(children)
            message += '], and ' + str(nr_terminal_children) + ' terminal' \
              + ' children: ['
            message += '.'.join(terminal_children)
            message += ']\n'
            print message

        for node in sorted(self.parent_of_node):
            if node <= 0:
                continue
            branch_length = 'undefined'
            if self.branch_length_of_node[node]:
                branch_length = self.branch_length_of_node[node]
            print 'node ' + str(node) + ' has branch-length ' + \
              str(branch_length) + '\n'


    def get_species_to_species_tree_distance(self, species1, species2):
        """
        Given two species, the cumulative branch-length needed to traverse
        from one to the other is computed and provided by this routine. The
        routine will return -1 in case of errors.
        """
        try:
            species1
            species2
        except NameError:
            return -1
        if species1 == species2:
            return 0
        if not self.parent_of_node[species1] or \
          not self.parent_of_node[species2]:
            return -1 # TK This is 'next' in Perl...
        parent_traversion_count = Autovivify()
        node1 = species1
        parent_traversion_count[node1] = 1
        try:
            node1
        except NameError:
            pass
        else:

            while node1 is not '-1':
                node1 = self.parent_of_node[node1]
                parent_traversion_count[node1] += 1
                try:
                    node1
                except NameError:
                    break

        node2 = species2
        parent_traversion_count[node2] = 1
        try:
            node2
        except NameError:
            pass
        else:

            while node2 is not '-1':
                node2 = self.parent_of_node[node2]
                parent_traversion_count[node2] += 1
                try:
                    node2
                except NameError:
                    break

        cumulative_tree_distance = 0
        node1 = species1

        while parent_traversion_count[node1] < 2:
            cumulative_tree_distance += self.branch_length_of_node[node1]
            node1 = self.parent_of_node[node1]

        node2 = species2

        while parent_traversion_count[node2] < 2:
            cumulative_tree_distance += self.branch_length_of_node[node2]
            node2 = self.parent_of_node[node2]

        return cumulative_tree_distance


    def optimize_species_display_order(self):
        """
        This makes the tree optically nicer by swapping subtrees such that
        the internal nodes are maximally 'staggered'. This only affects the
        displaying of the tree, not its typology or biological meaning. If
        necessary, manual intervention into the display order is also
        possible in the code below (note: this will fuck up if the
        underlying tree is changed!!)
        """
        self.optimized_display_order_counter = 1
        self.optimize_node_recursively(-1)


    def optimize_node_recursively(self, node):
        """TK"""
        if not (re.search('\A[\-\d]+\Z', str(node)) and float(node) < 0):
            self.species_display_order[node] = \
              self.optimized_display_order_counter
            self.optimized_display_order_counter += 1
            return
        # This node is an internal node, so sort its children
        # by 'terminal_size' and recurse down into them
        children_this_node = self.children_of_node[node]
        nr_terminal_children_per_child = {}

        for child in children_this_node:
            nr_children = len(self.terminal_children_of_node[child])
            nr_terminal_children_per_child[child] = nr_children

        for child in sorted(children_this_node, \
          key = cmp_to_key(optimize_node_recursively_sort_helper, \
          nr_terminal_children_per_child)):
            self.optimize_node_recursively(child)


    def find_last_common_ancestor(self, node1, node2):
        if node1 == node2:
            return node1
        if self.is_in_parental_line(node1, node2):
            return node1
        if self.is_in_parental_line(node2, node1):
            return node2

        while True:
            try:
                self.parent_of_node[node1]
            except NameError:
                return -1
            node1 = self.parent_of_node[node1]
            if self.is_in_parental_line(node1, node2):
                return node1
            if node1 == -1:
                return -1


    def is_in_parental_line(self, putative_parent_node, query_node):
        if putative_parent_node == query_node:
            return 1

        while True:
            query_node = self.parent_of_node[query_node]
            if putative_parent_node == query_node:
                return 1
            if query_node == -1:
                break

        return 0


def run_visualization(input_filename, input_filename_long, output_dir, \
  denominator, param_scale_bubble, text_of_denominator, bubble_type, \
  used_colors, color_mode, text_mode):
    text_of_denominator = text_of_denominator[denominator]
    print str(text_of_denominator)
    mltreemap_results = get_analysis_info(denominator)
    # TK Verified mltreemap_results
    mltreemap_results = get_picture_dimensions(mltreemap_results, \
      text_mode)
    # TK Verified mltreemap_results
    tree = read_tree_topology(mltreemap_results)
    # TK Verified tree
    mltreemap_results = \
      produce_terminal_children_of_strings_of_reference(tree, \
      mltreemap_results, used_colors)
    mltreemap_results = get_support_data(tree, mltreemap_results)
    mltreemap_results = read_RAxML_out(tree, mltreemap_results, \
      input_filename_long)
    # TK Verified mltreemap_results
#    sys.exit(mltreemap_results)
#    sys.exit(str(mltreemap_results['x_coordinate_of_tree_start']))
#    sys.exit(mltreemap_results['name_of_species'])
#    sys.exit(mltreemap_results['group_info'])
#    sys.exit(mltreemap_results['counts_per_species'])
#    sys.exit(mltreemap_results['image'])
#    sys.exit(str(mltreemap_results['x_coordinate_of_tree_end']))
#    sys.exit(mltreemap_results['tax_ids_file'])
#    sys.exit(mltreemap_results['nodes_of_terminal_children_string_reference'])
#    sys.exit(mltreemap_results['counts_per_node'])
#    sys.exit(str(mltreemap_results['x_coordinate_of_label_end']))
#    sys.exit(mltreemap_results['highest_count_per_species'])
#    sys.exit(mltreemap_results['tree_file'])
#    sys.exit(mltreemap_results['descriptions_file'])
#    sys.exit(str(mltreemap_results['species_count']))
    # TK DEBUG draw_the_color_legend--never tested before
    if color_mode > 0:
        draw_the_color_legend(denominator, output_dir, color_mode)
    do_the_drawing(tree, mltreemap_results, param_scale_bubble, \
      input_filename, output_dir, text_of_denominator, bubble_type, \
      used_colors, color_mode, text_mode)
    do_the_drawing_circular(tree, mltreemap_results, \
      param_scale_bubble, input_filename, output_dir, \
      text_of_denominator, bubble_type, used_colors, color_mode, \
      text_mode)
    return 'success'

def get_analysis_info(denominator):
    """TK"""
    mltreemap_results = Autovivify()
    tree_name = ''
    try:
        input = open('tree_data/drawing_info.txt', 'r')
    except IOError:
        sys.exit('ERROR: Can\'t open drawing_info.txt!\n')
    flag = 0

    for line in input:
        if flag < 3:
            flag += 1
            continue
        line = line.strip()
        dn, tn = line.split('\t')
        if dn is denominator:
            tree_name = tn

    input.close()
    if denominator is 'p':
        mltreemap_results['tree_file'] = 'MLTreeMap_reference.tree'
        mltreemap_results['tax_ids_file'] = 'tax_ids_tree_of_life.txt'
        mltreemap_results['descriptions_file'] = \
          'domain_and_color_descriptions.txt'
    elif denominator is 'g':
        mltreemap_results['tree_file'] = 'geba.tree'
        mltreemap_results['tax_ids_file'] = 'tax_ids_geba_tree.txt'
        mltreemap_results['descriptions_file'] = \
          'domain_and_color_descriptions_geba_tree.txt'
    else:
        if not tree_name:
            sys.exit('ERROR: ' + str(denominator) + ' is not ' + \
              'recognized!\n')
        mltreemap_results['tree_file'] = str(tree_name) + '_tree.txt'
        mltreemap_results['tax_ids_file'] = 'tax_ids_' + \
          str(tree_name) + '.txt'
        mltreemap_results['descriptions_file'] = \
          'domain_and_color_descriptions_' + str(tree_name) + '.txt'
    return mltreemap_results


def get_picture_dimensions(mltreemap_results, text_mode):
    """TK"""
    tree_file = mltreemap_results['tree_file']
    try:
        input = open('tree_data/' + str(tree_file), 'r')
    except IOError:
        sys.exit('ERROR: Can\'t open ' + str(tree_file) + '!\n')
    species_count = 0

    for line in input:
        line = line.strip()
        species_count += len(re.findall('\d+:', line))

    input.close()
    # Relations are as follows:
    # 267 species equal 3000 height,
    # while width remains constant at 1000
    image_width = 1000
    tree_height = (species_count * (image_width * 3.0)) / 267
    y_offset = image_width / 20
    image_height = tree_height + y_offset * 2.0
    mltreemap_results['image']['width'] = image_width
    mltreemap_results['image']['tree_height'] = tree_height
    mltreemap_results['image']['image_height'] = image_height
    mltreemap_results['image']['y_offset'] = y_offset
    # Define tree and label margins
    mltreemap_results['x_coordinate_of_label_start'] = int(0.6*image_width)
    mltreemap_results['x_coordinate_of_label_end'] = int(0.9*image_width)
    mltreemap_results['x_coordinate_of_tree_end'] = int(0.7*image_width)
    if not text_mode:
        # If we don't want the whole text to be displayed,
        # we have to reset some values
        mltreemap_results['x_coordinate_of_label_start'] = 0.7 * \
          image_width
        mltreemap_results['x_coordinate_of_label_end'] = 0.73 * \
          image_width
    return mltreemap_results


def produce_terminal_children_of_strings_of_reference(tree, \
  mltreemap_results, used_colors):
    """TK"""
    species_count = 0

    for n in sorted(tree.terminal_children_of_node):

        for color in sorted(used_colors):
            # Set all counts to 0
            mltreemap_results['counts_per_species'][n][color] = 0

        if n > 0:
            species_count += 1
        terminal_children_string = ""

        for term_child in sorted(tree.terminal_children_of_node[n]):
            terminal_children_string += "@" + str(term_child)

        mltreemap_results\
          ['nodes_of_terminal_children_string_reference']\
          [terminal_children_string] = n

    mltreemap_results['species_count'] = species_count
    return mltreemap_results


def read_RAxML_out(tree, mltreemap_results, input_filename_long):
    """TK"""
    nodes = Autovivify()
    input_filename = input_filename_long
    print 'input file: ' + str(input_filename)
    try:
        input = open(input_filename, 'r')
    except IOError:
        sys.exit('ERROR: Can\'t open ' + str(input_filename) + '\n')
    placements = Autovivify()
    terminal_children_of_placements = Autovivify()
    mltreemap_results['highest_count_per_species'] = 0
    placement = 0

    for line in input:
        line = line.strip()
        placement += 1
        color = 0
        if re.search('colorcode_(.+)\Z', line):
            color = re.search('colorcode_(.+)\Z', line).group(1)
        if color == 0:
            sys.exit('ERROR: Color info could not be read!!!!\n')
        try:
            color
        except NameError:
            sys.exit('ERROR: Color info could not be read!!!!\n')
        temp = re.findall('\((\d+)\)', line)
        for _ in temp:
            terminal_children_of_placements[placement][_] = 1

        terminal_children_string_of_placement = ''

        for terminal_child_of_placement in \
          sorted(terminal_children_of_placements[placement]):
            terminal_children_string_of_placement += '@' + \
              terminal_child_of_placement

        try:
            mltreemap_results\
              ['nodes_of_terminal_children_string_reference']\
              [terminal_children_string_of_placement]
        except NameError:
            pass
            # Attention! This circumvents a problem with the unrooted
            # RAxML tree!
            # continue
            # sys.exit('ERROR: A subtree\n(' + \
            # str(terminal_children_string_of_placement) + ')\n, as' + \
            # ' written in the parsed RAxML file, does not exist in' + \
            # ' the reference tree...\n')
        else:
            node = mltreemap_results\
              ['nodes_of_terminal_children_string_reference']\
              [terminal_children_string_of_placement]
            _pw = re.search('Placement weight (.+)%', line)
            if _pw:
                bootstrap = _pw.group(1)
                nodes[node][color] = bootstrap
                node_weight = float(bootstrap) / 100
                mltreemap_results = \
                  distribute_node_weight_by_topology(tree, node, \
                  node_weight, mltreemap_results, color)

    input.close()

    for node in sorted(nodes):
        is_first_color = 1

        for color in sorted(nodes[node]):
            mltreemap_results['counts_per_node'][node]['colors']\
              [color] = float(nodes[node][color]) / 100
            if is_first_color == 1:
                mltreemap_results['counts_per_node'][node]\
                  ['total_weight'] = 0
            is_first_color = 0
            mltreemap_results['counts_per_node'][node]['total_weight'] \
              += float(nodes[node][color]) / 100

    return mltreemap_results


def distribute_node_weight_by_topology(tree, node, node_weight, \
  mltreemap_results, color):
    """TK"""
    terminal_children = tree.terminal_children_of_node[node]
    nr_of_children = len(terminal_children)
    fractional_weight = node_weight / nr_of_children

    for child in terminal_children:
        mltreemap_results['counts_per_species'][child][color] += \
          fractional_weight
        total_weight_of_child = 0

        for color in \
          sorted(mltreemap_results['counts_per_species'][child]):
            color_weight = \
              mltreemap_results['counts_per_species'][child][color]
            total_weight_of_child += color_weight

        if total_weight_of_child > \
          mltreemap_results['highest_count_per_species']:
            mltreemap_results['highest_count_per_species'] = \
              total_weight_of_child

    return mltreemap_results


def get_support_data(tree, mltreemap_results):
    """TK"""
    # Get species names
    try:
        input = open('tree_data/' + \
          str(mltreemap_results['tax_ids_file']), 'r')
    except IOError:
        sys.exit('ERROR: Cannot read \'tree_data/' + \
          str(mltreemap_results['tax_ids_file']) + '\'!\n')

    for line in input:
        line = line.strip()
        if re.search('\A\#', line):
            continue
        species, _rest = line.split(None, 1)
        rest = _rest.split(None)
        name = ' '.join(rest)
        mltreemap_results['name_of_species'][species] = name

    input.close()
    try:
        input = open('tree_data/' + \
          str(mltreemap_results['descriptions_file']), 'r')
    except IOError:
        sys.exit('ERROR: Can\'t read tree_data/' + \
          str(mltreemap_results['descriptions_file']) + '\n')

    for line in input:
        line = line.strip()
        if re.search('\A\#', line) or line is '':
            continue
        start_taxon, end_taxon, background_red, background_green, \
          background_blue, subtree_name = line.split()
        if not start_taxon:
            continue
        group_color = 'rgb(' + str(background_red) + ',' + \
          str(background_green) + ',' + str(background_blue) + ')'
        y_coord_of_node_min = tree.y_coord_of_node_min[start_taxon]
        mltreemap_results['group_info'][y_coord_of_node_min]\
          [start_taxon][end_taxon]['color'] = str(group_color)
        mltreemap_results['group_info'][y_coord_of_node_min]\
          [start_taxon][end_taxon]['name'] = str(subtree_name)

    input.close()
    return mltreemap_results


def draw_the_color_legend(self, denominator, output_dir, color_mode):
    """TK"""
    try:
        input = open(str(output_dir) + str(denominator) + \
          '_color_legend.txt', 'r')
    except IOError:
        sys.exit('ERROR: Can\'t open ' + str(output_dir) + \
          str(denominator) + '_color_legend.txt\n')
    nr_of_datasets = 0
    color_info = Autovivify()
    namelength_max = 0

    for line in input:
        line = line.strip()
        color, filename = line.split('\t')
        try:
            color
            filename
        except NameError:
            sys.exit('ERROR: Something is wrong with ' + \
              'color_legend.txt!\n')
        color_info[color] = filename
        if len(filename) > namelength_max:
            namelength_max = len(filename)
        nr_of_datasets += 1

    input.close()
    output_file_name = str(output_dir) + str(denominator) + \
      '_color_legend.svg'
    try:
        input = open(output_file_name, 'w')
    except IOError:
        sys.exit('ERROR: Can\'t create ' + output_file_name + '!\n')
    legend_width = 7 * namelength_max
    y_offset = 12
    legend_height = y_offset + y_offset * nr_of_datasets
    legend = svgwrite.Drawing(filename=str(output_dir) + \
      str(denominator) + '_color_legend_.svg', size = \
      (str(legend_width) + 'px',str(legend_height) + 'px'))
    y_pos = y_offset

    for color in sorted(color_info):
        filename = color_info[color]
        legend.add(legend.rect(insert = (10, y_pos - y_offset), \
          size = (str(y_offset), str(y_offset)), rx = None, ry = None, \
          fill = str(color)))
        legend.add(legend.text(str(filename), (y_offset + 10, y_pos), \
          style = 'font-family:Verdana; font-size:10'))
        y_pos += y_offset

    legend.save()
    output.write(legend.tostring())
    output.close()


def do_the_drawing(tree, mltreemap_results, param_scale_bubble, \
  input_file_name, output_dir, text_of_denominator, bubble_type, \
  used_colors, color_mode, text_mode):
    """TK"""
    image_width = mltreemap_results['image']['width']
    tree_height = mltreemap_results['image']['tree_height']
    image_height = mltreemap_results['image']['image_height']
    y_offset = mltreemap_results['image']['y_offset']
    species_count = mltreemap_results['species_count']
    print 'Creating the picture...'
    output_file_name = str(output_dir) + str(input_file_name) + \
      '_image_straight.svg'
    image = svgwrite.Drawing(filename = str(output_dir) + \
      str(input_file_name) + '_image_straight.svg', \
      size = (str(image_width) + 'px', str(image_height) + 'px'))
    image.add(image.text(str(text_of_denominator), \
      (int(image_width * 0.01), y_offset / 2), style = \
      'font-family:Verdana; font-size:' + str(image_width / 50)))
    placements = ''
    image = draw_group_colors(image, tree, mltreemap_results)
    # TK Verified Visually, looks fine here
    image = draw_edges(image, tree, mltreemap_results)
    # TK Verified Visually, looks fine here
    image = draw_guide_lines_and_leaf_names(image, tree, \
      mltreemap_results, text_mode)
    # TK Verified Visually, looks fine here
    image, placements = draw_placement_bubbles(image, tree, \
      mltreemap_results, param_scale_bubble, bubble_type, 0)
    # TK Verified Visually, looks fine here
    image = draw_percents_and_placement_bars(image, placements, tree, \
      mltreemap_results, used_colors, color_mode, text_mode)
    # TK Verified Visually, looks fine here
    image.save()


def do_the_drawing_circular(tree, mltreemap_results, \
  param_scale_bubble, input_file_name, output_dir, \
  text_of_denominator, bubble_type, used_colors, color_mode, text_mode):
    """TK"""
    image_width = mltreemap_results['image']['width']
    tree_height = mltreemap_results['image']['tree_height']
    image_height = mltreemap_results['image']['image_height']
    y_offset = mltreemap_results['image']['y_offset']
    species_count = mltreemap_results['species_count']
    # For large trees, the image_diameter_circular needs to be adjusted
    image_diameter_circular = 2 * image_width
    mltreemap_results['image_circular']['diameter'] = \
      image_diameter_circular
    print 'Creating the picture...'
    output_file_name = str(output_dir) + str(input_file_name) + \
      '_image_circular.svg'
    image = svgwrite.Drawing(filename = str(output_dir) + \
      str(input_file_name) + '_image_circular.svg', size = \
      (str(image_diameter_circular) + 'px', \
      str(image_diameter_circular) + 'px'))
    image.add(image.text(str(text_of_denominator), \
      (image_width * 0.01, y_offset / 2), style = \
      'font-family:Verdana; font-size:' + str(image_width / 50)))
    placements = ''
    tree, mltreemap_results = prepare_coordinates_circular(tree, \
      mltreemap_results)
    image = draw_group_colors_circular(image, tree, mltreemap_results)
    # TK Verified Viusally, looks fine here
    image = draw_edges_circular(image, tree, mltreemap_results)
    # TK Verified Visually, looks fine here
    image = draw_guide_lines_and_leaf_names_circular(image, tree, \
      mltreemap_results, bubble_type, text_mode)
    # TK Verified Visually, looks fine here
    image, placements = draw_placement_bubbles(image, tree, \
      mltreemap_results, param_scale_bubble, bubble_type, 1)
    # TK Verified Visually, looks fine here
    image = draw_percents_and_placement_bars_circular(image, \
      placements, tree, mltreemap_results, bubble_type, used_colors, \
      color_mode, text_mode)
    # TK Verified Visually, looks fine here
    image.save()


def prepare_coordinates_circular(tree, mltreemap_results):

    for node in sorted(tree.y_coord_of_node):
        y_coord_of_node = tree.y_coord_of_node[node]
        x_coord_of_node = tree.x_coord_of_node[node]
        branch_length = tree.branch_length_of_node[node]
        x_coord_of_parent_node = x_coord_of_node - branch_length
        # Transform the coordinates of the node itself
        px, py = calculate_coordinates_circular(mltreemap_results, \
          x_coord_of_node, y_coord_of_node)
        tree.x_coord_of_node_circular[node] = px
        tree.y_coord_of_node_circular[node] = py
        # Transform the coordinates of the connecting point
        # of the node to the parent
        px, py = calculate_coordinates_circular(mltreemap_results, \
          x_coord_of_parent_node, y_coord_of_node)
        tree.x_coord_of_node_connector_circular[node] = px
        tree.y_coord_of_node_connector_circular[node] = py
        # Save the connecting point information with parent information
        if node != -1:
            parent_of_node = tree.parent_of_node[node]
            tree.connecting_points[parent_of_node][node]\
              ['x_coordinate_of_connector'] = px
            tree.connecting_points[parent_of_node][node]\
              ['y_coordinate_of_connector'] = py
        # Done
        if node > 0:
            y_coord_of_node_min = tree.y_coord_of_node_min[node]
            y_coord_of_node_max = tree.y_coord_of_node_max[node]
            px, py = calculate_coordinates_circular(mltreemap_results, \
              x_coord_of_node, y_coord_of_node_min)
            tree.x_coord_of_node_min_circular[node] = px
            tree.y_coord_of_node_min_circular[node] = py
            px, py = calculate_coordinates_circular(mltreemap_results, \
              x_coord_of_node, y_coord_of_node_max)
            tree.x_coord_of_node_max_circular[node] = px
            tree.y_coord_of_node_max_circular[node] = py

    return tree, mltreemap_results


def calculate_coordinates_circular(mltreemap_results, \
  x_coord_of_point, y_coord_of_point):
    """
    TK
    Note: Some re-programming could merge this process with
    the one that draws the pie-chart bubbles
    """
    image_diameter_circular = \
      mltreemap_results['image_circular']['diameter']
    tree_height = mltreemap_results['image']['tree_height']
    y_offset = mltreemap_results['image']['y_offset']
    pi = math.pi

    # la
    # ----------P(px|py)
    # | /
    # |alpha /
    # | /
    #lb| /x-pos
    # |/
    # M

    # Proportion: tree_height equals 2 pi * 0.95
    alpha = ((y_coord_of_point - y_offset) * 2 * pi * 0.95) / \
      tree_height
    center_x = image_diameter_circular / 2
    center_y = image_diameter_circular / 2
    lb = math.cos(alpha) * x_coord_of_point
    la = math.sin(alpha) * x_coord_of_point
    px = center_x + la
    py = center_y - lb
    return px, py


def draw_group_colors_circular(image, tree, mltreemap_results):
    image_width = mltreemap_results['image']['width']
    tree_height = mltreemap_results['image']['tree_height']
    image_height = mltreemap_results['image']['image_height']
    y_offset = mltreemap_results['image']['y_offset']
    y_scaling_factor = tree_height # Note: This works because
                                   # y values range from 0 to 1.
    x_coordinate_of_label_start = mltreemap_results\
      ['x_coordinate_of_label_start']
    x_coordinate_of_label_end = mltreemap_results\
      ['x_coordinate_of_label_end']
    groups = image.add(image.g(id='groups'))
    fontsize = image_width / 60.0
    pi = math.pi

    # The group color band looks somewhat as follows:
    #
    # C-----------D
    # \ /
    # \ /
    # \ /
    # A---B
    #
    # A(xa,ya), B(xb,yb), C(xc,yc), D(xd,yd)

    is_first_label = 1
    text_y = y_offset / 2
    x_pos_of_text = image_width * 2 - (image_width / 4)

    for y_coord_of_node_min in sorted(mltreemap_results['group_info']):

        for start_taxon in sorted(mltreemap_results['group_info']\
          [y_coord_of_node_min]):

            for end_taxon in sorted(mltreemap_results['group_info']\
              [y_coord_of_node_min][start_taxon]):
                color = mltreemap_results['group_info']\
                  [y_coord_of_node_min][start_taxon][end_taxon]['color']
                name = mltreemap_results['group_info']\
                  [y_coord_of_node_min][start_taxon][end_taxon]['name']
                if name is '#':
                    name = ''
                ya_linear = y_coord_of_node_min # == yc_linear
                yb_linear = tree.y_coord_of_node_max[end_taxon]
                  # == yd_linear
                xa_linear = x_coordinate_of_label_start # == xb_linear
                xc_linear = x_coordinate_of_label_end # == xd_linear
                draw_trapezoid(mltreemap_results, groups, image, color, \
                  xa_linear, xc_linear, ya_linear, yb_linear)
                # Prepare and write the group labels
                if not name:
                    continue
                if is_first_label:
                    text = 'Group names (clockwise):'
                    groups.add(image.text(str(text), (x_pos_of_text, \
                      text_y), style = \
                      'font-family:Verdana; font-size:' + \
                      str(fontsize) + '; fill: rgb(0,0,0)'))
                    is_first_label = 0
                name = name.replace('_', ' ')
                if re.search('rgb\((.+),(.+),(.+)\)', color):
                    temp = re.search('rgb\((.+),(.+),(.+)\)', color)
                    red = int(temp.group(1)) - 80
                    green = int(temp.group(2)) - 80
                    blue = int(temp.group(3)) - 80
                    if red < 0:
                        red = 0
                    if green < 0:
                        green = 0
                    if blue < 0:
                        blue = 0
                    color = 'rgb(' + str(red) + ',' + str(green) + ',' \
                      + str(blue) + ')' # Make the colors darker...
                                        # Otherwise, the writing is
                                        # almost invisible
                else:
                    sys.exit('ERRROR: Parsing error with ' + \
                      str(color) + '\n')
                text_y += fontsize * 1.2
                groups.add(image.text(str(name), (x_pos_of_text, \
                  int(text_y)), style = 'font-family: Verdana; ' + \
                  'font-size: ' + str(fontsize) + '; fill:' + \
                  str(color)))

    return image


def draw_edges_circular(image, tree, mltreemap_results):
    """TK"""
    image_width = mltreemap_results['image']['width']
    tree_height = mltreemap_results['image']['tree_height']
    y_offset = mltreemap_results['image']['y_offset']
    image_diameter_circular = mltreemap_results['image_circular']\
      ['diameter']
    allready_drawn_connector = Autovivify()
    edges = image.add(image.g(id = 'edges'))
    stroke_width = image_width / 1000

    for node in sorted(tree.y_coord_of_node):
        edge_color = 'rgb(100,100,100)'
        if node != -1:
            # First get all the necessary information
            x_coord_of_node = tree.x_coord_of_node_circular[node]
            y_coord_of_node = tree.y_coord_of_node_circular[node]
            x_coord_of_node_connector = \
              tree.x_coord_of_node_connector_circular[node]
            y_coord_of_node_connector = \
              tree.y_coord_of_node_connector_circular[node]
            # Second, we draw the line from the node to
            # the connection position of its parent
            edges.add(image.line((x_coord_of_node, y_coord_of_node), \
              (x_coord_of_node_connector, y_coord_of_node_connector), \
              stroke_linecap = 'round', stroke = str(edge_color), \
              stroke_width = str(stroke_width), stroke_opacity = 1))
        # Third, we draw the line between the connecting
        # points (if not already done)
        if allready_drawn_connector[node]:
            continue
        allready_drawn_connector[node] = 1
        sweep_flag = 1
        large_arc_flag = 0
        x1_coordinate_of_connector = ''
        y1_coordinate_of_connector = ''
        x2_coordinate_of_connector = ''
        y2_coordinate_of_connector = ''
        count = 0

        for child in sorted(tree.connecting_points[node], reverse=True):
            if count == 0:
                x1_coordinate_of_connector = tree.connecting_points\
                  [node][child]['x_coordinate_of_connector']
                y1_coordinate_of_connector = tree.connecting_points\
                  [node][child]['y_coordinate_of_connector']
            else:
                x2_coordinate_of_connector = tree.connecting_points\
                  [node][child]['x_coordinate_of_connector']
                y2_coordinate_of_connector = tree.connecting_points\
                  [node][child]['y_coordinate_of_connector']
            count += 1

        # The points (x1/y1) (x2/y2) have to be sorted otherwise
        # the arc function won't work properly
        if not x1_coordinate_of_connector:
            continue
        if y1_coordinate_of_connector > y2_coordinate_of_connector:
            temp_coordinate = y1_coordinate_of_connector
            y1_coordinate_of_connector = y2_coordinate_of_connector
            y2_coordinate_of_connector = temp_coordinate
            temp_coordinate = x1_coordinate_of_connector
            x1_coordinate_of_connector = x2_coordinate_of_connector
            x2_coordinate_of_connector = temp_coordinate
        if x1_coordinate_of_connector < (image_diameter_circular / 2):
            sweep_flag = 0
        # Sorting done
        radius_of_node = tree.x_coord_of_node[node]
        edges.add(image.path(d = 'M ' + \
          str(x1_coordinate_of_connector) + ',' + \
          str(y1_coordinate_of_connector) + ' A ' + \
          str(radius_of_node) + ',' + str(radius_of_node) + ' 0 ' + \
          str(large_arc_flag) + ',' + str(sweep_flag) + ' ' + \
          str(x2_coordinate_of_connector) + ',' + \
          str(y2_coordinate_of_connector), opacity = 1, stroke_width = \
          str(stroke_width), stroke = str(edge_color), fill = 'none'))

    return image


def draw_guide_lines_and_leaf_names_circular(image, tree, \
  mltreemap_results, bubble_type, text_mode):
    """TK"""
    image_width = mltreemap_results['image']['width']
    tree_height = mltreemap_results['image']['tree_height']
    y_offset = mltreemap_results['image']['y_offset']
    species_count = mltreemap_results['species_count']
    x_coordinate_of_label_start = mltreemap_results\
      ['x_coordinate_of_label_start']
    guide_lines_and_labels = image.add(image.g(id = \
      'guide_lines_and_labels'))
    edge_color = 'rgb(220,220,220)'
    stroke_width = image_width / 1000
    fontsize = image_width / 100
    # The picture has been optimized for 267 species.
    # If we have more, downsize the font.
    if species_count > 267:
        fontsize *= 267 / species_count
    x_gap = image_width / 100

    for node in sorted(tree.y_coord_of_node):
        if node <= 0:
            continue
        x_coord_of_node = tree.x_coord_of_node[node]
        y_coord_of_node = tree.y_coord_of_node[node]
        x_coord_of_text = x_coordinate_of_label_start + x_gap
        y_coord_of_text = y_coord_of_node
        max_text_length = 37
        if x_coord_of_node + x_gap < x_coordinate_of_label_start:
            x1_pos_linear = x_coord_of_node + x_gap / 2
            x2_pos_linear = x_coordinate_of_label_start - x_gap / 2
            x1, y1 = calculate_coordinates_circular(mltreemap_results, \
              x1_pos_linear, y_coord_of_node)
            x2, y2 = calculate_coordinates_circular(mltreemap_results, \
              x2_pos_linear, y_coord_of_node)
            if text_mode:
                guide_lines_and_labels.add(image.line((x1, y1), \
                  (x2, y2), stroke_linecap = 'round', stroke = \
                  str(edge_color), stroke_width = str(stroke_width), \
                  stroke_opacity = 1))
        else:
            if x_coord_of_text < x_coord_of_node + x_gap:
                x_coord_of_text = x_coord_of_node + x_gap
            max_text_length = 32
        node_name = mltreemap_results['name_of_species'][node]
        temp = re.compile('\A(.{%s}).'%str(max_text_length))
        if temp.search(str(node_name)):
            node_name = str(temp.search(str(node_name)).group(1)) + '...'
        rot_angle = (((y_coord_of_text - y_offset) * 360 * 0.95) / \
          tree_height) - 90 # Proportion: tree_height == 360 * 0.95
        x_text, y_text = calculate_coordinates_circular(\
          mltreemap_results, x_coord_of_text, y_coord_of_text)
        try:
            node_name
        except NameError:
            sys.exit('ERROR: ' + str(node) + ' has no name!\n')
        if rot_angle + 90 <= 180:
            if text_mode:
                guide_lines_and_labels.add(image.text(str(node_name), \
                  (x_text, y_text), transform = 'rotate(' + \
                  str(rot_angle) + ' ' + str(x_text) + ',' + \
                  str(y_text) + ')', style = \
                  'font-family:Verdana; font-size:' + str(fontsize)))
        else:
            rot_angle += 180
            if text_mode:
                # TK line 775 best translation guess below. DEBUG?
                guide_lines_and_labels.add(image.text(str(node_name), \
                  (x_text, y_text), transform = 'rotate(' + \
                  str(rot_angle) + ' ' + str(x_text) + ',' + \
                  str(y_text) + ')', style = \
                  'font-family:Verdana; font-size:' + str(fontsize) + \
                  '; text-anchor:end'))

    return image


def draw_percents_and_placement_bars_circular(image, placements, tree, \
  mltreemap_results, bubble_type, used_colors, color_mode, text_mode):
    """TK"""
    image_width = mltreemap_results['image']['width']
    tree_height = mltreemap_results['image']['tree_height']
    species_count = mltreemap_results['species_count']

    # The placement bar looks somewhat as follows:
    # 
    #  C-----------D
    #   \         /
    #    \       /
    #     \     /
    #      A---B
    # 
    # A(xa,ya), B(xb,yb), C(xc,yc), D(xd,yd)

    for node in sorted(tree.y_coord_of_node):
        if node <= 0:
            continue

        # First, get all necessary information
        y_coord_of_node = tree.y_coord_of_node[node]
        x_coordinate_of_label_end = mltreemap_results\
          ['x_coordinate_of_label_end']
        y_offset = mltreemap_results['image']['y_offset']
        highest_fraction_raw = mltreemap_results\
          ['highest_count_per_species'] * 100
        fraction_raw = 0

        # Prepare the text
        text = ''
        all_fractions_0 = 1
        fraction_total = 0

        for color in sorted(used_colors):
            # text as generated here is only used if color_mode == 1
            # (i.e. show different datasets in different colors and 
            # percentages for all of them).
            if mltreemap_results['counts_per_species'][node][color]:
                fraction_raw = mltreemap_results['counts_per_species']\
                  [node][color] * 100
            fraction = int(fraction_raw * 100 + 0.5) / 100.0
            fraction_total += fraction
            if fraction > 0:
                all_fractions_0 = 0
            text += '(' + str('{0:.2f}'.format(fraction)) + '%)'

        fraction_total = '%.2f' % fraction_total
        if color_mode == 11:
            text = '(Total: ' + str(fraction_total) + '%)'
        if text_mode:
            fontsize = image_width / 100.0
            # The picture has been optimized for 267 species.
            # If we have more, downsize the font
            if species_count > 267:
                fontsize *= 267 / species_count
            y_coord_of_text = y_coord_of_node
            x_gap = image_width / 200.0
            x_coord_of_text = x_coordinate_of_label_end - x_gap
            if all_fractions_0 != 0:
                continue
            rot_angle = (((y_coord_of_text - y_offset) * 360 * 0.95) / \
              tree_height) - 90 # Proportion: tree_height == 360 * 0.95
            x_text, y_text = calculate_coordinates_circular(\
              mltreemap_results, x_coord_of_text, y_coord_of_text)
            if rot_angle + 90 <= 180:
                placements.add(image.text(str(text), (x_text, \
                  y_text), transform = 'rotate(' + str(rot_angle) + \
                  ' ' + str(x_text) + ',' + str(y_text) + ')', \
                  style = 'font-family:Verdana; font-size:' + \
                  str(fontsize) + '; text-anchor:end'))
            else:
                rot_angle += 180
                # TK line 855 this is my best guess; DEBUG
                placements.add(image.text(str(text), (x_text, \
                  y_text), transform = 'rotate(' + str(rot_angle) + \
                  ' ' + str(x_text) + ',' + str(y_text) + ')', \
                  style = 'font-family:Verdana; font-size:' + \
                  str(fontsize)))
        only_one_color = 1
        # Prepare the placement bars
        x_offset = 0

        for color in sorted(used_colors):
            if mltreemap_results['counts_per_species'][node][color]:
                fraction_raw = mltreemap_results['counts_per_species']\
                  [node][color] * 100
            if fraction_raw <= 0:
                continue
            y_min = tree.y_coord_of_node_min[node]
            y_max = tree.y_coord_of_node_max[node]
            height_tot = y_max - y_min
            height = height_tot * 0.9
            # ya_linear == yc_linear
            ya_linear = y_min + ((height_tot - height) / 2)
            # yb_linear == yd_linear
            yb_linear = ya_linear + height
            start_y = y_min + ((height_tot - height) / 2)
            x_gap2 = image_width / 200
            start_x = x_coordinate_of_label_end + x_gap2
            max_length = (image_width - image_width / 100) - start_x
            start_x += x_offset
            # The highest fraction is max_length
            fractional_length = max_length * (fraction_raw / \
              highest_fraction_raw)
            x_offset += fractional_length
            xa_linear = start_x # == xb_linear
            xc_linear = start_x + fractional_length # == xd_linear
            draw_trapezoid(mltreemap_results, placements, image, color, \
              xa_linear, xc_linear, ya_linear, yb_linear)

    return image


def draw_trapezoid(mltreemap_results, svg, image, color, xa_linear, \
  xc_linear, ya_linear, yb_linear):
    """TK"""
    # The trapezoid looks somewhat as follows:
    #
    # C-----------D
    #  \         /
    #   \       /
    #    \     /
    #     A---B
    #
    # A(xa,ya), B(xb,yb), C(xc,yc), D(xd, yd)

    tree_height = mltreemap_results['image']['tree_height']
    x_coordinate_of_label_start = mltreemap_results\
      ['x_coordinate_of_label_start']
    x_coordinate_of_label_end = mltreemap_results\
      ['x_coordinate_of_label_end']
    large_arc_flag = 0
    sweep_flag = 1
    sweep_flag2 = 0
    if (yb_linear - ya_linear) >= (tree_height / 0.95) / 2:
        # ie. The group color band will cover more than 180*
        large_arc_flag = 1
    xa, ya = calculate_coordinates_circular(mltreemap_results,\
      xa_linear, ya_linear)
    xb, yb = calculate_coordinates_circular(mltreemap_results,\
      xa_linear, yb_linear)
    xc, yc = calculate_coordinates_circular(mltreemap_results,\
      xc_linear, ya_linear)
    xd, yd = calculate_coordinates_circular(mltreemap_results,\
      xc_linear, yb_linear)
    radius_AB = x_coordinate_of_label_start
    radius_CD = x_coordinate_of_label_end
    # TK Verified xb yb radius_CD
    svg.add(image.path(d = 'M ' + str(xa) + ',' + str(ya) + ' A ' + \
      str(radius_AB) + ',' + str(radius_AB) + ' 0 ' + \
      str(large_arc_flag) + ', ' + str(sweep_flag) + ' ' + str(xb) + \
      ',' + str(yb) + ' L ' + str(xd) + ',' + str(yd) + ' A ' + \
      str(radius_CD) + ',' + str(radius_CD) + ' 0 ' + \
      str(large_arc_flag) + ', ' + str(sweep_flag2) + ', ' + str(xc) + \
      ',' + str(yc) + ' z', opacity = 1, fill = str(color)))


def draw_group_colors(image, tree, mltreemap_results):
    image_width = mltreemap_results['image']['width']
    tree_height = mltreemap_results['image']['tree_height']
    y_scaling_factor = tree_height # Note: This works because y values 
                                   #       range from 0-1
    x_coordinate_of_label_start = mltreemap_results\
      ['x_coordinate_of_label_start']
    x_coordinate_of_label_end = mltreemap_results\
      ['x_coordinate_of_label_end']
    groups = image.add(image.g(id = 'groups'))
    rx = image_width / 400
    ry = rx
    fontsize = image_width / 125

    for y_coord_of_node_min in sorted(mltreemap_results['group_info']):

        for start_taxon in sorted(mltreemap_results['group_info']\
          [y_coord_of_node_min]):

            for end_taxon in sorted(mltreemap_results['group_info']\
              [y_coord_of_node_min][start_taxon]):
                color = mltreemap_results['group_info']\
                  [y_coord_of_node_min][start_taxon][end_taxon]['color']
                name = mltreemap_results['group_info']\
                  [y_coord_of_node_min][start_taxon][end_taxon]['name']
                if name is '#':
                    name = ''
                start_y = y_coord_of_node_min
                end_y = tree.y_coord_of_node_max[end_taxon]
                width = x_coordinate_of_label_end - \
                  x_coordinate_of_label_start
                height = end_y - start_y
                groups.add(image.rect(insert = \
                  (x_coordinate_of_label_start, start_y), size = \
                  (str(width), str(height)), rx = rx, ry = ry, fill = \
                  str(color)))
                # Prepare and write the group labels
                x_pos_of_text = image_width - (image_width / 100)
                name = name.replace('_', ' ')
                temp = re.search('rgb\((.+),(.+),(.+)\)', color)
                if temp:
                    red = int(temp.group(1)) - 80
                    green = int(temp.group(2)) - 80
                    blue = int(temp.group(3)) - 80
                    if red < 0:
                        red = 0
                    if green < 0:
                        green = 0
                    if blue < 0:
                        blue = 0
                    color = 'rgb(' + str(red) + ',' + str(green) + ',' \
                      + str(blue) + ')'
                else:
                    sys.exit('ERROR: Parsing error with ' + str(color) \
                      + '\n')
                groups.add(image.text(str(name), (x_pos_of_text, \
                  start_y), transform = 'rotate(90 ' + \
                  str(x_pos_of_text) + ',' + str(start_y) + ')', \
                  style = 'font-family:Verdana; font-size:' + \
                  str(fontsize) + '; fill:' + str(color)))

    return image


def draw_edges(image, tree, mltreemap_results):
    """TK"""
    image_width = mltreemap_results['image']['width']
    tree_height = mltreemap_results['image']['tree_height']
    y_offset = mltreemap_results['image']['y_offset']
    allready_drawn_parents = Autovivify()
    edges = image.add(image.g(id = 'edges'))
    stroke_width = image_width / 1000

    # The tree looks schematically as follows and 
    # thus has horizontal and vertical lines:
    #     |-----
    # ----|
    #     |-------

    for node in sorted(tree.y_coord_of_node):
        if node == -1:
            continue
        edge_color = 'rgb(100,100,100)'
        # First, get all the necessary information
        parent_of_node = tree.parent_of_node[node]
        branch_length = tree.branch_length_of_node[node]
        x_coord_of_node = tree.x_coord_of_node[node]
        y_coord_of_node = tree.y_coord_of_node[node]
        x_coord_of_parent_node = tree.x_coord_of_node[parent_of_node]
        y_coord_of_parent_node_center = tree.y_coord_of_node\
          [parent_of_node]
        y_distance_parent_to_node = y_coord_of_node - \
          y_coord_of_parent_node_center

        # Second, we draw the horizontal lines from
        # the node to the x-position of its parent.
        edges.add(image.line(start = (x_coord_of_node, \
          y_coord_of_node), end = (x_coord_of_parent_node, \
          y_coord_of_node), stroke_linecap = 'round', stroke = \
          str(edge_color), stroke_width = str(stroke_width), \
          stroke_opacity = 1))

        # Third, we draw the vertical line (if not already done)
        if allready_drawn_parents[parent_of_node]:
            continue
        allready_drawn_parents[parent_of_node] = 1
        edges.add(image.line(start = (x_coord_of_parent_node, \
          y_coord_of_node), end = (x_coord_of_parent_node, \
          y_coord_of_node - 2 * y_distance_parent_to_node), \
          stroke_linecap = 'round', stroke = str(edge_color), \
          stroke_width = str(stroke_width), stroke_opacity = 1))

    return image


def draw_guide_lines_and_leaf_names(image, tree, mltreemap_results, \
  text_mode):
    """TK"""
    image_width = mltreemap_results['image']['width']
    tree_height = mltreemap_results['image']['tree_height']
    x_coordinate_of_label_start = mltreemap_results\
      ['x_coordinate_of_label_start']
    guide_lines_and_labels = image.add(image.g(\
      id = 'guide_lines_and_labels'))
    edge_color = 'rgb(220,220,220)'
    stroke_width = image_width / 1000
    fontsize = image_width / 125
    y_offset2 = 0.3 * fontsize
    x_gap = image_width / 100

    for node in sorted(tree.y_coord_of_node):
        if node <= 0:
            continue
        x_coord_of_node = tree.x_coord_of_node[node]
        y_coord_of_node = tree.y_coord_of_node[node]
        x_coord_of_text = x_coordinate_of_label_start + x_gap
        max_text_length = 47
        if x_coord_of_node + x_gap < x_coordinate_of_label_start:
            if text_mode:
                guide_lines_and_labels.add(image.line(start = \
                  (x_coord_of_node + x_gap / 2, y_coord_of_node), \
                  end = (x_coordinate_of_label_start - x_gap / 2, \
                  y_coord_of_node), stroke_linecap = 'round', stroke = \
                  str(edge_color), stroke_width = str(stroke_width), \
                  stroke_opacity = 1))
        else:
            if x_coord_of_text < x_coord_of_node + x_gap:
                x_coord_of_text = x_coord_of_node + x_gap
                max_text_length = 42
        node_name = mltreemap_results['name_of_species'][node]
        temp = re.compile('\A(.{%s}).' % max_text_length)
        if temp.search(node_name):
            node_name = str(temp.search(node_name).group(1)) + '...'
        try:
            node_name
        except NameError:
            sys.exit('ERROR: ' + str(node) + ' has no name!\n')
        if text_mode:
            guide_lines_and_labels.add(image.text(str(node_name), \
              (x_coord_of_text, y_coord_of_node + y_offset2), style = \
              'font-family:Verdana; font-size:' + str(fontsize)))

    return image


def draw_placement_bubbles(image, tree, mltreemap_results, \
  param_scale_bubble, bubble_type, is_circular_image):
    """
    Note: This subroutine could easily be merged into 'draw_edges'. The 
    only reason why it exists is that I want the tree drawing to be 
    separated from 'results drawing'
    """
    pi = math.pi
    placements = image.add(image.g(id = 'placements'))

    # Prepare the code for the placement bubble, initial diameter 40, 
    # do not change, the bubble size is adjusted later
    radius = 20
    placements = create_aquabubblebody(placements, image, radius)
    placements = create_aquabubblebrilliance(placements, image, radius)
    # Done

    placement_bubbles = placements.add(image.g(id = 'placement_bubbles'))

    # We want to print the big bubbles first so that they cannot 
    # completely cover smaller ones
    total_weights = Autovivify()

    for node in sorted(tree.y_coord_of_node):
        if node == -1:
            continue
        if not mltreemap_results['counts_per_node'][node]['total_weight']:
            continue
        total_node_weight = mltreemap_results['counts_per_node'][node]\
          ['total_weight']
        total_weights[total_node_weight][node] = 1

    # Okay, now draw the bubble

    for total_node_weight in sorted(total_weights, reverse=True):

        for node in sorted(total_weights[total_node_weight]):
            # If total_node_weight == 1, 
            # radius_real == 20 * param_scale_bubble
            radius_real = math.sqrt(total_node_weight * 400) * \
              param_scale_bubble
            transform_factor = radius_real / radius

            # To print empty trees, uncomment below code
            # if transform_factor == 0:
            #     transform_factor = 1

            branch_length_of_node = tree.branch_length_of_node[node]
            x_coord_of_node_placement_center = tree.x_coord_of_node\
              [node] - branch_length_of_node / 2
            y_coord_of_node = tree.y_coord_of_node[node]

            # Transform the coordinates if we do the circular image
            if is_circular_image:
                x_coord_of_node_placement_center, y_coord_of_node = \
                  calculate_coordinates_circular(mltreemap_results, \
                  x_coord_of_node_placement_center, y_coord_of_node)
            # Done

            if not transform_factor:
                transform_factor = 1
            x_coord_of_node_placement = \
              (x_coord_of_node_placement_center - radius_real) * \
              (1 / transform_factor)
            y_coord_of_node_placement = \
              (y_coord_of_node - radius_real) * \
              (1 / transform_factor)

            # Color the bubble
            previous_end_angle = pi / 2

            for color in sorted(mltreemap_results['counts_per_node']\
              [node]['colors']):
                if mltreemap_results['counts_per_node'][node]['colors']\
                  [color]:
                    fraction = mltreemap_results['counts_per_node']\
                      [node]['colors'][color]
                if fraction < 0:
                    continue
                if total_node_weight <= 0:
                    continue
                only_one_hit = 0
                if fraction == total_node_weight:
                    only_one_hit = 1
                fract_of_bubble = fraction / total_node_weight
                start_angle = previous_end_angle
                if only_one_hit:
                    placement_bubbles.add(image.circle(center = \
                      (x_coord_of_node_placement_center, \
                      y_coord_of_node), r = radius_real, opacity = 1, \
                      fill = str(color)))
                else:
                    placement_bubbles, previous_end_angle = \
                      draw_pie_chart(placement_bubbles, \
                      x_coord_of_node_placement_center, \
                      y_coord_of_node, radius_real, start_angle, \
                      fract_of_bubble, color)

            if bubble_type:
                # Draw the bubble
                placement_bubbles.add(image.use('#body', insert = \
                  (x_coord_of_node_placement, \
                  y_coord_of_node_placement), transform = 'scale(' + \
                  str(transform_factor) + ')'))
                placement_bubbles.add(image.use('#brilliance', insert = \
                  (x_coord_of_node_placement, \
                  y_coord_of_node_placement), transform = 'scale(' + \
                  str(transform_factor) + ')'))

    return image, placements


def draw_pie_chart(placement_bubbles, center_x, center_y, radius, \
  alpha_r, fract_of_bubble, color):
    """TK"""
    # Get the coordinates

    # C (center_x|center_y) = center of circle; P (end_x|end_y), 
    # Q (start_x|start_y); alphatot = alpha_r + alpha
    # 
    # P---------Q
    # |alpha  / |
    # |     /   |
    # |alp/hatot|
    # | /alpha_r|
    # C----------
    # 

    pi = math.pi
    alpha = fract_of_bubble * 2 * pi
    alphatot = alpha + alpha_r
    x_start = center_x + math.cos(alpha_r) * radius
    y_start = center_y - math.sin(alpha_r) * radius
    x_end = center_x + math.cos(alphatot) * radius
    y_end = center_y - math.sin(alphatot) * radius

    # Next get the large_arc_flag
    large_arc_flag = 0
    if alpha > pi:
        large_arc_flag = 1

    # Draw the pie
    placement_bubbles.add(image.path(d = 'M ' + str(center_x) + ',' + \
      str(center_y) + ' L ' + str(x_start) + ',' + str(y_start) + \
      ' A ' + str(radius) + ',' + str(radius) + ' 0 ' + \
      str(large_arc_flag) + ', 0 ' + str(x_end) + ',' + str(y_end) + \
      ' z', opacity = 1, fill = str(color)))

    return placement_bubbles, alphatot


def draw_percents_and_placement_bars(image, placements, tree, \
  mltreemap_results, used_colors, color_mode, text_mode):
    """TK"""
    image_width = mltreemap_results['image']['width']
    tree_height = mltreemap_results['image']['tree_height']

    for node in sorted(tree.y_coord_of_node):
        if node <= 0:
            continue

        # First, get all the necessary information
        y_coord_of_node = tree.y_coord_of_node[node]
        x_coordinate_of_label_end = mltreemap_results\
          ['x_coordinate_of_label_end']
        highest_fraction_raw = mltreemap_results\
          ['highest_count_per_species'] * 100
        fraction_raw = 0

        # Prepare the text
        text = ''
        all_fractions_0 = 1
        fraction_total = 0.0

        for color in sorted(used_colors):
            # text as generated here is only used if color_mode == 1
            # (ie. show different datasets in different colors and
            # percentages for all of them)
            try:
                mltreemap_results['counts_per_species'][node][color]
            except NameError:
                pass
            else:
                fraction_raw = mltreemap_results['counts_per_species']\
                  [node][color] * 100
            fraction = (int(fraction_raw * 100 + 0.5)) / 100.0
            fraction_total += fraction
            if fraction > 0:
                all_fractions_0 = 0
            text += '(' + str('{0:.2f}'.format(fraction)) + '%)'

        fraction_total = int(fraction_total * 100 + 0.5) / 100.0
        if color_mode == 1:
            text = '(Total: ' + str(fraction_total) + '%)'
        if text_mode:
            fontsize = image_width / 125
            y_offset2 = 0.3 * fontsize
            x_gap = image_width / 400.0
            if all_fractions_0 != 0:
                continue
            placements.add(image.text(str(text), insert = \
              (x_coordinate_of_label_end - x_gap, y_coord_of_node + \
              y_offset2), style = 'font-family:Verdana; font-size:' + \
              str(fontsize) + '; text-anchor:end'))

        only_one_color = 1

        # Prepare the placement bars

        x_offset = 0

        for color in sorted(mltreemap_results['counts_per_species']\
          [node]):
            try:
                mltreemap_results['counts_per_species'][node][color]
            except NameError:
                pass
            else:
                fraction_raw = mltreemap_results['counts_per_species']\
                  [node][color] * 100
            if fraction_raw <= 0:
                continue
            y_min = tree.y_coord_of_node_min[node]
            y_max = tree.y_coord_of_node_max[node]
            height_tot = y_max - y_min
            height = height_tot * 0.9
            start_y = y_min + ((height_tot - height) / 2)
            x_gap2 = image_width / 200
            start_x = x_coordinate_of_label_end + x_gap2
            max_length = (image_width - image_width / 100) - start_x
            start_x += x_offset
            fractional_length = max_length * (fraction_raw / \
              highest_fraction_raw) # The highest fraction is max_length
            x_offset += fractional_length
            placements.add(image.rect(insert = (start_x, start_y), \
              size = (fractional_length, height), fill = str(color)))

    return image


def create_aquabubblebody(placements, image, radius):
    """TK"""
    svg = placements

    # Now the aquabubble body:

    # Part I
    bubble_gradient1 = image.defs.add(image.linearGradient(start = \
      ('50%', '100%'), end = ('50%', '0%'), id = 'bubble_gradient1'))
    svg.add(bubble_gradient1)
    bubble_gradient1.add_stop_color(offset = 0.2472, color = '#FAFAFA')
    bubble_gradient1.add_stop_color(offset = 0.3381, color = '#D0D0D0')
    bubble_gradient1.add_stop_color(offset = 0.4517, color = '#A2A2A2')
    bubble_gradient1.add_stop_color(offset = 0.5658, color = '#7C7C7C')
    bubble_gradient1.add_stop_color(offset = 0.6785, color = '#5E5E5E')
    bubble_gradient1.add_stop_color(offset = 0.7893, color = '#494949')
    bubble_gradient1.add_stop_color(offset = 0.8975, color = '#3C3C3C')
    bubble_gradient1.add_stop_color(offset = 1, color = '#383838')
    image.defs.add(image.circle(center = (radius, radius), r = radius, \
      id = 'bubble_part1', opacity = 0.4, fill = 'url(#bubble_gradient1)'))
    # Part I done

    # Part II
    bubble_gradient2 = image.defs.add(image.radialGradient(center = \
      ('50%', '50%'), r = '50%', id = 'bubble_gradient2'))
    svg.add(bubble_gradient2)
    bubble_gradient2.add_stop_color(offset = 0, color = '#FFFFFF')
    bubble_gradient2.add_stop_color(offset = 0.3726, color = '#FDFDFD')
    bubble_gradient2.add_stop_color(offset = 0.5069, color = '#F6F6F6')
    bubble_gradient2.add_stop_color(offset = 0.6026, color = '#EBEBEB')
    bubble_gradient2.add_stop_color(offset = 0.68, color = '#DADADA')
    bubble_gradient2.add_stop_color(offset = 0.7463, color = '#C4C4C4')
    bubble_gradient2.add_stop_color(offset = 0.805, color = '#A8A8A8')
    bubble_gradient2.add_stop_color(offset = 0.8581, color = '#888888')
    bubble_gradient2.add_stop_color(offset = 0.9069, color = '#626262')
    bubble_gradient2.add_stop_color(offset = 0.9523, color = '#373737')
    bubble_gradient2.add_stop_color(offset = 0.9926, color = '#090909')
    bubble_gradient2.add_stop_color(offset = 1, color = '#000000')

    image.defs.add(image.circle(center = (radius, radius), r = radius, \
      id ='bubble_part2', opacity = 0.1, fill = 'url(#bubble_gradient2)'))
    group1_svg = image.defs.add(image.g(id = 'body'))
    group1_svg.add(image.use('#bubble_part1'))
    group1_svg.add(image.use('#bubble_part2'))
    # Part II done

    # aquabubblebody done

    return placements


def create_aquabubblebrilliance(placements, image, radius):
    """TK"""
    svg = placements

    cy = radius / 2.5
    rx = radius / 2
    ry = radius / 4

    # Now the brilliance effect
    # This way of generating a opacity gradient generates valid SVG code
    # but cannot be interpreted by Adobe Illustrator CS3
    brilliance_gradient = image.defs.add(image.linearGradient(start = ('50%', '0%'), end = ('50%', '100%'), id = 'brilliance_gradient'))
    svg.add(brilliance_gradient)
    brilliance_gradient.add_stop_color(offset = 0, color = '#FFFFFF', \
      opacity = 1)
    brilliance_gradient.add_stop_color(offset = 0.1, color = '#FFFFFF',\
      opacity = 0.99)
    brilliance_gradient.add_stop_color(offset = 1, color = '#FFFFFF', \
      opacity = 0)

    image.defs.add(image.ellipse(center = (radius, cy), r = (rx, ry), \
      id = 'bubble_brilliance', opacity = 0.9, \
      fill = 'url(#brilliance_gradient)'))

    group1_svg = image.defs.add(image.g(id = 'brilliance'))
    svg.add(group1_svg)
    group1_svg.add(image.use('#bubble_brilliance'))

    return placements


def read_tree_topology(mltreemap_results):
    """TK"""
    print 'reading tree topology...'
    tree_file = 'tree_data/' + str(mltreemap_results['tree_file'])
    tree = NEWICK_tree()
    tree.read_tree_topology_from_file(tree_file)
    # TK Verified tree
    tree.optimize_species_display_order()

    # Do some hard coded manipulation
    tree.branch_length_of_node[-1] = 0.01
    # Done

    compute_node_positions(tree)
    # TK Verified all tree attributes
    tree = scale_node_positions(tree, mltreemap_results)
    # TK Verified tree.branch_length_of_node, tree.x_coord_of_node,
    #             tree.y_coord_of_node, tree.y_coord_of_node_min,
    #             tree.y_coord_of_node_max

#    sys.exit(str(tree.node_id_counter))
#    sys.exit(tree.parent_of_node)
#    sys.exit(tree.children_of_node)
#    sys.exit(tree.terminal_children_of_node)
#    sys.exit(tree.branch_length_of_node)
#    sys.exit(tree.terminal_nodes)
#    sys.exit(tree.bootstrap_support_of_node)
#    sys.exit(tree.x_coord_of_node)
#    sys.exit(tree.nodes_of_x_coords)
#    sys.exit(tree.y_coord_of_node)
#    sys.exit(tree.y_coord_of_node_min)
#    sys.exit(tree.y_coord_of_node_max)
#    sys.exit(tree.species_display_order)
 

    return tree


def scale_node_positions(tree, mltreemap_results):
    """TK"""
    image_width = mltreemap_results['image']['width']
    tree_height = mltreemap_results['image']['tree_height']
    y_offset = mltreemap_results['image']['y_offset']
    y_scaling_factor = tree_height # Note: This works because y values 
                                   #       range from 0 - 1

    # TK Verified image_width, tree_height, y_offset, y_scaling_factor

    # Assign the highest possible tree x coordinate 
    # and calculate the x_scaling_factor
    x_scaling_factor = 0.0

    for x_val in sorted(tree.nodes_of_x_coords, reverse=True):
        x_scaling_factor = float(mltreemap_results\
          ['x_coordinate_of_tree_end']) / float(x_val)
        try:
            x_scaling_factor
        except NameError:
            sys.exit('ERROR: x scaling factor could ' + \
              'not be calculated!\n')
        break

    if not x_scaling_factor:
        sys.exit('ERROR: x scaling factor could not be determined!\n')

    for node in sorted(tree.y_coord_of_node):
        y_coord_of_node = float(tree.y_coord_of_node[node]) * \
          float(y_scaling_factor) + float(y_offset)
        x_coord_of_node = float(tree.x_coord_of_node[node]) * \
          float(x_scaling_factor)
        branch_length_of_node = float(tree.branch_length_of_node[node]) * \
          float(x_scaling_factor)
        tree.y_coord_of_node[node] = y_coord_of_node
        tree.x_coord_of_node[node] = x_coord_of_node
        tree.branch_length_of_node[node] = branch_length_of_node
        if node > 0:
            y_coord_of_node_min = tree.y_coord_of_node_min[node] * \
              y_scaling_factor + y_offset
            y_coord_of_node_max = tree.y_coord_of_node_max[node] * \
              y_scaling_factor + y_offset
            tree.y_coord_of_node_min[node] = y_coord_of_node_min
            tree.y_coord_of_node_max[node] = y_coord_of_node_max

    return tree


def compute_node_positions_sort_helper(a, b):
    a_sort = 1
    if tree.species_display_order[a]:
        a_sort = tree.species_display_order[a]
    b_sort = 1
    if tree.species_display_order[b]:
        b_sort = tree.species_display_order[b]
    return cmp(a_sort, b_sort)


def compute_node_positions(tree):
    """
    This routine performs the actual placement of the tree (branches, 
    leaves) in terms of x,y-coordinate with which they will later be 
    put on the canvas. This routine does not perform the actual drawing.
    """
    y_position = 0.0

    for node in sorted(tree.parent_of_node, \
      key = lambda cmpNode:tree.species_display_order[cmpNode] if \
      cmpNode in tree.species_display_order else 1):
        if node <= 0:
            continue
        fraction = 1.0 / len(tree.terminal_nodes)
        # TK Verified fraction
        tree.y_coord_of_node_min[node] = y_position
        tree.y_coord_of_node[node] = y_position + (fraction * 0.5)
        tree.y_coord_of_node_max[node] = y_position + fraction
        y_position += fraction

    # TK Verified tree.y_coord_of_node
    # TK Verified tree.y_coord_of_node_min
    # TK Verified tree.y_coord_of_node_max
    assign_y_coord_internal_node(tree, -1)
    # TK Verified tree.y_coord_of_node
    # TK Verified tree.y_coord_of_node_min
    # TK Verified tree.y_coord_of_node_max
    assign_x_coord_to_node(tree, -1, 0.01) # Note: 0.01 is the hardcoded
                                           # branch length of node -1
    # TK Verified tree.x_coord_of_node
    # TK Verified tree.nodes_of_x_coords

    # CVM WATCH: The vertical position of node -2 is meddled with here 
    # in order to avoid the 'kink' at the root


def assign_y_coord_internal_node(tree, this_node):
    """
    This recursive routine assigns the Y-position for all internal 
    nodes. It also assigns positions needed for the 'taxon-illustration'
    bit. Call this routine only after you have already assigned 
    Y-positions to the terminal (leaf) nodes.
    """
    if tree.y_coord_of_node[this_node]:
        position = tree.y_coord_of_node[this_node]
        return position
    if this_node > 0: # For security. Should never happen.
        return 0
    min_position = 100000
    max_position = 0

    # TK Verified children_of_node[this_node]
    for child in tree.children_of_node[this_node]:
        position = assign_y_coord_internal_node(tree, child)
        if position < min_position:
            min_position = position
        if position > max_position and int(child) < 1000000:
            max_position = position

#    if this_node == -1:
#        print 'DEBUG'
#        sys.exit('min: ' + str(min_position))
#        sys.exit('max: ' + str(max_position))
    this_position = (min_position + max_position) / 2
    tree.y_coord_of_node[this_node] = this_position

    return this_position


def assign_x_coord_to_node(tree, this_node, current_x_position):
    """TK"""
    tree.x_coord_of_node[this_node] = current_x_position
    tree.nodes_of_x_coords[current_x_position] = this_node
    if this_node < 0:

        for child in tree.children_of_node[this_node]:
            branch_length = 0.02
            try:
                tree.branch_length_of_node[child]
            except NameError:
                print 'WARNING: No branch-length for node ' + \
                  str(this_node) + '!\n'
            else:
                branch_length = tree.branch_length_of_node[child]
            assign_x_coord_to_node(tree, child, \
              float(current_x_position) + float(branch_length))


def main(argv):
    print 'Start MLTreeMap Imagemaker TSN v. 0.0'
    #mltreemap_imagemaker = mltreemap_imagemaker()
    args = getParser()
    user_options = read_user_input(args)
    input_files = get_input_files(user_options)
    # TK Verified input_files
    concatenated_input_files, text_of_denominator, used_colors = concatenate_files(user_options, input_files)
    # TK Verified concatenated_input_files
    # TK Verified text_of_denominator
    # TK Verified text_of_denominator
    run_the_imagemaker(user_options, concatenated_input_files, text_of_denominator, used_colors)
    print 'Done.'


if __name__ == "__main__":
    main(sys.argv[1:])
