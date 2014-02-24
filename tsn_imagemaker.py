#/usr/bin/python

# MLTreeMap Imagemaker TSN v. 0.0

# Include the necessary libraries
try:

except:


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


def getParser():
    """TK"""
    parser = argparse.ArgumentParser(description='example input:\n\
                      python tsn_imagemaker.py -i example_input/\n\n\
                      mandatory input option:\n\
                      -i: either a) a path to a directory with input files \
                      or b) to a specific input file.\n\
                      note that in case a) all input files of the same type \
                      (i.e. p, g, n, h, r, m, d) will be concatenated.\n\n\
                      optional input options:\n\
                      -o: output directory (default: output/)\n\
                      -d: use different colors for different datasets \
                      (0: don\'t use this mode (default),\
                       1: use this mode,\
                       11: use this mode but concatenate percent values \
                      over datasets). \
                      See the user guide for a more detailed description.\n\
                      -t: display trees inclusive text and labels \
                      (1, default) \
                      or without text (0)\n\
                      -b: parameter for the size of the placement bubbles \
                      (default: 1)\n\
                      -r: display 16s and 18s rRNA hits in different trees \
                      (2, default) \
                      or one tree of life (1)\n\
                      -h: high quality placement bubbles (1, default) or \
                      low quality bubbles (0). \
                      High quality bubbles work perfectly with inkscape \
                      but can cause trouble with Adobe Illustrator.')
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
    parser.add_argument('-h', '--hq_bubbles', default=1, choices=[0,1], \
                        help='high quality bubbles work perfectly with \
                        inkscape but can cause trouble with Adobe \
                        Illustrator')
    parser.add_argument('-i', '--input', required=True, default=0, type=int, choices=[0, 's'], help='')
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

    input2 = PATHDELIM + input2

    # If input is an individual file, get file name
    if re.search(r'.*/((.)_.*RAxML_.+.txt)\Z', input2):
        try:
            in = open(input2, 'read')
        except: #TK check for exception type
            sys.exit('ERROR: Unable to open ' + str(input2) + '!\n')
        temp = re.search(r'.*/((.)_.*RAxML_.+.txt)\Z', input2)
        filename = temp.group(1)
        denominator = temp.group(2)
        if denominator in ['b','a'] and rRNA_display_option == 1:
            denominator = 'A'
        input_files[denominator][input] = filename

    # Else, if input is a directory, get file names of contents
    else:
        try:
            files = [f for f in os.listdir(input2) \
                     if os.path.isfile(join(input2,f))]
        except: # TK
            sys.exit('ERROR: Your input directory ' + str(input2) + \
                     'does not exist!\n')

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
    except: # TK
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
        colors[count] = color_tag
        count += 1

    nr_of_colors = len(colors)

    for denominator in sorted(input_files):
        output_filename = str(denominator) + '_concatenated_RAxML_outputs.txt'
        output_filename_long = str(output_dir) + str(output_filename)
        try:
            output = open(output_filename_long, 'w')
        except: # TK
            sys.exit('ERROR: Can\'t create ' + str(output_filename_long) + '\n')
        if args.colors > 0:
            try:
                output2 = open(str(output_dir) + str(denominator) + '_color_legend.txt', 'w')
            except: # TK
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
                file = open(filename_long, 'w')
            except: # TK
                sys.exit('ERROR: Can\'t open ' + str(filename_long))

            for fileLine in file:
                fileLine = fileLine.strip()
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
                try:
                    percentages_of_texts[text]
                except NameError:
                    percentages_of_texts[text] = weight
                else:
                    percentages_of_texts[text] += weight
            file.close()

        if args.colors > 0:
            output2.close()
        check = 0

        for text in sorted(percentages_of_texts):
            weight = percentages_of_texts[text]
            relative_weight = weight / nr_of_files
            relative_weight = (int(relative_weight * 10000 + 0.5)) / 10000
            output.write('Placement weight ' + str(relative_weight) + '%:' + str(text) + '\n')
            check += relative_weight

        print str(denominator) + ' files: sum of percentages = ' + str(check) + '\n'
        output.close()

    return concatenated_input_files, text_of_denominator, used_colors


def main(argv):
    print 'Start MLTreeMap Imagemaker TSN v. 0.0'


if __name__ == "__main__":
    main(sys.argv[1:])
