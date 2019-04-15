__author__ = 'Connor Morgan-Lang'

import sys
import re
import os
import logging

import _fasta_reader
from .utilities import median, reformat_string
from .external_command_interface import launch_write_command


# No bioinformatic software would be complete without a contribution from Heng Li.
# Adapted from his readfq generator
def generate_fasta(fasta_handler):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:
            for line in fasta_handler:  # search for the start of the next record
                if line[0] == '>':  # fasta header line
                    last = line[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:], [], None
        for line in fasta_handler:  # read the sequence
            if line[0] == '>':
                last = line[:-1]
                break
            seqs.append(line[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs)  # yield a fasta record
            if not last:
                break
        else:
            seq, seqs = ''.join(seqs), []
            for line in fasta_handler:  # read the quality
                seqs.append(line[:-1])
            if last:  # reach EOF before reading enough quality
                yield name, seq  # yield a fasta record instead
                break


def read_fasta_to_dict(fasta_file):
    """
    Reads any fasta file using a generator function (generate_fasta) into a dictionary collection

    :param fasta_file: Path to a FASTA file to be read into a dict
    :return: Dict where headers/record names are keys and sequences are the values
    """
    fasta_dict = dict()
    try:
        fasta_handler = open(fasta_file, 'r')
    except IOError:
        logging.error("Unable to open " + fasta_file + " for reading!\n")
        sys.exit(5)
    for record in generate_fasta(fasta_handler):
        name, sequence = record
        fasta_dict[name] = sequence.upper()
    return fasta_dict


class Header:
    def __init__(self, header):
        self.original = header
        self.formatted = ""
        self.treesapp_name = ""
        self.post_align = ""
        self.first_split = ""

    def get_info(self):
        info_string = "TreeSAPP ID = '" + self.treesapp_name + "'\tPrefix = '" + self.first_split + "'\n"
        info_string += "Original =  " + self.original + "\nFormatted = " + self.formatted
        return info_string


def register_headers(header_list):
    acc = 1
    header_registry = dict()
    for header in header_list:
        new_header = Header(header)
        new_header.formatted = reformat_string(header)
        new_header.first_split = header.split()[0]
        header_registry[str(acc)] = new_header
        acc += 1
    return header_registry


class FASTA:
    def __init__(self, file_name):
        self.file = file_name
        self.fasta_dict = dict()
        self.header_registry = dict()
        self.amendments = dict()

    def synchronize_seqs_n_headers(self):
        if len(self.fasta_dict.keys()) != len(self.header_registry):
            sync_header_registry = dict()
            excluded_headers = list()
            for num_id in self.header_registry:
                header = self.header_registry[num_id]
                if header.formatted not in self.fasta_dict and header.original not in self.fasta_dict:
                    excluded_headers.append(self.header_registry[num_id].original)
                else:
                    sync_header_registry[num_id] = self.header_registry[num_id]
            self.header_registry = sync_header_registry

    def summarize_fasta_sequences(self):
        num_headers = 0
        longest = 0
        shortest = 0
        sequence_lengths = []
        for name in self.fasta_dict:
            sequence = self.fasta_dict[name]
            # Calculate stats
            num_headers += 1
            if longest == 0 and shortest == 0:
                longest = len(sequence)
                shortest = len(sequence)
            elif len(sequence) > longest:
                longest = len(sequence)
            elif len(sequence) < shortest:
                shortest = len(sequence)
            else:
                pass
            sequence_lengths.append(len(sequence))

        stats_string = "\tNumber of sequences: " + str(num_headers) + "\n"
        stats_string += "\tLongest sequence length: " + str(longest) + "\n"
        stats_string += "\tShortest sequence length: " + str(shortest) + "\n"
        stats_string += "\tMean sequence length: " + str(round(sum(sequence_lengths) / num_headers, 1)) + "\n"
        stats_string += "\tMedian sequence length: " + str(median(sequence_lengths)) + "\n"

        logging.info(stats_string)
        return

    def deduplicate_fasta_sequences(self):
        """
        Removes exact duplicate sequences from a FASTA-formatted dictionary of sequence records
        """
        dedup_dict = dict()
        if len(set(self.fasta_dict.values())) == len(self.fasta_dict):
            return
        else:
            for header in self.fasta_dict.keys():
                if self.fasta_dict[header] not in dedup_dict.values():
                    dedup_dict[header] = self.fasta_dict[header]
            self.fasta_dict = dedup_dict
        self.synchronize_seqs_n_headers()
        return

    def use_first_split(self):
        first_split_dict = dict()
        for acc in self.header_registry:
            header = self.header_registry[acc]  # type: Header
            if header.original in self.fasta_dict:
                first_split_dict[header.first_split] = self.fasta_dict[header.original]
            elif header.formatted in self.fasta_dict:
                first_split_dict[header.first_split] = self.fasta_dict[header.formatted]
            else:
                pass
        self.fasta_dict = first_split_dict
        return

    def update(self, fasta, file=True):
        # Format the inputs
        if file:
            if not os.path.isfile(fasta):
                logging.error("File '" + fasta + "' does not exist!\n")
                sys.exit(13)
            new_fasta = read_fasta_to_dict(fasta)
            new_fasta_headers = register_headers(get_headers(fasta))
        else:
            new_fasta = fasta
            new_fasta_headers = register_headers(fasta.keys())

        self.amendments.update(new_fasta_headers)
        # Load the new fasta and headers
        for seq_name in new_fasta:
            self.fasta_dict[seq_name] = new_fasta[seq_name]
        acc = max([int(x) for x in self.header_registry.keys()])
        for num_id in sorted(new_fasta_headers, key=int):
            acc += 1
            self.header_registry[str(acc)] = new_fasta_headers[num_id]
        return


def format_read_fasta(fasta_input, molecule, output_dir, max_header_length=110, min_seq_length=10):
    """
    Reads a FASTA file, ensuring each sequence and sequence name is valid.

    :param fasta_input: Absolute path of the FASTA file to be read
    :param molecule: Molecule type of the sequences ['prot', 'dna', 'rrna']
    :param output_dir: Path to a directory for writing the log file to
    :param max_header_length: The length of the header string before all characters after this length are removed
    :param min_seq_length: All sequences shorter than this will not be included in the returned list.
    :return: A Python dictionary with headers as keys and sequences as values
    """

    if sys.version_info > (2, 9):
        py_version = 3
    else:
        py_version = 2
        from itertools import izip

    fasta_list = _fasta_reader._read_format_fasta(fasta_input,
                                                  min_seq_length,
                                                  output_dir,
                                                  molecule,
                                                  max_header_length)
    if not fasta_list:
        sys.exit(5)
    tmp_iterable = iter(fasta_list)
    if py_version == 2:
        formatted_fasta_dict = dict(izip(tmp_iterable, tmp_iterable))
    elif py_version == 3:
        formatted_fasta_dict = dict(zip(tmp_iterable, tmp_iterable))
    else:
        logging.error("Unexpected Python version detected: " + str(py_version), )
        sys.exit(5)

    for header in formatted_fasta_dict.keys():
        if len(header) > max_header_length:
            logging.error(header + " is too long (" + str(len(header)) + ")!\n" +
                          "There is a bug in _read_format_fasta - please report!\n")
            sys.exit(5)

    return formatted_fasta_dict


def get_headers(fasta_file):
    """
    Reads a FASTA file and returns a list of all headers it found in the file. No reformatting or filtering performed.
    :param fasta_file: Path to the FASTA file to be read.
    :return:
    """
    original_headers = list()
    try:
        fasta = open(fasta_file, 'r')
    except IOError:
        logging.error("Unable to open the FASTA file '" + fasta_file + "' for reading!")
        sys.exit(5)

    n_headers = 0
    for name, _ in generate_fasta(fasta):
        n_headers += 1
        original_headers.append('>' + str(name))

    fasta.close()
    logging.debug("Read " + str(n_headers) + " headers from " + fasta_file + ".\n")

    return original_headers


def write_new_fasta(fasta_dict, fasta_name, max_seqs=None, headers=None):
    """
    Function for writing sequences stored in dictionary to file in FASTA format; optional filtering with headers list

    :param fasta_dict: A dictionary containing headers as keys and sequences as values
    :param fasta_name: Name of the FASTA file to write to
    :param max_seqs: If not None, the maximum number of sequences to write to a single FASTA file
    :param headers: Optional list of sequence headers. Only fasta_dict keys in headers will be written
    :return:
    """
    split_files = list()
    file_counter = 1
    sequence_accumulator = 0
    fasta_string = ""

    # Check for '>' leading sequence names. Strip them if present.
    if headers:
        side_chevy = headers[0][0] == '>'
        for header in headers:
            state = header[0] == '>'
            if state is not side_chevy:
                logging.error("Inconsistent header names in headers list object\n")
                sys.exit(5)
        if side_chevy:
            headers = [header[1:] for header in headers]

    if max_seqs is not None:
        fasta_name = fasta_name + '_' + str(file_counter) + ".fasta"

    try:
        fa_out = open(fasta_name, 'w')
    except IOError:
        logging.error("Unable to open " + fasta_name + " for writing!\n")
        sys.exit(5)

    for name in sorted(fasta_dict.keys()):
        seq = fasta_dict[name]
        if name[0] != '>':
            name = '>' + name
        sequence_accumulator += 1
        if max_seqs and sequence_accumulator > max_seqs:
            # If input is to be split and number of sequences per file has been exceeded begin writing to new file
            fa_out.write(fasta_string)
            fa_out.close()
            fasta_string = ""
            split_files.append(fasta_name)
            file_counter += 1
            sequence_accumulator = 1
            fasta_name = re.sub("_\d+.fasta$", '_' + str(file_counter) + ".fasta", fasta_name)
            fa_out = open(fasta_name, 'w')

        # Only write the records specified in `headers`, or all records if `headers` isn't provided
        if headers is None:
            fasta_string += name + "\n" + seq + "\n"
        elif name[1:] in headers:
            fasta_string += name + "\n" + seq + "\n"
        else:
            sequence_accumulator -= 1

    fa_out.write(fasta_string)
    fa_out.close()
    split_files.append(fasta_name)

    return split_files


def get_header_format(header, code_name=""):
    """
    Used to decipher which formatting style was used and parse information, ideally reliably
    HOW TO ADD A NEW REGULAR EXPRESSION:
        1. create a new compiled regex pattern, like below
        2. add the name of the compiled regex pattern to the header_regexes dictionary
        3. if the regex groups are new and complicated (parsing more than the accession and organism info),
        alter return_sequence_info_groups in create_treesapp_ref_data to add another case

    :param header: A sequences header from a FASTA file
    :param code_name:
    :return:
    """
    # The regular expressions with the accession and organism name grouped
    # Protein databases:
    gi_re = re.compile(">gi\|(\d+)\|[a-z]+\|[_A-Z0-9.]+\|.* RecName: Full=([A-Za-z1-9 _\-]+);?.*$")  # a
    gi_prepend_proper_re = re.compile(">gi\|\d+\|[a-z]{2,4}\|([_A-Z0-9.]+)\| (.*) \[(.*)\]$")  # a, d, o
    gi_prepend_mess_re = re.compile(">gi\|(\d+)\|[a-z]{2,4}\|.*\|([\w\s.,\-\()]+)$")  # a
    dbj_re = re.compile(">dbj\|(.*)\|.*\[(.*)\]")  # a, o
    emb_re = re.compile(">emb\|(.*)\|.*\[(.*)\]")
    gb_re = re.compile(">gb\|(.*)\|.*\[(.*)\]")
    ref_re = re.compile(">ref\|(.*)\|.*\[(.*)\]")
    pdb_re = re.compile(">pdb\|(.*)\|.+$")  # a
    pir_re = re.compile(">pir\|.*\|(\w+).* - (.*)$")  # a, o
    presf_re = re.compile(">prf\|.*\|([A-Z0-9]+)\s+.*$")  # a
    sp_re = re.compile(">sp\|(.*)\|.*Full=.*;?.*$")  # a
    fungene_gi_bad = re.compile("^>[0-9]+\s+coded_by=.+,organism=.+,definition=.+$")
    mltree_re = re.compile("^>(\d+)_" + re.escape(code_name) + "$")
    pfam_re = re.compile("^>([A-Za-z0-9_|]+)/[0-9]+-[0-9]+$")  # a

    # Nucleotide databases:
    # silva_arb_re = re.compile("^>([A-Z0-9]+)\.([0-9]+)\.([0-9]+)_(.*)$")
    # refseq_nuc_re = re.compile("^>([A-Z]+_[0-9]+\.[0-9])_.+$")  # a
    # nr_re = re.compile("^>([A-Z0-9]+\.[0-9])_.*$")  # a

    # Ambiguous:
    # genbank_exact_genome = re.compile("^>([A-Z]{1,2}[0-9]{5,6}\.?[0-9]?) .* \[(.*)\]$")  # a, o
    accession_only = re.compile(r"^>([A-Z]+_?[0-9]+\.?[0-9]?)$")  # a
    ncbi_ambiguous = re.compile(r"^>([A-Za-z0-9.\-_]+)\s+.*(?<!])$")  # a
    ncbi_org = re.compile(r"^>([A-Z0-9.\-_]+\.?[0-9]?)\s+(?!lineage=).*\[.*\]$")  # a

    # Custom fasta header with taxonomy:
    # First group = contig/sequence name, second = full taxonomic lineage, third = description for tree
    # There are no character restrictions on the first and third groups
    # The lineage must be formatted like:
    #   cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria
    custom_tax = re.compile(r"^>(.*) lineage=([A-Za-z ]+; .*) \[(.*)\]$")  # a, l, o

    header_regexes = {"prot": {dbj_re: "dbj",
                               emb_re: "emb",
                               gb_re: "gb",
                               pdb_re: "pdb",
                               pir_re: "pir",
                               ref_re: "ref",
                               sp_re: "sp",
                               gi_re: "gi_re",
                               gi_prepend_proper_re: "gi_proper",
                               gi_prepend_mess_re: "gi_mess",
                               pfam_re: "pfam",
                               presf_re: "prf"},
                      "dna": {mltree_re: "mltree"},
                      "ambig": {accession_only: "bare",
                                ncbi_ambiguous: "ncbi_ambig",
                                ncbi_org: "ncbi_org",
                                custom_tax: "custom"}
                      }

    if fungene_gi_bad.match(header):
        logging.warning(header + " uses GI numbers which are now unsupported by the NCBI! " +
                        "Consider switching to Accession.Version identifiers instead.\n")

    header_format_re = None
    header_db = None
    header_molecule = None
    format_matches = list()
    for molecule in header_regexes:
        for regex in header_regexes[molecule]:
            if regex.match(header):
                header_format_re = regex
                header_db = header_regexes[molecule][regex]
                header_molecule = molecule
                format_matches.append(header_db)
            else:
                pass
    if len(format_matches) > 1:
        logging.error("Header '" + header + "' matches multiple potential formats:\n\t" +
                      ", ".join(format_matches) + "\n" +
                      "TreeSAPP is unable to parse necessary information properly.\n")
        sys.exit(5)

    if header_format_re is None:
        logging.error("Unable to parse header '" + header + "'\n")
        sys.exit(5)

    return header_format_re, header_db, header_molecule


def summarize_fasta_sequences(fasta_file):
    try:
        fasta_handler = open(fasta_file, 'r')
    except IOError:
        logging.error("Unable to open " + fasta_file + " for reading!\n")
        sys.exit(5)

    num_headers = 0
    longest = 0
    shortest = 0
    sequence = None
    sequence_lengths = []
    line = fasta_handler.readline()
    while line:
        if line[0] == '>':
            num_headers += 1
            if sequence is not None:
                if longest == 0 and shortest == 0:
                    longest = len(sequence)
                    shortest = len(sequence)
                if len(sequence) > longest:
                    longest = len(sequence)
                if len(sequence) < shortest:
                    shortest = len(sequence)
                else:
                    pass
                sequence_lengths.append(len(sequence))
            sequence = ""
        else:
            sequence += line.strip()
        line = fasta_handler.readline()
    # Log the last sequence in the file
    if longest == 0 and shortest == 0:
        longest = len(sequence)
        shortest = len(sequence)
    if len(sequence) > longest:
        longest = len(sequence)
    if len(sequence) < shortest:
        shortest = len(sequence)
    sequence_lengths.append(len(sequence))

    stats_string = "\tNumber of sequences: " + str(num_headers) + "\n"
    stats_string += "\tLongest sequence length: " + str(longest) + "\n"
    stats_string += "\tShortest sequence length: " + str(shortest) + "\n"
    stats_string += "\tMean sequence length: " + str(round(sum(sequence_lengths)/num_headers, 1)) + "\n"
    stats_string += "\tMedian sequence length: " + str(median(sequence_lengths)) + "\n"

    logging.info(stats_string)
    return


def trim_multiple_alignment(executable, mfa_file, molecule, tool="BMGE"):
    """
    Trims the multiple sequence alignment using either BMGE or trimAl

    :param executable: Path to the executable of `tool`
    :param mfa_file: Name of a MSA file
    :param molecule: prot | dna
    :param tool: Name of the software to use for trimming [BMGE|trimAl]
    Returns file name of the trimmed multiple alignment file in FASTA format
    """
    f_ext = mfa_file.split('.')[-1]
    if not re.match("mfa|fasta|phy|fa", f_ext):
        logging.error("Unsupported file format: '" + f_ext + "'\n")
        sys.exit(5)

    trimmed_msa_file = re.sub('.' + re.escape(f_ext), '-' + re.escape(tool) + ".fasta", mfa_file)
    if tool == "trimAl":
        trim_command = [executable]
        trim_command += ['-in', mfa_file,
                         '-out', trimmed_msa_file,
                         '-gt', str(0.02)]
    elif tool == "BMGE":
        if molecule == "prot":
            bmge_settings = ["-t", "AA", "-m", "BLOSUM30"]
        else:
            bmge_settings = ["-t", "DNA", "-m", "DNAPAM100:2"]
        trim_command = ["java", "-jar", executable]
        trim_command += bmge_settings
        trim_command += ["-g", "0.99:0.33"]  # Specifying the gap rate per_sequence:per_character
        trim_command += ['-i', mfa_file,
                         '-of', trimmed_msa_file]
    else:
        logging.error("Unsupported trimming software requested: '" + tool + "'")
        sys.exit(5)

    logging.debug("STAGE: Multiple sequence alignment trimming\n" +
                  "\tINPUT: " + mfa_file + "\n" +
                  "\tTOOL: " + tool + "\n" +
                  "\tCOMMAND:\n" + " ".join(trim_command) + "\n" +
                  "\tOUTPUT:\n")

    stdout, return_code = launch_write_command(trim_command)
    if return_code != 0:
        logging.error(tool + " did not complete successfully!\n" +
                      tool + "output:\n" + stdout + "\n")
        sys.exit(5)

    return trimmed_msa_file