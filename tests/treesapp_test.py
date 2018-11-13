import pytest
import unittest
import os
import argparse

import sys, inspect
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))))
import file_parsers
import fasta
import treesapp

def create_parser(treesapp, targets, reftree):
    args = argparse.Namespace()
    args.alignment_mode = 'd'
    args.reftree = reftree
    args.targets = [targets]
    args.treesapp = treesapp
    args.check_trees = False
    args.fasta_input = 'test_data/marker_test_suite.faa'
    args.output = '/home/travis/build/hallamlab/marker_test/'
    args.skip = 'n'
    args.molecule = 'prot'
    return args

class TreeSAPPTest(unittest.TestCase):

    def test_get_alignment_data(self):
        args = create_parser('/home/travis/build/hallamlab/TreeSAPP/', 'M0701', 'p')
        marker_build_dict = file_parsers.parse_ref_build_params(args)
        alignment_dimensions_dict = treesapp.get_alignment_dims(args, marker_build_dict)
        assert(alignment_dimensions_dict['M0701'] == (214, 837))

    def test_hmmsearch_orfs_parse_domain_tables(self):
        args = create_parser('/home/travis/build/hallamlab/TreeSAPP/', 'M0701', 'p')
        marker_build_dict = file_parsers.parse_ref_build_params(args)
        marker_build_dict = file_parsers.parse_cog_list(args, marker_build_dict)
        hmm_domtbl_files = treesapp.hmmsearch_orfs(args, marker_build_dict)
        assert(hmm_domtbl_files[0] == '/home/travis/build/hallamlab/marker_test/various_outputs/McRA_to_ORFs_domtbl.txt')

        args = create_parser('~', 'ALL', 'p')
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            treesapp.hmmsearch_orfs(args, marker_build_dict)
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 3

    def test_extract_hmm_matches(self):
        args = create_parser('/home/travis/build/hallamlab/TreeSAPP/', 'M0701', 'p')
        args.formatted_input_file = args.output_dir_var + args.fasta_input + "_formatted.fasta"
        marker_build_dict = fasta.parse_ref_build_params(args)
        marker_build_dict = fasta.parse_cog_list(args, marker_build_dict)
        formatted_fasta_dict = fasta.format_read_fasta(args.fasta_input, "prot", args.output)

        fasta.write_new_fasta(formatted_fasta_dict, args.formatted_input_file)

        hmm_domtbl_files = treesapp.hmmsearch_orfs(args, marker_build_dict)
        hmm_matches = treesapp.parse_domain_tables(args, hmm_domtbl_files)
        homolog_seq_files, numeric_contig_index = treesapp.extract_hmm_matches(args, hmm_matches, formatted_fasta_dict)


        assert(homolog_seq_files == '/home/travis/build/hallamlab/marker_test/various_outputs/McrA_hmm_purified_group0.faa')
        assert('McrA' in numeric_contig_index.keys())
        assert(len(numeric_contig_index['McrA']) == 12)
        assert(numeric_contig_index['McrA'].keys() == [-12, -2, -10, -9, -8, -7, -6, -5, -4, -3, -1, -11])

