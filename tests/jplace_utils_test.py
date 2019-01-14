import pytest
import unittest
import os
import argparse
import glob

from .treesapp_test import create_parser

import sys, inspect
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))))
import jplace_utils
import treesapp

TEST_DATA_PATH='/home/travis/build/hallamlab/TreeSAPP/tests/test_data'                              TREESAPP_PATH = '/home/travis/build/hallamlab/TreeSAPP/'

class JplaceUtilsTests(unittest.TestCase):

    def test_organize_jplace_files(self):
        args = create_parser(TREESAPP_PATH, 'M0701', 'p')
        jplace_files = glob.glob(TEST_DIR + 'test_data' + os.sep + '*.jplace')

        jplace_collection = jplace_utils.organize_jplace_files(jplace_files)

        assert(jplace_collection['M0701'] == jplace_files)

        jplace_files = [TEST_DIR + '/*.jplace']
        with pytest.raises(Exception):
            jplace_utils.organize_jplace_files(jplace_files)

    def test_sub_indices_for_seq_names_jplace(self):
        short_numeric_contig_index = {'McrA': {-12: 'PHP46140.1_methyl-coenzyme_M_reductase_subunit_alpha_Methanosarcinales_archaeon_ex4572_44_8_595', -2: 'OYT62528.1_hypothetical_protein_B6U67_04395_Methanosarcinales_archaeon_ex4484_138_1_471', -10: 'AAU83782.1_methyl_coenzyme_M_reductase_subunit_alpha_uncultured_archaeon_GZfos33H6_1_570'}}

        args = create_parser(TREESAPP_PATH, 'M0701', 'p')
        marker_build_dict = treesapp.parse_ref_build_params(args)
        marker_build_dict = treesapp.parse_cog_list(args, marker_build_dict)

        
    def test_jplace_parser(self):
        jplace_files = glob.glob(TEST_DIR + 'test_data' + os.sep + '*.jplace')
        itol = jplace_utils.jplace_parser(jplace_files[0])
        
