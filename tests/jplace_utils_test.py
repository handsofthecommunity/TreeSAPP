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
import file_parsers

TEST_DATA_PATH ='/home/travis/build/hallamlab/TreeSAPP/tests/test_data'
TREESAPP_PATH = '/home/travis/build/hallamlab/TreeSAPP/'
TEST_DIR = '/home/travis/build/hallamlab/TreeSAPP/tests/'

class JplaceUtilsTests(unittest.TestCase):

    @pytest.mark.dependency(depends={"test_jplace_parser"})    
    def test_write_jplace(self):
        jplace_file = TEST_DIR + "test_data/test_write_jplace.jplace"
        jplace_original = TEST_DIR + 'test_data/tmp.jplace'
        itol = jplace_utils.jplace_parser(TEST_DATA_PATH + '/tmp.jplace')

        jplace_utils.write_jplace(itol, jplace_file)

        itol_check = jplace_utils.jplace_parser(TEST_DATA_PATH + '/ref_jplace.jplace')

        assert(itol.version == itol_check.version)
        assert(len(itol.placements) == len(itol_check.placements))

        assert(itol.metadata == itol_check.metadata)
        assert(len(itol.fields) == len(itol_check.fields))
        assert(itol.tree == itol_check.tree)


    @pytest.mark.dependency(name="test_jplace_parser")
    def test_jplace_parser(self):
        jplace_test_file = TEST_DIR + 'test_data/tmp.jplace'
        itol = jplace_utils.jplace_parser(jplace_test_file)

        assert(itol.version == 2)
        assert(len(itol.placements) == 12)

        for val in itol.placements:
            assert('n' in val.keys())
            assert('p' in val.keys())
