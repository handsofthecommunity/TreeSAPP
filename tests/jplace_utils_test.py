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

TEST_DATA_PATH ='/home/travis/build/hallamlab/TreeSAPP/tests/test_data'
TREESAPP_PATH = '/home/travis/build/hallamlab/TreeSAPP/'
TEST_DIR = '/home/travis/build/hallamlab/TreeSAPP/tests'

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

        #NOTE: placements do not match
        
   def test_sub_indices_for_seq_names_jplace(self):
        short_numeric_contig_index = {'McrA': {-12: 'PHP46140.1_methyl-coenzyme_M_reductase_subunit_alpha_Methanosarcinales_archaeon_ex4572_44_4_595', -2: 'OYT62528.1_hypothetical_protein_B6U67_04395_Methanosarcinales_archaeon_ex4484_138_1_471', -10: 'AAM30936.1_Methyl-coenzyme_M_reductase__alpha_subunit_Methanosarcina_mazei_Go1_1_570' , -9: 'AAU83782.1_methyl_coenzyme_M_reductase_subunit_alpha_uncultured_archaeon_GZfos33H6_3_570', -8 : 'ADN36741.1_methyl-coenzyme_M_reductase__alpha_subunit_Methanolacinia_petrolearia_DSM_11571_2_568', -7: 'AUD55425.1_methyl-coenzyme_M_reductase_alpha_subunit__partial_uncultured_euryarchaeote_1_563', -6: 'OFV67773.1_methyl_coenzyme_M_reductase_subunit_alpha_Candidatus_Syntrophoarchaeum_caldarius_1_561', -5: 'KUE73676.1_methyl-coenzyme_M_reductase_subunit_alpha_Candidatus_Methanomethylophilus_sp._1R26_1_554', -4: 'PKL66143.1_coenzyme-B_sulfoethylthiotransferase_subunit_alpha_Methanobacteriales_archaeon_HGW-Methanobacteri_1_552', -3: 'AAU82491.1_methyl_coenzyme_M_reductase_I_subunit_alpha_uncultured_archaeon_GZfos18B6_1_501', -1: 'AFD09581.1_methyl-coenzyme_M_reductase_alpha_subunit__partial_uncultured_Methanomicrobiales_archaeon_1_254', -11: 'PKL62129.1_methyl-coenzyme_M_reductase_subunit_alpha_Methanomicrobiales_archaeon_HGW-Methanomicrobiales-2_1_580'}}
        
        args = create_parser(TREESAPP_PATH, 'M0701', 'p')
        marker_build_dict = treesapp.parse_ref_build_params(args)
        marker_build_dict = treesapp.parse_cog_list(args, marker_build_dict)

        jplace_path = TEST_DIR + "test_data/tmp.jplace"
        jplace_data = jplace_utils.jplace_parser(jplace_path)

        jplace_ref_data = TEST_DIR + "test_data/ref_jplace.jplace"
        jplace_out = args.output_dir_var + 'RAxML_portableTree.M0701_hmm_purified_group0-BMGE-qcd.phy.jplace'
        copyfile(jplace_ref_data, jplace_out)
 
        jplace_utils.sub_indices_for_seq_names_jplace(args, short_numeric_contig_index, marker_build_dict)       
        test_data = jplace_utils.jplace_parser(jplace_out)
        for i in range(len(test_data.placements)):
            assert(jplace_data.placements[i] == test_data.placements[i])

    @pytest.mark.dependency(name="test_jplace_parser")
    def test_jplace_parser(self):
        jplace_test_file = TEST_DIR + 'test_data/tmp.jplace'
        itol = jplace_utils.jplace_parser(jplace_test_file)

        assert(itol.version == 2)
        assert(len(itol.placements) == 12)

        for val in itol.placements:
            assert('n' in val.keys())
            assert('p' in val.keys())
