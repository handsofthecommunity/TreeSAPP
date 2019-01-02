import pytest
import unittest
import os
import argparse

import sys, inspect
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))))
import file_parsers

TEST_PATH = '/home/travis/build/halamlab/TreeSAPP/tests'
TREESAPP_PATH = '/home/travis/build/hallamlab/TreeSAPP/'

def create_parser(treesapp,targets, reftree):
        args = argparse.Namespace()
        args.alignment_mode = 'd'
        args.reftree = reftree
        args.targets = [targets]
        args.treesapp = treesapp
        args.check_trees = False
        args.skip = 'n'
        args.molecule = 'prot'
        return args

class ParserTest(unittest.TestCase):

    def check_parse_ref_build_params_out(self, expected_vals, targets, marker_build_dict):
        for i in range(0, len(targets)):
                assert (marker_build_dict[targets[i]].cog == expected_vals[i][0])
                assert (marker_build_dict[targets[i]].denominator == expected_vals[i][1])
                assert (marker_build_dict[targets[i]].molecule == expected_vals[i][2])
                assert (marker_build_dict[targets[i]].model == expected_vals[i][3])
                assert (marker_build_dict[targets[i]].pid == expected_vals[i][4])
                assert (marker_build_dict[targets[i]].lowest_confident_rank == expected_vals[i][5])
                assert (marker_build_dict[targets[i]].update == expected_vals[i][11])
                assert (marker_build_dict[targets[i]].description == expected_vals[i][6])
                assert (marker_build_dict[targets[i]].kind == expected_vals[i][7])
                assert (marker_build_dict[targets[i]].tree_tool == expected_vals[i][8]) 
                assert (marker_build_dict[targets[i]].num_reps == expected_vals[i][9]) 
                for j in range (0, len(expected_vals[i][10])):
                        assert(expected_vals[i][10][j] in marker_build_dict[targets[i]].pfit)

    def test_parse_ref_build_params(self):
        expected_out = [['McrA', 'M0701', 'prot', 'PROTGAMMALG', '0.97', 'Classes', '04_Dec_2018','functional', 'FastTree', '211', [-4.08046639871, 6.03601100802], '04_Dec_2018'], ['p_amoA', 'N0102', 'prot', 'PROTGAMMALG', '0.97', 'Families', '04_Dec_2018', 'functional', 'FastTree', '80', [-2.83232814805, 5.67790899421], '04_Dec_2018'], ['narG', 'D0101', 'prot', 'PROTGAMMALG', '0.80', 'Phyla', '04_Dec_2018', 'functional', 'FastTree', '307', [-4.5658136261, 6.43765586015], '04_Dec_2018']]

        args = create_parser(TREESAPP_PATH, 'M0701', 'p')
        marker_build_dict = file_parsers.parse_ref_build_params(args)
        self.check_parse_ref_build_params_out(expected_out, ['M0701'], marker_build_dict)

        
        args = create_parser(TREESAPP_PATH, 'ALL', 'p')
        marker_build_dict = file_parsers.parse_ref_build_params(args)
        self.check_parse_ref_build_params_out(expected_out, ['M0701', 'N0102', 'D0101'], marker_build_dict)


    def test_exit(self):
        args = create_parser('~/', 'ALL', 'p')
        with pytest.raises(SystemExit) as pytest_wrapped_e:
                file_parsers.parse_ref_build_params(args)
                assert pytest_wrapped_e.type == SystemExit
                assert pytest_wrapped_e.value.code == 5
                

    def test_parse_cog_list(self):
        expected_out = [['McrA', 'M0701', 'prot', 'PROTGAMMALG', '0.97', 'Classes', 'Methyl coenzyme M reductase alpha subunit','functional', 'FastTree', '211', [-4.08046639871, 6.03601100802], '04_Dec_2018'], ['p_amoA', 'N0102', 'prot', 'PROTGAMMALG', '0.97', 'Families', 'Ammonia monooxygenase (Archaea) and particulate methane monoxygenase', 'functional', 'FastTree', '80', [-2.83232814805, 5.67790899421], '04_Dec_2018'], ['narG', 'D0101', 'prot', 'PROTGAMMALG', '0.80', 'Phyla', 'nitrate reductase / nitrite oxidoreductase, alpha subunit', 'functional', 'FastTree', '307', [-4.5658136261, 6.43765586015], '04_Dec_2018']]

        args = create_parser(TREESAPP_PATH, 'M0701', 'p')
        marker_build_dict = file_parsers.parse_cog_list(args, file_parsers.parse_ref_build_params(args))
        self.check_parse_ref_build_params_out(expected_out, args.targets, marker_build_dict)

        args = create_parser(TREESAPP_PATH, 'ALL', 'p')
        args.targets = ['M0701', 'N0102', 'D0101']
        marker_build_dict = file_parsers.parse_cog_list(args, file_parsers.parse_ref_build_params(args))
        self.check_parse_ref_build_params_out(expected_out, args.targets, marker_build_dict)


    def test_read_species_translation_files(self):
        args = create_parser(TREESAPP_PATH, 'M0701', 'p')
        marker_build_dict = file_parsers.parse_ref_build_params(args)
        marker_build_dict = file_parsers.parse_cog_list(args, marker_build_dict)
        tree_numbers_translation = file_parsers.read_species_translation_files(args, marker_build_dict)
        assert(tree_numbers_translation['M0701'][0].lineage == 'cellular organisms; Archaea')
        assert tree_numbers_translation['M0701'][0].complete
        assert(tree_numbers_translation['M0701'][0].number == '1')
        assert(tree_numbers_translation['M0701'][0].description == "NM4 | NM423249334157")


    
