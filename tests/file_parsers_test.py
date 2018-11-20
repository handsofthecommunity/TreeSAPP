import pytest
import unittest
import os
import argparse

import sys, inspect
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))))
import file_parsers

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
                assert (marker_build_dict[targets[i]].update == expected_vals[i][6])
                assert (marker_build_dict[targets[i]].description == expected_vals[i][7])
                assert (marker_build_dict[targets[i]].kind == expected_vals[i][8])


    def test_parse_ref_build_params(self):
        expected_out = [['McrA', 'M0701', 'prot', 'PROTGAMMALG', '0.97', 'Classes', '01_Aug_2018', '', ''], ['hzs', 'A0100', 'prot', 'PROTGAMMAJTTDCMUT', '0.99', 'Species', '27_Apr_2018', '', ''], ['narG', 'D0101', 'prot', 'PROTGAMMALG', '80', 'Genera', '23_Jan_2018', '', '']]

        args = create_parser('/home/travis/build/hallamlab/TreeSAPP/', 'M0701', 'p')
        marker_build_dict = file_parsers.parse_ref_build_params(args)
        self.check_parse_ref_build_params_out(expected_out, ['M0701'], marker_build_dict)

        
        args = create_parser('/home/travis/build/hallamlab/TreeSAPP/', 'ALL', 'p')
        marker_build_dict = file_parsers.parse_ref_build_params(args)
        self.check_parse_ref_build_params_out(expected_out, ['M0701', 'A0100', 'D0101'], marker_build_dict)


    def test_exit(self):
        args = create_parser('~/', 'ALL', 'p')
        with pytest.raises(SystemExit) as pytest_wrapped_e:
                file_parsers.parse_ref_build_params(args)
                assert pytest_wrapped_e.type == SystemExit
                assert pytest_wrapped_e.value.code == 5
                

    def test_parse_cog_list(self):
        expected_out = [['McrA', 'M0701', 'prot', 'PROTGAMMALG', '0.97', 'Classes', '01_Aug_2018', 'Methyl coenzyme M reductase alpha subunit', 'functional_cogs'], ['hzs', 'A0100', 'prot', 'PROTGAMMAJTTDCMUT', '0.99', 'Species', '27_Apr_2018', 'Hydrazine synthase (hzs)', 'functional_cogs'], ['narG', 'D0101', 'prot', 'PROTGAMMALG', '80', 'Genera', '23_Jan_2018', 'nitrate reductase / nitrite oxidoreductase, alpha subunit', 'functional_cogs']]

        args = create_parser('/home/travis/build/hallamlab/TreeSAPP/', 'M0701', 'p')
        marker_build_dict = file_parsers.parse_cog_list(args, file_parsers.parse_ref_build_params(args))
        self.check_parse_ref_build_params_out(expected_out, args.targets, marker_build_dict)

        args = create_parser('/home/travis/build/hallamlab/TreeSAPP/', 'ALL', 'p')
        args.targets = ['M0701', 'A0100', 'D0101']
        marker_build_dict = file_parsers.parse_cog_list(args, file_parsers.parse_ref_build_params(args))
        self.check_parse_ref_build_params_out(expected_out, args.targets, marker_build_dict)


    def test_read_species_translation_files(self):
        args = create_parser('/home/travis/build/hallamlab/TreeSAPP/', 'M0701', 'p')
        marker_build_dict = file_parsers.parse_ref_build_params(args)
        marker_build_dict = file_parsers.parse_cog_list(args, marker_build_dict)
        tree_numbers_translation = file_parsers.read_species_translation_files(args, marker_build_dict)
        assert(tree_numbers_translation['M0701'][0].lineage == 'cellular organisms; Archaea')
        assert tree_numbers_translation['M0701'][0].complete
        assert(tree_numbers_translation['M0701'][0].number == '1')
        assert(tree_numbers_translation['M0701'][0].description == "NM4 | NM423249334157")


