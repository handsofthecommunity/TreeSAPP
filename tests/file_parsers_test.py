import pytest
import argparse
import file_parsers
import unittest
import os

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
            assert (len(marker_build_dict[targets[i]].pfit) == 0)
            assert(marker_build_dict[targets[i]].description == expected_vals[i][7])
            assert (marker_build_dict[targets[i]].kind == expected_vals[i][8])


    def test_parse_ref_build_params(self):
        expected_out = [['McrA', 'M0701', 'prot', 'PROTGAMMALG', 0.97, 'Classes', '01_Aug_2018', '', ''],
                        ['hzs', 'A0100', 'prot', 'PROTGAMMAJTTDCMUT', 0.99, 'Species', '27_Apr_2018', '', ''],
                        ['narG', 'D0101', 'prot', 'PROTGAMMALG', 80, 'Genera', '15_Jan_2018', '', '']]

        args = create_parser('~/github/TreeSAPP', 'M0701', 'p')
        marker_build_dict = file_parsers.parse_ref_build_params(args)
        self.check_parse_ref_build_params_out(expected_out, ['M0701'], marker_build_dict)

        args = create_parser('~/github/TreeSAPP', 'ALL', 'p')
        marker_build_dict = file_parsers.parse_ref_build_params(args)
        self.check_parse_ref_build_params_out(expected_out, ['M0701', 'A0100', 'D0101'], marker_build_dict)


    def test_exit(self):
        args = create_parser('~/', 'ALL', 'p')
        with pytest.raises(exit(5)):
            file_parsers.parse_ref_build_params(args)

        marker_build_dict = file_parsers.parse_ref_build_params(args)
        with pytest.raises(exit(5)):
            file_parsers.parse_cog_list(args, marker_build_dict)

        args = create_parser('~/', 'ALL', 'z')
        with pytest.raises(exit(5)):
            file_parsers.parse_cog_list(args, marker_build_dict)


    def test_parse_cog_list(self):
        expected_out = [['McrA', 'M0701', 'prot', 'PROTGAMMALG', 0.97, 'Classes', '01_Aug_2018', 'Methyl coenzyme M reductase alpha subunit', 'functional_cogs'],
                        ['hzs', 'A0100', 'prot', 'PROTGAMMAJTTDCMUT', 0.99, 'Species', '27_Apr_2018', 'Hydrazine synthase (hzs)', 'functional_cogs'],
                        ['narG', 'D0101', 'prot', 'PROTGAMMALG', 80, 'Genera', '15_Jan_2018', 'Methyl coenzyme M reductase alpha subunit', 'functional_cogs']]

        args = create_parser('~/github/TreeSAPP', 'M0701', 'p')
        marker_build_dict = file_parsers.parse_cog_list(args, file_parsers.parse_ref_build_params(args))
        self.check_parse_ref_build_params_out(expected_out, args.targets, marker_build_dict)

        args = create_parser('~/github/TreeSAPP', 'ALL', 'p')
        marker_build_dict = file_parsers.parse_cog_list(args, file_parsers.parse_ref_build_params(args))
        self.check_parse_ref_build_params_out(expected_out, args.targets, marker_build_dict)


    def test_read_species_translation_files(self):
        args = create_parser('~/github/TreeSAPP', 'M0701', 'p')
        marker_build_dict = file_parsers.parse_ref_build_params(args)
        marker_build_dict = file_parsers.parse_cog_list(args, marker_build_dict)
        tree_numbers_translation = file_parsers.read_species_translation_files(args, marker_build_dict)
        assert(tree_numbers_translation['M0701'] == file_parsers.tax_ids_file_to_leaves(os.sep.join([args.treesapp, "data", "tree_data"]) + os.sep + 'tax_ids_' + str(marker_build_dict['M0701'].cog) + '.txt'))
        assert(tree_numbers_translation['M0701'][0].lineage == 'cellular organisms; Archaea')
        assert tree_numbers_translation['M0701'][0].complete
        assert(tree_numbers_translation['M0701'][0].number == 1)
        assert(tree_numbers_translation['M0701'][0].translation == "NM4 | NM423249334157")

    def test_tax_ids_file_to_leaves(self):
        args = create_parser('~/github/TreeSAPP', 'M0701', 'p')
        marker_build_dict = file_parsers.parse_ref_build_params(args)
        marker_build_dict = file_parsers.parse_cog_list(args, marker_build_dict)
        tree_numbers_translation = file_parsers.read_species_translation_files(args, marker_build_dict)
        tree_leaves = file_parsers.tax_ids_file_to_leaves(tree_numbers_translation)
        assert(tree_leaves[0].lineage == 'cellular organisms; Archaea')
        assert tree_leaves[0].complete
        assert(tree_leaves[0].number == 1)
        assert(tree_leaves[0].translation == "NM4 | NM423249334157")




