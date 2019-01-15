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
    args.skip = 'n'
    args.molecule = 'prot'
    return args

TEST_PATH = '/home/travis/build/hallamlab/TreeSAPP/tests'
TREESAPP_PATH = '/home/travis/build/hallamlab/TreeSAPP/'

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
    
    def test_read_colours_file(self):
        marker_subgroups = dict()
        internal_nodes = dict()
        annotation_file = TREESAPP_PATH + 'data/iTOL_datasets/McrA_colours_style.txt'
        args = create_parser(TREESAPP_PATH, 'M0701', 'p')

        marker_subgroups['McrA'] = dict()
        internal_nodes['McrA'] = dict()
        marker_subgroups['McrA'], internal_nodes['McrA'] = file_parsers.read_colours_file(args, annotation_file)

        expected_keys = ['Bathyarchaeota', 'Methanosarcinales', 'Verstraetearchaeota', 'Methanomicrobiales', 'Methanococcales', 'Methanocellales', 'Methanomassiliicoccales', 'Methanococcales-MrtA', 'ANME-1', 'Syntrophoarchaeum', 'Methanopyrus', 'Methanobacteriales', 'Methanonatronarchaeia', 'Methanobacteriales-MrtA']

        for key in expected_keys:
            assert(key in marker_subgroups['McrA'].keys())

        assert(marker_subgroups['McrA']['Bathyarchaeota'][0] == ('204', '203'))
        assert(not internal_nodes['McrA'])

    def test_tax_ids_file_to_leaves(self):
        marker_tax_ids = TREESAPP_PATH + 'data/tree_data/tax_ids_McrA.txt'
        ref_tree_leaves = file_parsers.tax_ids_file_to_leaves(marker_tax_ids)

        assert(len(ref_tree_leaves) == 214)
        assert(ref_tree_leaves[0].description == 'NM4 | NM423249334157')
        assert(ref_tree_leaves[0].lineage == 'cellular organisms; Archaea')

        for i in range(214):
            assert(ref_tree_leaves[i].number == str(i + 1))
            assert(ref_tree_leaves[i].complete)

        marker_tax_ids = TEST_PATH + '/test_data/empty.fasta'

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            ref_tree_leaves = file_parsers.tax_ids_file_to_leaves(marker_tax_ids)
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 5

        marker_tax_ids = TEST_PATH + '/test_data/tax_ids_McrA.txt' 
        with pytest.raises(ValueError) as excinfo:
            ref_tree_leaves = file_parsers.tax_ids_file_to_leaves(marker_tax_ids)
            assert ecxinfo.value.code == 5

TREESAPP_PATH = '/home/travis/build/hallamlab/TreeSAPP/'
TEST_DATA_PATH = '/home/travis/build/hallamlab/TreeSAPP/tests/test_data/'

class FastaTest(unittest.TestCase):

    def short_fasta():
        fasta_file = TEST_DATA_PATH + '/short_fasta.fa'
        return fasta_file, open(fasta_file, 'r')

    def test_read_fasta_to_dict(self):
        fasta_file, fasta_handler = short_fasta()
        assert(fasta.read_fasta_to_dict(fasta_file) == {results[0][0]: results[0][1], results[1][0]:results[1][1]})

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            fasta.read_fasta_to_dict(' ')
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 5
                                                            
    def test_generate_fasta(self):
        fasta_file, fasta_handler = short_fasta()
        i = 0;
        for record in fasta.generate_fasta(fasta_handler):
            assert(record == results[i])               
            i+=1
    

    def test_format_read_fasta(self):
        args = create_parser(TREESAPP_PATH, 'M0701', 'p')
        args.fasta_input = 'tests/test_data/short_fasta_valid.faa'
        formatted_fasta_dict = fasta.format_read_fasta(args.fasta_input, "prot", args.output)
        assert('>k127_1003429_914638_1_#_2_#_1513_#_1_#_ID=914638_1_partial=10_start_type=Edge_rbs_motif=None_rbs_spacer=None' in formatted_fasta_dict.keys())
        assert('>k127_35937_flag_381292_3_#_288_#_416_#_-1_#_ID=381292_3_partial=01_start_type=Edge_rbs_motif=None_rbs_spacer' in formatted_fasta_dict.keys())
        assert( '>Prodigal_Seq_6_6_3_#_3683_#_4678_#_-1_#_ID=6_3_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_nosZ' in formatted_fasta_dict.keys())

        args.fasta_input = 'tests/test_data/fasta_invalid.faa'

        with pytest.raises(SystemExit) as pytest_wrapped_e:
             fasta.format_read_fasta(args.fasta_input, "prot", args.output)
             assert pytest_wrapped_e.type == SystemExit
             assert pytest_wrapped_e.value.code == 5

    def test_get_headers(self):
       ref_headers = fasta.get_headers(TEST_DATA_PATH + '/short_fasta.fa')
       assert(ref_headers == ['>213_McrA', '>214_McrA'])

       with pytest.raises(SystemExit) as pytest_wrapped_e:
             fasta.get_headers('')
             assert pytest_wrapped_e.type == SystemExit
             assert pytest_wrapped_e.value.code == 5

    def test_generate_fasta(self):
        assert(True)

    def test_write_new_fasta(self):
        args, formatted_fasta_dict = get_formatted_fasta_dict()
        fasta_name = TEST_DATA_PATH + '/test_new_fasta.fasta'

        fasta.write_new_fasta(formatted_fasta_dict, fasta_name)
        args.fasta_input = 'tests/test_data/test_new_fasta.fasta'
        
        formatted_fasta_dict = fasta.format_read_fasta(args.fasta_input, "prot", args.output)
        assert('>k127_1003429_914638_1_#_2_#_1513_#_1_#_ID=914638_1_partial=10_start_type=Edge_rbs_motif=None_rbs_spacer=None' in formatted_fasta_dict.keys())
        assert('>k127_35937_flag_381292_3_#_288_#_416_#_-1_#_ID=381292_3_partial=01_start_type=Edge_rbs_motif=None_rbs_spacer' in formatted_fasta_dict.keys())
        assert( '>Prodigal_Seq_6_6_3_#_3683_#_4678_#_-1_#_ID=6_3_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_nosZ' in formatted_fasta_dict.keys())

    def test_deduplicate(self):
        fasta_dict = fasta.read_fasta_to_dict(TEST_DATA_PATH + "/dup_fasta.fa")
        fasta_dict = fasta.deduplicate_fasta_sequences(fasta_dict)

        for i in range(0, len(results)):
            assert(results[i][0] in fasta_dict.keys())
            assert(results[i][1] in fasta_dict.values())

        
        fasta_dict = fasta.read_fasta_to_dict(TEST_DATA_PATH + "/short_fasta.fa")
        fasta_dict = fasta.deduplicate_fasta_sequences(fasta_dict)

        for i in range(0, len(results)):
            assert(results[i][0] in fasta_dict.keys())
            assert(results[i][1] in fasta_dict.values())
