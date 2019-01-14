import pytest
import unittest
import os
import argparse

import sys, inspect
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))))

import file_parsers
import fasta
import treesapp

# def create_parser(treesapp, targets, reftree):
#     args = argparse.Namespace()
#     args.alignment_mode = 'd'
#     args.reftree = reftree
#     args.targets = [targets]
#     args.treesapp = treesapp
#     args.check_trees = False
#     args.skip = 'n'
#     args.molecule = 'prot'
#     return args

def create_parser(treesapp, targets, reftree):
    args = argparse.Namespace()
    args.alignment_mode = 'd'
    args.reftree = reftree
    args.targets = [targets]
    args.treesapp = treesapp
    args.check_trees = False
    args.fasta_input = 'test_data/marker_test_suite.faa'
    args.output = '/home/travis/build/hallamlab/marker_test/'
    args.output_dir_var = '/home/travis/build/hallamlab/marker_test/various_outputs'
    args.skip = 'n'
    args.molecule = 'prot'
    args.executables = {'BMGE.jar': '/home/travis/build/hallamlab/TreeSAPP/sub_binaries/BMGE.jar', 'hmmalign': '/usr/bin/hmmalign', 'usearch': HOME_DIR + 'sub_binaries/usearch', 'hmmsearch': '/usr/bin/hmmsearch', 'trimal': '/usr/bin/trimal', 'raxmlHPC': '/usr/bin/raxmlHPC', 'hmmbuild': '/usr/bin/hmmbuild', 'prodigal': '/usr/local/bin/prodigal', 'papara': '/usr/bin/papara'}
    args.reference_data_prefix=''
    args.num_threads = 3
    args.output_dir_final = '/home/travis/build/hallamlab/marker_test/final_outputs/'
    args.formatted_input_file = ''
    args.composition = 'meta'
    args.overwrite = True
    args.delete = False
    args.reclassify = False
    args.min_e = 0.0001
    args.min_acc = 0.7
    args.min_likelihood = 0.2
    args.min_seq_length = 30
    args.perc_aligned = 15
    return args

TEST_PATH = '/home/travis/build/hallamlab/TreeSAPP/tests'
TREESAPP_PATH = '/home/travis/build/hallamlab/TreeSAPP/'

HOME_DIR = '/home/travis/build/hallamlab/TreeSAPP/'
TREESAPP_TEST_DIR = '/home/travis/build/hallamlab/TreeSAPP/tests/test_data/'

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


def create_itol():
    itol = jplace_parser(HOME_DIR + 'tests/test_data/RAxML_portableTree.M0701_hmm_purified_group0-BMGE-qcd.phy.jplace')
    return itol

def init_itol():
    itol = create_itol()
    itol.correct_decoding()
    itol.create_jplace_node_map()
    return itol

class TreeSappTest(unittest.TestCase):

    def create_short_qc_ma_dict():
        qc_ma_dict = {'M0701': {'tests/test_data/expected_phy_file': {'186': '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------KHAGVIQMADILPARRARGPNEPGGIKFGHFGDMIQADPVKATLEVVGAGAMLFDQIWLGSYMSGGYATAAYTDNILDDYCYYGLDYTQEVLNDIATEVTLYGMEQYEQYPTTLESHFGGSQRASVLAAASGISCSLATANSNAGLNGWYMSMLAHKEGWSRLGFFGYDLQDQCGSTNSMSIRPDEGCIGELRGPNYPNYAMNVGHQGEYAAIASAAHYGRQDAWVLSPLIKIAFADPSLKFDFSEPRREFARGAIREFMPAGERSLIIP', '70': 'DDLHFVNNAAIQQMVDDIKRTVIVGMDTAHAVLEKRLGVEVTPETINEYMEAINHALPGGAVVQEHMVEVHPGLVEDCYAKIFTGDDNLADELDKRILIDINKEFPEEQLKSYIGNRTYQVNRVPTIVVRTCDGGTVSRWSAMQIGMSFISAYKLCAGEAAIADFSYAAKHADVIEMGTIMPARRARGPNEPGGVAFGTFADIVQTDPANVSLEVIAGAAALYDQVWLGSYMSGGYATAAYTDDILDDFVYYGMEYTMDVVRDISTEVTLYSLEQYEEYPTLLEDHFGGSQRAAVAAAAAGCSTAFATGNSNAGINGWYLSQILHKEAHSRLGFYGYDLQDQCGASNSLSIRSDEGLIHELRGPNYPNYAMNVGHQPEYAGIAQAPHAARGDAFCTNPLIKVAFADKDLAFDFTSPRKSIAAGALREFMPEGERDLIIP', '59': 'DDLHYVNNAAIQQAWDDIRRTVIVGLNTAHNVLEKRLGIEVTPETITHYLETVNHAMPGAAVVQEHMVETDPLIVQDSYVKVFTGDDELADEIDSAFVLDINKEFPEEALKAEVGGAIWQAVRIPSIVGRVCDGGNTSRWSAMQIGMSMISAYNQCAGEGATGDFAYASKHAEVIHMGTYLPVRRARAENELGGVPFGFMADICQGDPVRVSLEVVALGAALYDQIWLGSYMSGGYATAAYTDNVLDDFTYYGKDYNMDTVLDVGTEVAFYALEQYEEYPALLETHFGGSQRASVVSAAAGCSTAFATGNAQTGLSAWYLAMYLHKEQHSRLGFYGFDLQDQCGAANVFSIRNDEGLPLEMRGPNYPNYAMNVGHQGEYAGIAQAPHAARGDAWAFNPLVKIAFADKNLCFDFSKVREEFAKGALREFEPAGERTAITP'}}}
        return qc_ma_dict

    def arguments():
        args = create_parser(HOME_DIR, 'M0701', 'p')
        args.formatted_input_file = args.output_dir_var + 'marker_test_suite.faa'  + "_formatted.fasta"
        fasta.write_new_fasta(formatted_fasta_dict, args.formatted_input_file)
        marker_build_dict = treesapp.parse_ref_build_params(args)
        marker_build_dict = treesapp.parse_cog_list(args, marker_build_dict)
        formatted_fasta_dict = fasta.format_read_fasta(args.fasta_input, "prot", args.output)
        homolog_seq_files, numeric_contig_index = treesapp.extract_hmm_matches(args, hmm_matches, formatted_fasta_dict)
        return args, marker_build_dict, formatted_fasta_dict, homolog_seq_files, numeric_contig_index


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
        marker_build_dict = treesapp.parse_ref_build_params(args)
        marker_build_dict = treesapp.parse_cog_list(args, marker_build_dict)
        formatted_fasta_dict = fasta.format_read_fasta(args.fasta_input, "prot", args.output)

        fasta.write_new_fasta(formatted_fasta_dict, args.formatted_input_file)

        hmm_domtbl_files = treesapp.hmmsearch_orfs(args, marker_build_dict)
        hmm_matches = treesapp.parse_domain_tables(args, hmm_domtbl_files)
        assert(len(hmm_matches['McrA']) == 12)
        homolog_seq_files, numeric_contig_index = treesapp.extract_hmm_matches(args, hmm_matches, formatted_fasta_dict)
        
        assert(homolog_seq_files == ['/home/travis/build/hallamlab/marker_test/various_outputs/McrA_hmm_purified_group0.faa'])
        assert('McrA' in numeric_contig_index.keys())
        assert(len(numeric_contig_index['McrA']) == 12)

        expected_out = [-12, -2, -10, -9, -8, -7, -6, -5, -4, -3, -1, -11];

        for i in range (0, len(expected_out)):
            assert(expected_out[i] in numeric_contig_index['McrA'].keys())


        #---------------------------------------------------------------------------------
        args = create_parser('~', 'ALL', 'p')
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            treesapp.hmmsearch_orfs(args, marker_build_dict)
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 3
        


    def test_create_ref_phy_files(self):
        args = create_parser(HOME_DIR, 'M0701', 'p')
        args.formatted_input_file = args.output_dir_var + 'marker_test_suite.faa'  + "_formatted.fasta"
        marker_build_dict = treesapp.parse_ref_build_params(args)
        marker_build_dict = treesapp.parse_cog_list(args, marker_build_dict)
        
        formatted_fasta_dict = fasta.format_read_fasta(args.fasta_input, "prot", args.output)
        hmm_domtbl_files = treesapp.hmmsearch_orfs(args, marker_build_dict)
        hmm_matches = treesapp.parse_domain_tables(args, hmm_domtbl_files) 
        fasta.write_new_fasta(formatted_fasta_dict, args.formatted_input_file)
        homolog_seq_files, numeric_contig_index = treesapp.extract_hmm_matches(args, hmm_matches, formatted_fasta_dict)
        ref_alignment_dimensions = treesapp.get_alignment_dims(args, marker_build_dict)
        treesapp.create_ref_phy_files(args, homolog_seq_files, marker_build_dict, ref_alignment_dimensions)
        marker = re.match("(.*)_hmm_purified.*", os.path.basename('/home/travis/build/hallamlab/marker_test/various_outputs/McrA_hmm_purified_group0.faa')).group(1)
        ref_alignment_phy = args.output_dir_var + marker + ".phy"

        assert(filecmp.cmp(ref_alignment_phy, HOME_DIR + 'tests/test_data/expected_ref_alignment.phy'))


    def test_multiple_alignments(self):
        args = create_parser(HOME_DIR, 'M0701', 'p')
        args.formatted_input_file = args.output_dir_var + 'marker_test_suite.faa'  + "_formatted.fasta"
        single_query_sequence_files = '../../marker_test/various_outputs/McrA_hmm_purified_group0.faa'
        marker_build_dict = treesapp.parse_ref_build_params(args)
        marker_build_dict = treesapp.parse_cog_list(args, marker_build_dict)
        tool = 'hmmalign'

        assert(treesapp.multiple_alignments(single_query_sequence_files, marker_build_dict, tool).key == 'M0701')
        
        #invalid tool name
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            treesapp.multiple_alignments(single_query_sequence_files, marker_build_dict, 'none')
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 3
    
        

    def test_prepare_and_run_papara(self):
        assert(True)

    def test_prepare_and_run_hmmalign(self):
        args = create_parser(HOME_DIR, 'M0701', 'p')
        single_query_fasta_files = ['/marker_test/various_outputs/McrA_hmm_purified_group0.faa']
        marker_build_dict = treesapp.parse_ref_build_params(args)
        marker_build_dict = treesapp.parse_cog_list(args, marker_build_dict)

        # result = dict()
        # result = treesapp.prepare_and_run_hmmalign(args, single_query_fasta_files, marker_build_dict)
        # assert('M0701' in result.keys and len(result.keys) == 1)
        
        #wrong file name
        single_query_fasta_files = ['']
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            treesapp.prepare_and_run_hmmalign(args, single_query_fasta_files, marker_build_dict)
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 3

        single_query_fasta_files = ['./McrB_hmm_purified_group0.faa']
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            treesapp.prepare_and_run_hmmalign(args, single_query_fasta_files, marker_build_dict)
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 3


    def test_filter_multiple_alignments(self):
        args = create_parser(HOME_DIR, 'M0701', 'p')

        marker_build_dict = treesapp.parse_ref_build_params(args)
        marker_build_dict = treesapp.parse_cog_list(args, marker_build_dict)

        if re.match(r'\A.*\/(.*)', args.fasta_input):
            input_multi_fasta = os.path.basename(args.fasta_input)
        else:
            input_multi_fasta = args.fasta_input
            
        args.formatted_input_file = args.output_dir_var + input_multi_fasta  + "_formatted.fasta"
        formatted_fasta_dict = fasta.format_read_fasta(args.fasta_input, "prot", args.output)

        ref_alignment_dimensions = treesapp.get_alignment_dims(args, marker_build_dict)
        hmm_domtbl_files = treesapp.hmmsearch_orfs(args, marker_build_dict)
        
        hmm_matches = treesapp.parse_domain_tables(args, hmm_domtbl_files)
        homolog_seq_files, numeric_contig_index = treesapp.extract_hmm_matches(args, hmm_matches, formatted_fasta_dict)

        concatenated_mfa_files = treesapp.multiple_alignments(args, homolog_seq_files, marker_build_dict, "hmmalign")

        mfa_files = treesapp.filter_multiple_alignments(args, concatenated_mfa_files, marker_build_dict, 'BMGE')


        ## test get sequence counts
        assert(treesapp.get_sequence_counts(concatenated_mfa_files, ref_alignment_dimensions, False, 'Fasta') == {'/home/travis/build/hallamlab/marker_test/various_outputs/McrA_hmm_purified_group0.mfa' : 844})

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            treesapp.get_sequence_counts(concatenated_mfa_files, ref_alignment_dimensions, False, 'NotFasta')
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 3


        #test filter multiple alignments
        assert( 'M0701' in treesapp.filter_multiple_alignments(args, concatenated_mfa_files, marker_build_dict, 'BMGE').keys())
        assert(['/home/travis/build/hallamlab/marker_test/various_outputs/McrA_hmm_purified_group0-BMGE.fasta'] in treesapp.filter_multiple_alignments(args, concatenated_mfa_files, marker_build_dict, 'BMGE').values())


        invalid_mfa_file = dict()
        invalid_mfa_file.update({'M0701' : ['/home/travis/build/hallamlab/marker_test/various_outputs/McrA_hmm_purified.txt']})
        
        #check for removed sequences
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            treesapp.check_for_removed_sequences(args, invalid_mfa_file, marker_build_dict)
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 3

        empty_mfa_file = dict()
        empty_mfa_file.update({'M0701' : [HOME_DIR + 'tests/test_data/empty.fasta']})
        
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            treesapp.check_for_removed_sequences(args, empty_mfa_file, marker_build_dict)
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 3


    def test_produce_phy_files(self):
        qc_ma_dict = create_short_qc_ma_dict()

        args = argparse.Namespace()
        args.molecule = 'prot'
        phy_file = treesapp.produce_phy_files(args, qc_ma_dict)

        assert('M0701' in phy_file.keys())
        assert([ HOME_DIR + 'tests/test_data/expected_phy_file'] in phy_file.values() and len(phy_file) == 1)
        with open(phy_file['M0701'][0]) as f:
            content = f.readlines()

        results = ['3', '439', '59', 'DDLHYVNNAAIQQAWDDIRRTVIVGLNTAHNVLEKRLGIEVTPETITHYL', '70',
                   'DDLHFVNNAAIQQMVDDIKRTVIVGMDTAHAVLEKRLGVEVTPETINEYM', '186', 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX', 'ETVNHAMPGAAVVQEHMVETDPLIVQDSYVKVFTGDDELADEIDSAFVLDINKEFPEEALKAEVGGAIWQAVRIPSIVGRVCDGGNTSRWSAMQIGMSMISAYNQCAGEGATGDFAYASKHAEVIHMGTYLPVRRARAENELGGVPFGFMADICQGDPVRVSLEVVALGAALYDQIWLGSYMSGGYATAAYTDNVLDDFTYYGKDYNMDTVLDVGTEVAFYALEQYEEYPALLETHFGGSQRASVVSAAAGCSTAFATGNAQTGLSAWYLAMYLHKEQHSRLGFYGFDLQDQCGAANVFSIRNDEGLPLEMRGPNYPNYAMNVGHQGEYAGIAQAPHAARGDAWAFNPLVKIAFADKNLCFDFSKVREEFAKGALREFEPAGERTAITP', 'EAINHALPGGAVVQEHMVEVHPGLVEDCYAKIFTGDDNLADELDKRILIDINKEFPEEQLKSYIGNRTYQVNRVPTIVVRTCDGGTVSRWSAMQIGMSFISAYKLCAGEAAIADFSYAAKHADVIEMGTIMPARRARGPNEPGGVAFGTFADIVQTDPANVSLEVIAGAAALYDQVWLGSYMSGGYATAAYTDDILDDFVYYGMEYTMDVVRDISTEVTLYSLEQYEEYPTLLEDHFGGSQRAAVAAAAAGCSTAFATGNSNAGINGWYLSQILHKEAHSRLGFYGYDLQDQCGASNSLSIRSDEGLIHELRGPNYPNYAMNVGHQPEYAGIAQAPHAARGDAFCTNPLIKVAFADKDLAFDFTSPRKSIAAGALREFMPEGERDLIIP', 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXKHAGVIQMADILPARRARGPNEPGGIKFGHFGDMIQADPVKATLEVVGAGAMLFDQIWLGSYMSGGYATAAYTDNILDDYCYYGLDYTQEVLNDIATEVTLYGMEQYEQYPTTLESHFGGSQRASVLAAASGISCSLATANSNAGLNGWYMSMLAHKEGWSRLGFFGYDLQDQCGSTNSMSIRPDEGCIGELRGPNYPNYAMNVGHQGEYAAIASAAHYGRQDAWVLSPLIKIAFADPSLKFDFSEPRREFARGAIREFMPAGERSLIIP']

        content = [x.strip() for x in content]
        for l in range(4):
            tmp = content[l].split()
            assert(results[l*2] == tmp[0])
            assert(results[l*2 + 1] == tmp[1])

        for y in range(8):
            if (y == 7):
                assert(content[5 + y*4] == results[8][50*y:])
                assert(content[5 + y*4 + 1] == results[9][50*y:])
                assert(content[5 + y*4 + 2] == results[10][50*y:])
            else :
                assert(content[5 + y*4] == results[8][50*y:(y+1)*50])
                assert(content[5 + y*4 + 1] == results[9][50*y:(y+1)*50])
                assert(content[5 + y*4 + 2] == results[10][50*y:(y+1)*50])


    def test_multiple_alignments(self):
      single_query_sequence_files = ['tests/test_data/McrA_hmm_purified_group0.faa']

      args = create_parser(HOME_DIR, 'M0701', 'p')

      marker_build_dict = treesapp.parse_ref_build_params(args)
      marker_build_dict = treesapp.parse_cog_list(args, marker_build_dict)

      assert(treesapp.multiple_alignments(args, single_query_sequence_files, marker_build_dict, "hmmalign") == '')
      
    def test_sub_indices_for_seq_names_jplace(self):
        short_numeric_contig_index = {'McrA': {-12: 'PHP46140.1_methyl-coenzyme_M_reductase_subunit_alpha_Methanosarcinales_archaeon_ex4572_44_8_595', -2: 'OYT62528.1_hypothetical_protein_B6U67_04395_Methanosarcinales_archaeon_ex4484_138_1_471', -10: 'AAU83782.1_methyl_coenzyme_M_reductase_subunit_alpha_uncultured_archaeon_GZfos33H6_1_570'}}

        args = create_parser(HOME_DIR, 'M0701', 'p')
        marker_build_dict = treesapp.parse_ref_build_params(args)
        marker_build_dict = treesapp.parse_cog_list(args, marker_build_dict)

    def test_validate_inputs(self):
        args = create_parser(HOME_DIR, 'M0701', 'p')

        marker_build_dict = treesapp.parse_ref_build_params(args)
        marker_build_dict = treesapp.parse_cog_list(args, marker_build_dict)

        # Correctly formatted reference tree
        # assert(validate_inputs(args, marker_build_dict))

        # # Incorrectly formatted reference tree
        # new_marker_build = MarkerBuild()
        # marker_build_dict['TEST'] = new_marker_build
        # copyfile(TEST_DIR + '/TEST_tree.txt', args.treesapp + "/data/tree_data/")
        
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            treesapp.validate_inputs(args, marker_build_dict)
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 3  
