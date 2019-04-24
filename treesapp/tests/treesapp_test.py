import pytest
import unittest
import os
import argparse
import re

import sys, inspect
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))))
import file_parsers
import fasta
import treesapp

HOME_DIR = '/home/travis/build/hallamlab/TreeSAPP/'
TREESAPP_TEST_DIR = '/home/travis/build/hallamlab/TreeSAPP/tests/test_data/'
TEST_DIR = '/home/travis/build/hallamlab/TreeSAPP/tests/test_data/'

def create_short_qc_ma_dict():
    qc_ma_dict = {'M0701': {'tests/test_data/expected_phy_file': {'186': '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------KHAGVIQMADILPARRARGPNEPGGIKFGHFGDMIQADPVKATLEVVGAGAMLFDQIWLGSYMSGGYATAAYTDNILDDYCYYGLDYTQEVLNDIATEVTLYGMEQYEQYPTTLESHFGGSQRASVLAAASGISCSLATANSNAGLNGWYMSMLAHKEGWSRLGFFGYDLQDQCGSTNSMSIRPDEGCIGELRGPNYPNYAMNVGHQGEYAAIASAAHYGRQDAWVLSPLIKIAFADPSLKFDFSEPRREFARGAIREFMPAGERSLIIP', '70': 'DDLHFVNNAAIQQMVDDIKRTVIVGMDTAHAVLEKRLGVEVTPETINEYMEAINHALPGGAVVQEHMVEVHPGLVEDCYAKIFTGDDNLADELDKRILIDINKEFPEEQLKSYIGNRTYQVNRVPTIVVRTCDGGTVSRWSAMQIGMSFISAYKLCAGEAAIADFSYAAKHADVIEMGTIMPARRARGPNEPGGVAFGTFADIVQTDPANVSLEVIAGAAALYDQVWLGSYMSGGYATAAYTDDILDDFVYYGMEYTMDVVRDISTEVTLYSLEQYEEYPTLLEDHFGGSQRAAVAAAAAGCSTAFATGNSNAGINGWYLSQILHKEAHSRLGFYGYDLQDQCGASNSLSIRSDEGLIHELRGPNYPNYAMNVGHQPEYAGIAQAPHAARGDAFCTNPLIKVAFADKDLAFDFTSPRKSIAAGALREFMPEGERDLIIP', '59': 'DDLHYVNNAAIQQAWDDIRRTVIVGLNTAHNVLEKRLGIEVTPETITHYLETVNHAMPGAAVVQEHMVETDPLIVQDSYVKVFTGDDELADEIDSAFVLDINKEFPEEALKAEVGGAIWQAVRIPSIVGRVCDGGNTSRWSAMQIGMSMISAYNQCAGEGATGDFAYASKHAEVIHMGTYLPVRRARAENELGGVPFGFMADICQGDPVRVSLEVVALGAALYDQIWLGSYMSGGYATAAYTDNVLDDFTYYGKDYNMDTVLDVGTEVAFYALEQYEEYPALLETHFGGSQRASVVSAAAGCSTAFATGNAQTGLSAWYLAMYLHKEQHSRLGFYGFDLQDQCGAANVFSIRNDEGLPLEMRGPNYPNYAMNVGHQGEYAGIAQAPHAARGDAWAFNPLVKIAFADKNLCFDFSKVREEFAKGALREFEPAGERTAITP'}}}
    return qc_ma_dict

def arguments():
    args = create_parser(HOME_DIR, 'M0701', 'p')
    args.formatted_input_file = args.output_dir_var + 'marker_test_suite.faa'  + "_formatted.fasta"
    fasta.write_new_fasta(formatted_fasta_dict, args.formatted_input_file)
    marker_build_dict = treesapp.parse_ref_build_params(args)
    formatted_fasta_dict = fasta.format_read_fasta(args.fasta_input, "prot", args.output)
    homolog_seq_files, numeric_contig_index = treesapp.extract_hmm_matches(args, hmm_matches, formatted_fasta_dict)
    return args, marker_build_dict, formatted_fasta_dict, homolog_seq_files, numeric_contig_index
    

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
    args.min_score = 15
    args.min_acc = 0.7
    args.min_ie = 0.1
    args.min_likelihood = 0.2
    args.min_seq_length = 30
    args.perc_aligned = 15
    return args

class TempTest(unittest.TestCase):

    def test_prepare_and_run_hmmalign(self):
        args = create_parser(HOME_DIR, 'M0701', 'p')
        single_query_fasta_files = ['/marker_test/various_outputs/McrA_hmm_purified_group0.faa']
        marker_build_dict = treesapp.parse_ref_build_params(args)

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
        assert(treesapp.get_sequence_counts(concatenated_mfa_files, ref_alignment_dimensions, False, 'Fasta') == {'/home/travis/build/hallamlab/marker_test/various_outputs/McrA_hmm_purified_group0.mfa' : 866})

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
        empty_mfa_file.update({'M0701' : [TEST_DIR + 'empty.fasta']})
            
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            treesapp.check_for_removed_sequences(args, empty_mfa_file, marker_build_dict)
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 3

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
        assert(['tests/test_data/expected_phy_file'] in phy_file.values() and len(phy_file) == 1)
        with open(phy_file['M0701'][0]) as f:
            content = f.readlines()

        results = ['3  439', '186       XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX', '59        DDLHYVNNAAIQQAWDDIRRTVIVGLNTAHNVLEKRLGIEVTPETITHYL','70        DDLHFVNNAAIQQMVDDIKRTVIVGMDTAHAVLEKRLGVEVTPETINEYM', '', 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX', 'ETVNHAMPGAAVVQEHMVETDPLIVQDSYVKVFTGDDELADEIDSAFVLD', 'EAINHALPGGAVVQEHMVEVHPGLVEDCYAKIFTGDDNLADELDKRILID', '', 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX', 'INKEFPEEALKAEVGGAIWQAVRIPSIVGRVCDGGNTSRWSAMQIGMSMI', 'INKEFPEEQLKSYIGNRTYQVNRVPTIVVRTCDGGTVSRWSAMQIGMSFI', '', 'XXXXXXXXXXXXXXXXXXXKHAGVIQMADILPARRARGPNEPGGIKFGHF', 'SAYNQCAGEGATGDFAYASKHAEVIHMGTYLPVRRARAENELGGVPFGFM', 'SAYKLCAGEAAIADFSYAAKHADVIEMGTIMPARRARGPNEPGGVAFGTF', '', 'GDMIQADPVKATLEVVGAGAMLFDQIWLGSYMSGGYATAAYTDNILDDYC', 'ADICQGDPVRVSLEVVALGAALYDQIWLGSYMSGGYATAAYTDNVLDDFT', 'ADIVQTDPANVSLEVIAGAAALYDQVWLGSYMSGGYATAAYTDDILDDFV', '', 'YYGLDYTQEVLNDIATEVTLYGMEQYEQYPTTLESHFGGSQRASVLAAAS', 'YYGKDYNMDTVLDVGTEVAFYALEQYEEYPALLETHFGGSQRASVVSAAA', 'YYGMEYTMDVVRDISTEVTLYSLEQYEEYPTLLEDHFGGSQRAAVAAAAA', '', 'GISCSLATANSNAGLNGWYMSMLAHKEGWSRLGFFGYDLQDQCGSTNSMS', 'GCSTAFATGNAQTGLSAWYLAMYLHKEQHSRLGFYGFDLQDQCGAANVFS', 'GCSTAFATGNSNAGINGWYLSQILHKEAHSRLGFYGYDLQDQCGASNSLS', '', 'IRPDEGCIGELRGPNYPNYAMNVGHQGEYAAIASAAHYGRQDAWVLSPLI', 'IRNDEGLPLEMRGPNYPNYAMNVGHQGEYAGIAQAPHAARGDAWAFNPLV']

        content = [x.strip() for x in content]

        for i in range(30):
            assert(content[i] == results[i])
            
     
    def test_multiple_alignments(self):
      single_query_sequence_files = [TREESAPP_TEST_DIR + 'McrA_hmm_purified_group0.faa']

      args = create_parser(HOME_DIR, 'M0701', 'p')

      marker_build_dict = treesapp.parse_ref_build_params(args)

      assert(treesapp.multiple_alignments(args, single_query_sequence_files, marker_build_dict, "hmmalign")['M0701'] == [TREESAPP_TEST_DIR + 'McrA_hmm_purified_group0.mfa'])
      
        

    def test_sub_indices_for_seq_names_jplace(self):
        short_numeric_contig_index = {'McrA': {-12: 'PHP46140.1_methyl-coenzyme_M_reductase_subunit_alpha_Methanosarcinales_archaeon_ex4572_44_8_595', -2: 'OYT62528.1_hypothetical_protein_B6U67_04395_Methanosarcinales_archaeon_ex4484_138_1_471', -10: 'AAU83782.1_methyl_coenzyme_M_reductase_subunit_alpha_uncultured_archaeon_GZfos33H6_1_570'}}

        args = create_parser(HOME_DIR, 'M0701', 'p')
        marker_build_dict = treesapp.parse_ref_build_params(args)
        
    def test_validate_inputs(self):
        args = create_parser(HOME_DIR, 'M0701', 'p')

        marker_build_dict = treesapp.parse_ref_build_params(args)

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

    def test_get_alignment_data(self):
        args = create_parser(HOME_DIR, 'M0701', 'p')
        marker_build_dict = file_parsers.parse_ref_build_params(args)
        alignment_dimensions_dict = treesapp.get_alignment_dims(args, marker_build_dict)
        assert(alignment_dimensions_dict['M0701'] == (211, 851))

           

class TreeSAPPTest(unittest.TestCase):

    def test_hmmsearch_orfs_parse_domain_tables(self):
        args = create_parser(HOME_DIR, 'M0701', 'p')
        marker_build_dict = file_parsers.parse_ref_build_params(args)
        hmm_domtbl_files = treesapp.hmmsearch_orfs(args, marker_build_dict)
        assert(hmm_domtbl_files[0] == '/home/travis/build/hallamlab/marker_test/various_outputs/McRA_to_ORFs_domtbl.txt')

        args = create_parser('~', 'ALL', 'p')
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            treesapp.hmmsearch_orfs(args, marker_build_dict)
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 3

    def test_extract_hmm_matches(self):
        args = create_parser(HOME_DIR, 'M0701', 'p')
        args.formatted_input_file = args.output_dir_var + args.fasta_input + "_formatted.fasta"
        marker_build_dict = treesapp.parse_ref_build_params(args)
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
        
