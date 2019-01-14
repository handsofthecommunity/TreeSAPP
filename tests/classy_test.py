import pytest
import unittest
import os

import sys, inspect
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))))
from classy import ItolJplace, TreeProtein, TreeLeafReference, ReferenceSequence, MarkerBuild, ReferencePackage, TaxonTest, Cluster, Header, register_headers
from jplace_utils import jplace_parser
from json import loads
import file_parsers
import treesapp

HOME_DIR = '/home/travis/build/hallamlab/TreeSAPP/'
TREESAPP_TEST_DIR = '/home/travis/build/hallamlab/TreeSAPP/tests/test_data/'

def create_itol():
    itol = jplace_parser(HOME_DIR + 'tests/test_data/RAxML_portableTree.M0701_hmm_purified_group0-BMGE-qcd.phy.jplace')
    return itol

def init_itol():
    itol = create_itol()
    itol.correct_decoding()
    itol.create_jplace_node_map()
    return itol

class TreeSappTest(unittest.TestCase):

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


def short_fasta():
    fasta_file = TEST_DATA_PATH + '/short_fasta.fa'
    return fasta_file, open(fasta_file, 'r')

results = [('213_McrA', 'M---------------------------------------------------------------------------------------------------------------------------------------------------------------AKKIEKTQKLFLKALKEKFA-------------EDPQS--TSTIFAREGLKQS--PRKMEFVKAGNAAA-MSR--GLSMYDPVRCHI---GGIPLGQRQLMTYEVSGT-G---------------------VFVEGDDLHFVNNAAMQQMWDDIRRTILVNMDLAHQTLQKRLGKEVTPETINEFLHVVNHAMPGA-AVVQEHMVETHPSLVDDCYVKVFTGDDELADDLEPQFVINVEKLFPG------K-QA----AQLKAAVGKSLWQAIRIPTIVSRTCDGGTTSRWSAMQLGMSFIGAYHMCAGEAATADLAYAAKHAGVIQMAE-ILPARRARGPNEPGGIKFGHFADMVQT-DRKYPH-----------------DPAKASLEVV-AAGTMLFDQIWLGSYMSGG-VGFTQ-YATAAYTDNILDDYTYYGMDY-IKDKYKVDWKNPG-EKDKV-KP-TQEVVNDIASE-VTLYGMEQYEQFPTALETHFGGSQRASVLAAASGLSTAIATGNSNAGLNGW-YLSMLLHKEGWSRLGFYGYDLQDQCGSANTESYRADEGCVGELRGANYPNYAMNVGHQGEYAAIAGAAHITRGDAWALNPLIKIAFADP-SLKFDFSEPRREFAKGAIREF-MPAGERALIIP-AR-----------------------'), ('214_McrA','----------------------------------------------------------------------------------------------------------------------------------------------------------------MAKIERTQKLFLKSLKEKFA------------G-DPTG-TTASYFTFGDMKQS--PRKMEFLEQGRRVS-MDR--GISQYDPRRAHL---GGIPLGQRQLMTYEVSTT-G---------------------VFVEGDDLHFVNNSAMQQCWDDIRRTVIVGMDLAHQTLQKRLGKEVTPETINEYLHVLNHAMPGA-AVVQEHMVETAPALVDDCYVKVFSGDDELVDDLEPQFVLNVDKLFPA------K-QA----EGLKAAVGKSLWQAVHIPTTVSRTCDGGTTSRWSAMQLGMSYIAAYRMCAGEAAVADLSFAAKHAGVIQMAS-HLPARRARGPNEPGGIGFGLFSDIIQA-NRKYPN-----------------DPARASLEVV-AAGTMLFDQIWLGSYMSGG-VGFTQ-YATAAYTDNILDEYTYYGMDY-LKDKYKVDWKNPS-PADRV-KA-SQDIVNDLATE-VSLNAMEQYEQFPTLMEDHFGGSQRAGVIAAACGLTCSIGTGNSNAGLNGW-YLSMLLHKEGWSRLGFFGYDLQDQCGSTNSLSIRPDEGAMGEVRGPNYPNYAMNVGHQGEYAAIVGGAHYGRGDGWCFDPRVAITFADP-ALKFDFAEPRREFAKGAIREF-MPAGERSLIIP-AR-----------------------')]

class FastaTests(unittest.TestCase):

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

    def test_get_header_format(self):
        assert(True)

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

            
class MarkerBuildTest(unittest.TestCase):

    def test_create_MarkerBuild(self):
        line = 'McrA\tM0701\tprot\tPROTGAMMALG\tfunctional\t0.97\t211\tFastTree\t-4.08046639871,6.03601100802\tClasses\t04_Dec_2018'
        mb = classy.MarkerBuild()

        assert(mb.cog == '')
        assert(mb.denominator == '')
        assert(mb.molecule == '')
        assert(mb.model == '')
        assert(mb.kind == '')
        assert(mb.pid == 1.0)
        assert(mb.num_reps == 0)
        assert(mb.tree_tool == "")
        assert(mb.lowest_confident_rank == '')
        assert(mb.update == '')
        assert(mb.description == '')
    
    def test_load_build_params(self):
        line = 'McrA\tM0701\tprot\tPROTGAMMALG\tfunctional\t0.97\t211\tFastTree\t-4.08046639871,6.03601100802\tClasses\t04_Dec_2018'
        mb = classy.MarkerBuild()
        mb.load_build_params(line)
        assert(mb.cog == 'McrA')
        assert(mb.denominator == 'M0701')
        assert(mb.molecule == 'prot')
        assert(mb.model == 'PROTGAMMALG')
        assert(mb.kind == 'functional')
        assert(mb.pid == '0.97')
        assert(mb.num_reps == '211')
        assert(mb.tree_tool == "FastTree")
        assert(mb.lowest_confident_rank == 'Classes')
        assert(mb.update == '04_Dec_2018')
        assert(mb.description == '04_Dec_2018')

        incomplete_line = 'McrA\tM0701\tprot\tPROTGAMMALG\tfunctional\t0.97\t211\tFastTree\t-4.08046639871,6.03601100802'
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            mb.load_build_params(incomplete_line)
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 17

    def test_load_pfit_params(self):
        line = 'McrA\tM0701\tprot\tPROTGAMMALG\tfunctional\t0.97\t211\tFastTree\t-4.08046639871,6.03601100802\tClasses\t04_Dec_2018'
        mb = classy.MarkerBuild()
        mb.load_build_params(line)
        mb.load_pfit_params(line)
        assert(-4.08046639871 in mb.pfit and 6.03601100802 in mb.pfit and len(mb.pfit) == 2)
        
    def test_check_rank(self):
        line = 'McrA\tM0701\tprot\tPROTGAMMALG\tfunctional\t0.97\t211\tFastTree\t-4.08046639871,6.03601100802\tClasses\t04_Dec_2018'
        mb = classy.MarkerBuild()
        assert(mb.load_build_params(line) == None)

        
        line_bad = 'McrA\tM0701\tprot\tPROTGAMMALG\tfunctional\t0.97\t211\tFastTree\t-4.08046639871,6.03601100802\tBurgers\t04_Dec_2018'
        mb.load_build_params(line_bad)
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            mb.check_rank()
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 17

class ItolJplaceTest(unittest.TestCase):

    def test_init(self):
        itol = ItolJplace()
        assert(itol.contig_name == "")
        assert(itol.name == "")
        assert(itol.abundance == None)
        assert(itol.node_map == dict())
        assert(itol.seq_len == 0)
        assert(itol.lineage_list == list())
        assert(itol.wtd == 0)
        assert(itol.lct == "")
        assert(itol.placements == list())
        assert(itol.lwr == 0)
        assert(itol.likelihood == 0)
        assert(itol.avg_evo_dist == 0.0)
        assert(itol.distances == '')
        assert(itol.classified)
        assert(itol.inode == "")
        assert(itol.tree == "")
        assert(itol.metadata == "")
        assert(itol.version == "")

    def test_summarize(self):
        itol = ItolJplace()
        assert(itol.summarize() == "\nInformation for query sequence ''\n0 sequence(s) grafted onto the  tree.\nJPlace fields:\n\t[]\nPlacement information:\n\tNone.\nNon-redundant lineages of child nodes:\n\tNone.\nLowest common taxonomy:\n\tNone.\n\n")

    def test_list_placements(self):
        itol = ItolJplace()
        assert(len(itol.list_placements()) == 0)

        itol = create_itol()
        itol.correct_decoding()
        actual_entries = itol.list_placements()
        entries = ['9', '201', '200', '50', '305', '304', '215', '309', '35', '77', '166', '50', '9', '263', '266', '267', '265', '264', '262', '268']
        for entry in entries:
            assert(entry in actual_entries)
            
    def test_correct_decoding(self):
        itol = ItolJplace()
        assert(itol.correct_decoding() == None)

        itol = create_itol()
        
    def test_name_placed_sequence(self):
        itol = ItolJplace()
        assert(itol.name_placed_sequence() == None)

    def test_get_field_position_from_jplace_fields(self):
        itol = ItolJplace()
        field_name = ""

        assert(itol.get_field_position_from_jplace_fields('') == None)
        itol = create_itol()
        itol.correct_decoding()
        assert(itol.get_field_position_from_jplace_fields("edge_num") == 0)
        assert(itol.get_field_position_from_jplace_fields("likelihood") == 1)
        assert(itol.get_field_position_from_jplace_fields("like_weight_ratio") == 2)
        assert(itol.get_field_position_from_jplace_fields("distal_length") == 3)
 
    def test_filter_min_weight_threshold(self):
        itol = create_itol()
        itol.correct_decoding()
        assert(itol.filter_min_weight_threshold() == None)

    def test_sum_rpkms_per_node(self):
        itol = init_itol()
        itol.placements = [itol.placements[0]]
        itol.abundance = 1.0
        leaf_rpkm_sums = dict()
        result = itol.sum_rpkms_per_node(leaf_rpkm_sums)
        assert(result['212'] == 1.75)
        assert(result['188'] == 1.75)
        assert(result['190'] == 1.75)
        assert(result['180'] == 1.75)
        assert(len(result) == 4)

    def test_get_jplace_element(self):
        itol = create_itol()
        itol.correct_decoding()
        assert(itol.get_jplace_element("distal_length") == 0.006343)
        assert(itol.get_jplace_element("pendant_length") == 0.066892)
        assert(itol.get_jplace_element("distal_length") == 0.006343)
        assert(itol.get_jplace_element("distal_length") == 0.006343)
        
    def test_filter_max_weight_placement(self):
        itol = ItolJplace()
        assert(itol.filter_max_weight_placement() == None)
        itol = create_itol()
        itol.correct_decoding()
        itol.filter_max_weight_placement()
        placement = loads(itol.placements[0], encoding='wtf-8')
        assert(placement['n'][0] == 'AFD09581.1_methyl-coenzyme_M_reductase_alpha_subunit__partial_uncultured_Methanomicrobiales_archaeon_1_254')

        assert(263 in placement['p'][0] and 0.003777 in placement['p'][0] and 0.071205 in placement['p'][0] and 0.215867 in placement['p'][0])

    def test_create_jplace_node_map(self):
        itol = ItolJplace()
        assert(itol.create_jplace_node_map() == None)
        itol = create_itol()
        itol.correct_decoding()
        assert(itol.create_jplace_node_map() == None)
        itol.create_jplace_node_map()
        for i in range(425):
            assert(i in itol.node_map.keys())
        
    def test_harmonize_placements(self):
        itol = ItolJplace()
        itol = create_itol()
        itol.name = 'McrA'
        itol.placements = [itol.placements[0]]
        itol.correct_decoding()
        itol.create_jplace_node_map()
        itol.harmonize_placements(HOME_DIR)
        results = loads(itol.placements[0], encoding='utf-8')
        assert(results['n'][0] == 'AFD09581.1_methyl-coenzyme_M_reductase_alpha_subunit__partial_uncultured_Methanomicrobiales_archaeon_1_254')
        assert(results['p'][0] == [261, -41251.840196, 0.97, 0, 0])
        
    def test_clear_object(self):
        itol = ItolJplace()
        assert(itol.clear_object() == None)
        itol = init_itol()
        itol.clear_object()
        assert(len(itol.placements) == 0)
        assert(len(itol.fields) == 0)
        assert(len(itol.node_map) == 0)
        assert(itol.contig_name == '')
        assert(itol.name == '')
        assert(itol.tree == "")
        assert(itol.metadata == "")
        assert(itol.version == "")
        assert(itol.lct == "")
        assert(len(itol.lineage_list) == 0)
        assert(itol.abundance == None)
        
class TreeProteinTest(unittest.TestCase):
    def test_transfer(self):
        itol = create_itol();
        itol.correct_decoding()
        treeprotein = TreeProtein()
        treeprotein.transfer(itol)
        assert(treeprotein.placements == itol.placements)
        assert(treeprotein.fields == itol.fields)        
        assert(treeprotein.version == itol.version)
        assert(treeprotein.metadata == itol.metadata)

    def test_megan_lca(self):
        itol = init_itol()
        treesapp = TreeProtein()
        treesapp.transfer(itol)
        ###########################
        
class TreeLeafReferenceTest(unittest.TestCase):
    def test_init(self):
        tlr = TreeLeafReference(1, 'NM4 | NM423249334157')
        assert(tlr.number == 1)
        assert(tlr.description == 'NM4 | NM423249334157')
        assert(tlr.lineage == "")
        assert(tlr.complete == False)
        tlr2 = TreeLeafReference(2, 'M. Mx2 | 2617404690')
        tlr3 = TreeLeafReference(18, 'Candidatus Methanomethylophilus sp. 1R26 | KUE73676')

class ReferenceSequenceTest(unittest.TestCase):
    def test_init(self):
        rs = ReferenceSequence()
        assert(rs.accession == '')
        assert(rs.description == '')
        assert(rs.organism == '')
        assert(rs.lineage == '')
        assert(rs.short_id == '')
        assert(rs.sequence == '')
        assert(rs.locus == '')
        assert(not rs.cluster_rep)
        assert(rs.cluster_rep_similarity == 0)
        assert(rs.cluster_lca == None)

    def test_get_info(self):
        rs = ReferenceSequence()
        assert(rs.get_info() == 'accession = , mltree_id = \ndescription = , locus = \norganism = \nlineage = \n')
        
class MarkerBuildTest(unittest.TestCase):
        
    def test_create_MarkerBuild(self):
        line = 'McrA\tM0701\tprot\tPROTGAMMALG\tfunctional\t0.97\t211\tFastTree\t-4.08046639871,6.03601100802\tClasses\t04_Dec_2018'
        mb = MarkerBuild()

        assert(mb.cog == '')
        assert(mb.denominator == '')
        assert(mb.molecule == '')
        assert(mb.model == '')
        assert(mb.kind == '')
        assert(mb.pid == 1.0)
        assert(mb.num_reps == 0)
        assert(mb.tree_tool == "")
        assert(mb.lowest_confident_rank == '')
        assert(mb.update == '')
        assert(mb.description == '')
    
    def test_load_build_params(self):
        line = 'McrA\tM0701\tprot\tPROTGAMMALG\tfunctional\t0.97\t211\tFastTree\t-4.08046639871,6.03601100802\tClasses\t04_Dec_2018'
        mb = MarkerBuild()
        mb.load_build_params(line)
        assert(mb.cog == 'McrA')
        assert(mb.denominator == 'M0701')
        assert(mb.molecule == 'prot')
        assert(mb.model == 'PROTGAMMALG')
        assert(mb.kind == 'functional')
        assert(mb.pid == '0.97')
        assert(mb.num_reps == '211')
        assert(mb.tree_tool == "FastTree")
        assert(mb.lowest_confident_rank == 'Classes')
        assert(mb.update == '04_Dec_2018')
        assert(mb.description == '04_Dec_2018')

        incomplete_line = 'McrA\tM0701\tprot\tPROTGAMMALG\tfunctional\t0.97\t211\tFastTree\t-4.08046639871,6.03601100802'
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            mb.load_build_params(incomplete_line)
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 17

    def test_load_pfit_params(self):
        line = 'McrA\tM0701\tprot\tPROTGAMMALG\tfunctional\t0.97\t211\tFastTree\t-4.08046639871,6.03601100802\tClasses\t04_Dec_2018'
        mb = MarkerBuild()
        mb.load_build_params(line)
        mb.load_pfit_params(line)
        assert(-4.08046639871 in mb.pfit and 6.03601100802 in mb.pfit and len(mb.pfit) == 2)
        
    def test_check_rank(self):
        line = 'McrA\tM0701\tprot\tPROTGAMMALG\tfunctional\t0.97\t211\tFastTree\t-4.08046639871,6.03601100802\tClasses\t04_Dec_2018'
        mb = MarkerBuild()
        assert(mb.load_build_params(line) == None)

        
        line_bad = 'McrA\tM0701\tprot\tPROTGAMMALG\tfunctional\t0.97\t211\tFastTree\t-4.08046639871,6.03601100802\tBurgers\t04_Dec_2018'
        mb.load_build_params(line_bad)
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            mb.check_rank()
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 17

            
class MarkerInfoTest(unittest.TestCase):
    def test_mi_init(self):
            mi = TreeLeafReference.MarkerInfo('marker_ph', 'R0016', 'test')
            assert(mi.marker == 'marker_ph')
            assert(mi.denominator == 'R0016')
            assert(mi.marker_class == '')
            assert(mi.description == 'test')
            assert(mi.analysis_type == "")
                                       
class ReferencePackageTest(unittest.TestCase):

    def test_rp_init(self):
        rp = ReferencePackage()
        assert(rp.prefix == "")
        assert(rp.msa == "")
        assert(rp.profile == "")
        assert(rp.tree == "")
        assert(rp.boot_tree == "")
        assert(rp.lineage_ids == "")
        assert(rp.num_seqs == 0)
        assert(len(rp.core_ref_files) == 0)
        
class TaxonTestTest(unittest.TestCase):

    def test_tt_init(self):
        tt = TaxonTest('; name')
        assert(tt.lineage == '; name')
        assert(tt.taxon == 'name')
        assert(len(tt.queries) == 0)
        assert(len(tt.classifieds) == 0)
        assert(len(tt.distances) == 0)
        assert(len(tt.assignments) == 0)
        assert(tt.taxonomic_tree == None)

class ClusterTest(unittest.TestCase):
    def test_cluster_init(self):
        cluster = Cluster('rep')
        assert(cluster.representative == 'rep')
        assert(len(cluster.members) == 0)
        assert(cluster.lca == '')
        
class HeaderTest(unittest.TestCase):
    def test_header_init(self):
        header = Header('>1_McrA')
        assert(header.original ==  '>1_McrA')
        assert(header.formatted == '')
        assert(header.treesapp_name == '')
        assert(header.post_align == '')
        assert(header.first_split == '')

    def test_get_info(self):
        header = Header('>2_McrA')
        result_info_str = header.get_info()
        
    def test_register_headers(self):
        ref_headers = ['>1_McrA', '>2_McrA', '>3_McrA', '>4_McrA', '>5_McrA']
        actual_hr = register_headers(ref_headers)
        for i in range(1,6):
            assert(str(i) in actual_hr.keys())
            val = actual_hr[str(i)]
            assert(val.original == '>' + str(i) + '_McrA')
            assert(val.formatted == '>' + str(i) + '_McrA')
            assert(val.first_split == '>' + str(i) + '_McrA')

    #MarkerTest, MyFormatter, TaxonTest
