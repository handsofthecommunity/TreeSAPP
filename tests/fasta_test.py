import pytest
import unittest
import os

from .treesapp_test import create_parser

import sys, inspect
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))))
import fasta

TREESAPP_PATH = '/home/travis/build/hallamlab/TreeSAPP/'
TEST_DATA_PATH = '/home/travis/build/hallamlab/TreeSAPP/tests/test_data/'

def short_fasta():
    fasta_file = TEST_DATA_PATH + '/short_fasta.fa'
    return fasta_file, open(fasta_file, 'r')

def get_formatted_fasta_dict():
    args = create_parser(TREESAPP_PATH, 'M0701', 'p')
    args.fasta_input = 'tests/test_data/short_fasta_valid.faa'
    formatted_fasta_dict = fasta.format_read_fasta(args.fasta_input, "prot", args.output)
    return args, formatted_fasta_dict

results = [('213_McrA', 'M---------------------------------------------------------------------------------------------------------------------------------------------------------------AKKIEKTQKLFLKALKEKFA-------------EDPQS--TSTIFAREGLKQS--PRKMEFVKAGNAAA-MSR--GLSMYDPVRCHI---GGIPLGQRQLMTYEVSGT-G---------------------VFVEGDDLHFVNNAAMQQMWDDIRRTILVNMDLAHQTLQKRLGKEVTPETINEFLHVVNHAMPGA-AVVQEHMVETHPSLVDDCYVKVFTGDDELADDLEPQFVINVEKLFPG------K-QA----AQLKAAVGKSLWQAIRIPTIVSRTCDGGTTSRWSAMQLGMSFIGAYHMCAGEAATADLAYAAKHAGVIQMAE-ILPARRARGPNEPGGIKFGHFADMVQT-DRKYPH-----------------DPAKASLEVV-AAGTMLFDQIWLGSYMSGG-VGFTQ-YATAAYTDNILDDYTYYGMDY-IKDKYKVDWKNPG-EKDKV-KP-TQEVVNDIASE-VTLYGMEQYEQFPTALETHFGGSQRASVLAAASGLSTAIATGNSNAGLNGW-YLSMLLHKEGWSRLGFYGYDLQDQCGSANTESYRADEGCVGELRGANYPNYAMNVGHQGEYAAIAGAAHITRGDAWALNPLIKIAFADP-SLKFDFSEPRREFAKGAIREF-MPAGERALIIP-AR-----------------------'), ('214_McrA','----------------------------------------------------------------------------------------------------------------------------------------------------------------MAKIERTQKLFLKSLKEKFA------------G-DPTG-TTASYFTFGDMKQS--PRKMEFLEQGRRVS-MDR--GISQYDPRRAHL---GGIPLGQRQLMTYEVSTT-G---------------------VFVEGDDLHFVNNSAMQQCWDDIRRTVIVGMDLAHQTLQKRLGKEVTPETINEYLHVLNHAMPGA-AVVQEHMVETAPALVDDCYVKVFSGDDELVDDLEPQFVLNVDKLFPA------K-QA----EGLKAAVGKSLWQAVHIPTTVSRTCDGGTTSRWSAMQLGMSYIAAYRMCAGEAAVADLSFAAKHAGVIQMAS-HLPARRARGPNEPGGIGFGLFSDIIQA-NRKYPN-----------------DPARASLEVV-AAGTMLFDQIWLGSYMSGG-VGFTQ-YATAAYTDNILDEYTYYGMDY-LKDKYKVDWKNPS-PADRV-KA-SQDIVNDLATE-VSLNAMEQYEQFPTLMEDHFGGSQRAGVIAAACGLTCSIGTGNSNAGLNGW-YLSMLLHKEGWSRLGFFGYDLQDQCGSTNSLSIRPDEGAMGEVRGPNYPNYAMNVGHQGEYAAIVGGAHYGRGDGWCFDPRVAITFADP-ALKFDFAEPRREFAKGAIREF-MPAGERSLIIP-AR-----------------------')]


class FastaTest(unittest.TestCase):

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

    def test_write_new_fasta(self):
        args, formatted_fasta_dict = get_formatted_fasta_dict()
        fasta_name = TEST_DATA_PATH + '/test_new_fasta.fasta'

        fasta.write_new_fasta(formatted_fasta_dict, fasta_name)
        args.fasta_input = 'tests/test_data/test_new_fasta.fasta'

        formatted_fasta_dict = fasta.format_read_fasta(args.fasta_input, "prot", args.output)
        assert('>k127_1003429_914638_1_#_2_#_1513_#_1_#_ID=914638_1_partial=10_start_type=Edge_rbs_motif=None_rbs_spacer=None' in formatted_fasta_dict.keys())
        assert('>k127_35937_flag_381292_3_#_288_#_416_#_-1_#_ID=381292_3_partial=01_start_type=Edge_rbs_motif=None_rbs_spacer' in formatted_fasta_dict.keys())
        assert( '>Prodigal_Seq_6_6_3_#_3683_#_4678_#_-1_#_ID=6_3_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_nosZ' in formatted_fasta_dict.keys())


    # def test_deduplicate(self):
    #     fasta_dict = fasta.read_fasta_to_dict(TEST_DATA_PATH + "/dup_fasta.fa")
    #     fasta_dict = fasta.deduplicate_fasta_sequences(fasta_dict)

    #     for i in range(0, len(results)):
    #         assert(results[i][0] in fasta_dict.keys())
    #         assert(results[i][1] in fasta_dict.values())

    #     fasta_dict = fasta.read_fasta_to_dict(TEST_DATA_PATH + "/short_fasta.fa")
    #     fasta_dict = fasta.deduplicate_fasta_sequences(fasta_dict)

    #     for i in range(0, len(results)):
    #         assert(results[i][0] in fasta_dict.keys())
    #         assert(results[i][1] in fasta_dict.values())
