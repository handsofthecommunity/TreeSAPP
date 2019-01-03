import pytest
import unittest
import os

import sys, inspect
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))))
from classy import ItolJplace, TreeProtein, TreeLeafReference, ReferenceSequence, MarkerBuild, ReferencePackage, TaxonTest, Cluster, Header, register_headers
from jplace_utils import jplace_parser

HOME_DIR = '/home/travis/build/hallamlab/TreeSAPP/'

def create_itol():
    itol = jplace_parser(HOME_DIR + 'tests/test_data/RAxML_portableTree.M0701_hmm_purified_group0-BMGE-qcd.phy.jplace')
    return itol

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
        
    def test_create_jplace_node_map(self):
        itol = ItolJplace()
        assert(itol.create_jplace_node_map() == None)
        filename = '/home/ace/github/TreeSAPP/tests/test_data/RAxML_portableTree.M0701_hmm_purified_group0-BMGE-qcd.phy.jplace'
        jplace_data = jplace_parser(filename)
        
    def test_harmonize_placements(self):
        itol = ItolJplace()
        treesapp_dir = '/home/ace/github/TreeSAPP/'

    def test_clear_object(self):
        itol = ItolJplace()
        assert(itol.clear_object() == None)
        

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
