import pytest
import unittest
import os

import sys, inspect
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))))
from classy import ItolJplace
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

        itol = create_itol()
        itol.correct_decoding()
        assert(itol.summarize() == '\nInformation for query sequence \'\'\n12 sequence(s) grafted onto the  tree.\nJPlace fields:\n\t[\'"edge_num"\', \'"likelihood"\', \'"like_weight_ratio"\', \'"distal_length"\', \'"pendant_length"\']\nPlacement information:\n\t[[9, -41260.14731, 1.0, 0.220188, 1e-06]]\n\t[[201, -41233.206425, 0.9706, 0.009032, 0.018712], [200, -41237.61184, 0.011852, 0.002155, 0.022476]]\n\t[[50, -41180.626427, 0.999769, 0.020085, 1e-06]]\n\t[[305, -41202.097508, 0.948026, 0.005416, 0.006097], [304, -41205.121452, 0.046083, 0.003828, 0.006218]]\n\t[[215, -41173.620848, 1.0, 0.068377, 1e-06]]\n\t[[309, -41255.874824, 1.0, 0.116764, 0.014456]]\n\t[[35, -41173.620869, 1.0, 0.138396, 1e-06]]\n\t[[77, -41173.620868, 1.0, 0.018389, 1e-06]]\n\t[[166, -41320.629566, 0.999834, 0.025274, 0.05397]]\n\t[[50, -41180.626427, 0.999769, 0.020085, 1e-06]]\n\t[[9, -41173.620668, 1.0, 0.223728, 1e-06]]\n\t[[263, -41251.840196, 0.215867, 0.003777, 0.071205], [266, -41252.081796, 0.169535, 1e-06, 0.069548], [267, -41252.102745, 0.166021, 0.011828, 0.069443], [265, -41252.112437, 0.164419, 1e-06, 0.069452], [264, -41252.664146, 0.0947, 0.003525, 0.07157], [262, -41252.670307, 0.094118, 1e-06, 0.073156], [268, -41252.974363, 0.069442, 0.006343, 0.066892]]\nNon-redundant lineages of child nodes:\n\tNone.\nLowest common taxonomy:\n\tNone.\n\n')
        
    def test_list_placements(self):
        itol = ItolJplace()
        assert(len(itol.list_placements()) == 0)

        itol = create_itol()
        itol.correct_decoding()
        assert(itol.list_placements() == ['9', '201', '200', '50', '305', '304', '215', '309', '35', '77', '166', '50', '9', '263', '266', '267', '265', '264', '262', '268'])
        
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
        itol = ItolJplace()
        element_value = None
        
    def test_filter_max_weight_placement(self):
        itol = ItolJplace()
        assert(itol.filter_max_weight_placement() == None)
        
    def test_create_jplace_node_map(self):
        itol = ItolJplace()
        assert(itol.create_jplace_node_map() == None)

    def test_harmonize_placements(self):
        itol = ItolJplace()

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
