import pytest
import unittest
import os

import sys, inspect
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))))
import classy

class MarkerBuildTest(unittest.TestCase):
        
    def test_create_MarkerBuild(self):
        line = 'hzs\tA0100\tprot\tPROTGAMMAJTTDCMUT\t0.99\t\tSpecies\t27_Apr_2018'
        mb = classy.MarkerBuild(line)
        assert(mb.cog == 'hzs')
        assert(mb.denominator == 'A0100')
        assert(mb.molecule == 'prot')
        assert(mb.model == 'PROTGAMMAJTTDCMUT')
        assert(mb.pid == '0.99')
        assert(mb.lowest_confident_rank == 'Species')
        assert(mb.update == '27_Apr_2018')
        assert(len(mb.pfit) == 0)

        line_bad = 'hzs\tA0100\tprot\tPROTGAMMAJTTDCMUT\t0.99\t\tSpecies'
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            classy.MarkerBuild(line_bad)
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 17
            
    def test_load_pfit_params(self):

        line_true = 'nifD\tN0401\tprot\tPROTGAMMALG\t85\t-6.9676601767537747,3954157334089663\tSpecies\t17_Jan_2018'
        mb_true = classy.MarkerBuild(line_true)
        mb_true.load_pfit_params(line_true)
        assert(mb_true.pfit == [-6.9676601767537747, 3954157334089663])

        line_false = 'nifD\tN0401\tprot\tPROTGAMMALG\t85\t\tSpecies\t17_Jan_2018'
        mb_false = classy.MarkerBuild(line_false)
        mb_false.load_pfit_params(line_false)
        assert(len(mb_false.pfit) == 0)

    # #FUNCTION NOT USED
    # def test_load_rank_distances(self):
    #     line_break = 'hzs\tA0100\tprot\tPROTGAMMAJTTDCMUT\t0.99\t\tSpecies\t27_Apr_2018'
    #     mb_break = classy.MarkerBuild(line_break)

    #     assert(mb_break.load_rank_distances(line_break))
    #     assert(len(mb_break.distances) == 0)

    #     line = 'nifD\tN0401\tprot\tPROTGAMMALG\t85\t-6.9676601767537747,3954157334089663\t\tSpecies\t17_Jan_2018'
    #     mb = classy.MarkerBuild(line)
    #     assert(not mb.load_rank_distances(line))

    def test_check_rank(self):
        line = 'hzs\tA0100\tprot\tPROTGAMMAJTTDCMUT\t0.99\t\tNotSpecies\t27_Apr_2018'
        mb = classy.MarkerBuild(line)
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            mb.check_rank()
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 17

        line = 'hzs\tA0100\tprot\tPROTGAMMAJTTDCMUT\t0.99\t\tSpecies\t27_Apr_2018'
        mb = classy.MarkerBuild(line)
        assert mb.check_rank() is None
        
