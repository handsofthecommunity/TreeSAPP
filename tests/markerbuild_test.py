import pytest
import unittest
import os

import sys, inspect
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))))
import classy

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
