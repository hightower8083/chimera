from __future__ import print_function,division
import os 

tests_succ = 0
tests_fail = 0

print(\
"""
#########################################
WELCOME TO CHIMERA TESTS!
THIS SCRIPT RUNS A FEW TYPICAL SCENARIOS
#########################################
TEST 1: FREE ELECTRON LASER AMPLIFICATION
#########################################
""")

msg=os.system('python ./doc/fel-testrun.py')
if msg==0:
	tests_succ +=1
else:
	tests_fail +=1

print(\
"""
#########################################
TEST 2: LASER PLASMA ACCELERATION
#########################################
""")

msg = os.system('python ./doc/lpa-testrun.py')
if msg==0:
	tests_succ +=1
else:
	tests_fail +=1

print(\
"""
#########################################
{:d} TESTS SUCCEEDED
{:d} TESTS FAILED
#########################################
""".format(tests_succ,tests_fail))
