from __future__ import print_function,division
import os, sys

print_msg_intro = """
#########################################
      WELCOME TO CHIMERA TESTS
THIS SCRIPT RUNS A FEW TYPICAL SCENARIOS
#########################################
"""

print_msg = """
#########################################
TEST {:d}: {}
#########################################
"""

print_msg_end ="""
#########################################
{:d} TESTS SUCCEEDED
{:d} TESTS FAILED
#########################################
"""


print (print_msg_intro)

tests_succ = 0
tests_fail = 0

print(print_msg.format(1,'FREE ELECTRON LASER AMPLIFICATION'))
msg=os.system('python ./fel-testrun.py')
if msg==0:
	tests_succ +=1
else:
	tests_fail +=1

print(print_msg.format(2,'LASER PLASMA ACCELERATION'))
msg = os.system('python ./lpa-testrun.py')
if msg==0:
	tests_succ +=1
else:
	tests_fail +=1
print(print_msg_end.format(tests_succ,tests_fail))

sys.exit(tests_fail)
