from __future__ import print_function
import sys
import unittest
import traceback
import os

# Append the path of the module to the syspath
sys.path.append(os.path.join('..','..','..'))
from factor.actions.imager import imager

class TestImager(unittest.TestCase):
    
    def setUp(self):
        self.im = imager("op1", "test_ms")
    
    def test_command(self):
        print(self.im._get_command())
    
if __name__ == '__main__':
    unittest.main()