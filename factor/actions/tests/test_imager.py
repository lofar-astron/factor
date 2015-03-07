from __future__ import print_function
import sys
import unittest
import traceback
import os

# Append the path of the module to the syspath
sys.path.append(os.path.join('..','..','..'))
from factor.actions.images import MakeImage

class TestImager(unittest.TestCase):

    def setUp(self):
        self.im = MakeImage({"op_name":"op1"}, "test_ms")


if __name__ == '__main__':
    unittest.main()
