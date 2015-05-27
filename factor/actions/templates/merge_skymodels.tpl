#!/usr/bin/env python
"""
Script to merge two sky models
"""
import os
import lsmtool
import sys

if __name__ == '__main__':
    import optparse
    opt = optparse.OptionParser(usage='%prog <model1_filename> <model2_filename> <outmodel_filename>')
    (options, args) = opt.parse_args()

    if len(args) != 3:
        opt.print_help()
        sys.exit(1)

    inmodel1 = args[0]
    inmodel2 = args[1]
    outmodel = args[2]

    s1 = lsmtool.load(inmodel1)
    s2 = lsmtool.load(inmodel2)

    s1.concatenate(s2, matchBy='{{ matchby }}', radius={{ radius }}, keep='{{ keep }}',
        inheritPatches=True)
    s1.group('every')
    s1.write(fileName=outmodel, clobber=True)
