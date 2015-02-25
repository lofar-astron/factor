import lsmtool
import sys

inmodel1 = sys.argv[1]
inmodel2 = sys.argv[2]
outmodel = sys.argv[3]

s1 = lsmtool.load(inmodel1)
s2 = lsmtool.load(inmodel2)

s1.merge(s2, matchBy={{ matchby }}, radius={{ radius }}, keep={{ keep }},
    inheritPatches=True)
s1.group('every')
s1.write(outfile=outmodel, clobber=True)
