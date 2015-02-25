import lsmtool
import sys

inmodel1 = sys.argv[1]
inmodel2 = sys.argv[2]
outputdir = sys.argv[3]

root1 = os.path.basename(inmodel1).split('model')[1]
root2 = os.path.basename(inmodel2).split('model')[1]
outskymodel = os.path.splitext(root1)[0] + '_merged.skymodel'
outmodel = os.path.join(outputdir, outskymodel)

s1 = lsmtool.load(inmodel1)
s2 = lsmtool.load(inmodel2)

s1.merge(s2, matchBy={{ matchby }}, radius={{ radius }}, keep={{ keep }},
    inheritPatches=True)
s1.group('every')
s1.write(outfile=outmodel)
