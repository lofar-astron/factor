import os
import sys

ms = sys.argv[6]
imageout = sys.argv[7]

mask = {{ mask }}
scales={{ scales }}
cell=['{{ cell }}', '{{ cell }}']
imsize=[{{ imsize }}, {{ imsize }}]
threshold='{{ threshold }}'
uvrange='{{ uvrange }}'
wplanes={{ wplanes }}
nterms={{ nterms }}
niter={{ niter }}
timer = '{{ timer }}'
nfacets = {{ nfacets }}
cycfactor = {{ cycfactor }}

clean(vis=ms,imagename=imageout,outlierfile="",field="",spw="",selectdata=True,timerange=timer,\
      uvrange=uvrange,antenna="",scan="",observation="",mode="mfs",gridmode="widefield",wprojplanes=wplanes,\
      facets=nfacets,cfcache="cfcache.dir",painc=360.0,epjtable="",interpolation="linear",        \
      niter=niter,gain=0.1,threshold=threshold,psfmode="clark",imagermode="csclean",        \
      ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=scales,negcomponent=-1,   \
      smallscalebias=0.6,interactive=False,mask=mask,nchan=-1,start=0,width=1,outframe="",  \
      veltype="radio",imsize=imsize,cell=cell,phasecenter="",restfreq="",stokes="I",        \
      weighting="briggs",robust=-0.25,uvtaper=False,outertaper=[''],innertaper=['1.0'],     \
      modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=True,noise="1.0Jy",\
      npixels=0,npercycle=100,cyclefactor=cycfactor,cyclespeedup=-1,nterms=nterms,reffreq="",          \
      chaniter=False,flatnoise=True,allowchunk=False)

os.system('touch {{ completed_file }}')
