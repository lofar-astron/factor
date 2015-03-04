import os
from taskinit import *

def ftw(vis=None,field=None,spw=None,model=None,nterms=None,reffreq=None,wprojplanes=None,complist=None,incremental=None, usescratch=None):
       """ Insert a source model into the MODEL_DATA column of a visibility set:

       A source model (souce.model image) or components list is converted into a
       model visibility that is inserted into the MODEL_DATA column.
       TASK created by R.J. van Weeren (Dec. 2013, adapted from ft-task)
       If you do not need Wprojection than simply use the task ft.


       Keyword arguments:
       vis -- Name of input visibility file
               default: none; example: vis='ngc5921.ms'
       field -- Field name list
               default: '' ==> all
               NOTE: each source must be specified in a multi-source vis.
               field = '1328+307'  specifies source '1328+307'
               field = '4' specified field with index 4
       spw -- Spw selection
               default: spw = '' (all spw)
       model -- Name of input model image
               default: None;
               example: model='/usr/lib/casapy/data/nrao/VLA/CalModels/3C286_X.im'
               Note: The model visibilities are scaled from the model frequency
                     to the observed frequency of the data.
       nterms -- Number of terms used to model the sky frequency dependence
                 default: 1
                 example : nterms=3  represents a 2nd order Taylor-polynomial in frequency
                           and is to be used along with 3 model-image names.
		           model=['xxx.image.tt0','xxx.image.tt1', 'xxx.image.tt2']
          reffreq -- Reference-frequency about which this Taylor-expansion is defined.
       wprojplanes -- the number of pre-computed w-planes used for
                   the W-Projection algorithm.
       complist -- Name of component list
               default: None; ; example: complist='test.cl'
               components tool not yet available
       incremental -- Add model visibility to the existing MODEL_DATA visibilties
               default: False; example: incremental=True

       """
       casalog.origin('ftw')

       #Python script
       try:
               # Check if datafile exists and open it
               if ((type(vis)==str) & (os.path.exists(vis))):
                       im.open(vis, usescratch=usescratch)
               else:
                       raise Exception, 'Visibility data set not found - please verify the name'

               # Select data
               im.selectvis(field=field,spw=spw)

               # Define image co-ordinates (all defaults)
               im.defineimage()

               # Check 'model'. The 'xml' allows a variant => do the checking here.
               if( (not type(model)==str) and (not (type(model)==list) ) ) :
		       raise Exception, 'The model image must be a string or a list of strings (or \'\' or [])';

               # If model is a single string, make it a list
               if( type(model)==str ):
                       model = [model];

               # Check that either a model or a complist has been given.
               if( (model==[] or model==['']) and complist=='' ):
                       raise Exception, 'Please specify a model image or component list to ft';

               #model is a list now. Check that all elements are strings. If so, check file existence too.
               if( type(model)==list ):
                       for onemodel in model:
                              if(not type(onemodel)==str):
                                    raise Exception, 'Model image names must be strings';
                              if( (not onemodel=='') and (not os.path.exists(onemodel)) ):
                                    raise Exception, 'Model image '+onemodel+' cannot be found';

               # Check complist : one string : name of complist file. Check existance on disk.
               if( (not complist=='') and (not os.path.exists(complist)) ):
                       raise Exception, 'Componentlist '+complist+' cannot be found';


               # If nterms>1, then check that len(model)=nterms [ no multifield for now ]
               # Call im.settaylorterms()
               #
               if (nterms > 1) :
		       if(type(model)==str or (not (type(model)==list and len(model)==nterms)) ):
			       raise Exception, 'For nterms>1, please provide a list of nterms model-image names';
		       # parse the reference-frequency field.
                       qat=qatool();
                       try:
		          rff=qat.canonical(reffreq);
		       except Exception, instance:
                          print '*** Error *** In conversion of reffreq=\'',reffreq,'\' to a numerical value';
                          raise Exception, instance
                       reffreqVal=rff['value'];  # This is the frequency in Hz
		       if(reffreqVal==0.0):   # if unspecified, set the default from the model image
			       ia.open(model[0]);
			       icsys = ia.coordsys();
			       ia.close();
                               reffreqVal=icsys.referencevalue(type='spectral')['numeric'][0];
			       casalog.post('Using reference frequency from model image : '+str(reffreqVal)+' Hz');
		       else:
		               casalog.post('Using reference frequency : '+str(reffreqVal)+' Hz');
		       # set nterms and ref-freq
		       im.settaylorterms(ntaylorterms=nterms,reffreq=reffreqVal)

               # Just checking...
	       if (nterms < 1) :
		       raise Exception, 'nterms must be greater than or equal to 1';
	       if (wprojplanes < 1) :
		       raise Exception, 'wprojplanes must be greater than or equal to 1';

	       # Use Wprojection

	       im.setoptions(ftmachine='wproject',wprojplanes=wprojplanes, padding=1.2)

               # Do the forward transform and close.
               im.ft(model=model,complist=complist,incremental=incremental)
               im.close()


               #write history
               ms.open(vis,nomodify=False)
               ms.writehistory(message='taskname = ftw',origin='ftw')
               ms.writehistory(message='vis         = "'+str(vis)+'"',origin='ftw')
               ms.writehistory(message='field       = "'+str(field)+'"',origin='ftw')
               ms.writehistory(message='spw         = "'+str(spw)+'"',origin='ftw')
               ms.writehistory(message='model       = "'+str(model)+'"',origin='ftw')
               ms.writehistory(message='complist    = "'+str(complist)+'"',origin='ftw')
               ms.writehistory(message='incremental = "'+str(incremental)+'"',origin='ftw')
               ms.close()

       except Exception, instance:
               print '*** Error ***',instance
               raise Exception, instance
