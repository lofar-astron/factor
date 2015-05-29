<?xml version="1.0" encoding="UTF-8"?>
<?xml-stylesheet type="text/xsl" ?>
<casaxml xmlns="http://casa.nrao.edu/schema/psetTypes.html"
      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
      xsi:schemaLocation="http://casa.nrao.edu/schema/casa.xsd
      file:///opt/casa/code/xmlcasa/xml/casa.xsd">

<task type="function" name="ftw" category="imaging, calibration">
	<shortdescription>Insert a source model  a visibility set:</shortdescription>
	<description>

       TASK created by R.J. van Weeren (Dec. 2013, adapted from ft-task)

       A source model (souce.model image) or components list is converted into
       model visibilities that is inserted into the MODEL_DATA column or alternatively
       is stored  in the header of the MS to be served on the fly when requested.  This is
       needed to use more complicated sources than setjy provides; e.g resolved source
       or off centered sources in gaincal. Wprojection is included to deal with large FOV.
       If you do not need Wprojection than simply use the task ft.


	</description>
	<input>
		<param type="string" name="vis" mustexist="true">
			<description>Name of input visibility file (MS)</description>
			<value></value>
		</param>

		<param type="string" name="field">
			<description>Field selection</description>
			<value></value>
		</param>

		<param type="string" name="spw">
			<description>Spw selection</description>
			<value></value>
		</param>

		<param type="any" name="model">
			<description>Name of input model image(s)</description>
			<any type="variant"/>
			<value type="string"></value>
		</param>

		<param type="int" name="nterms">
			<description>Number of terms used to model the sky frequency dependence</description>
			<value>1</value>
		</param>

		<param type="string" name="reffreq" subparam="true">
			<description>Reference frequency (e.g. \'1.5e+9\' or \'1.5GHz\')</description>
			<value></value>
		</param>

                <param type="int" name="wprojplanes">
                <description>Number of w-projection planes for predict </description>
                <value>1</value>
                </param>

		<param type="string" name="complist">
			<description>Name of component list</description>
			<value></value>
		</param>

		<param type="bool" name="incremental">
			<description>Add to the existing model visibility?</description>
			<value>False</value>
		</param>
		<param type="bool" name="usescratch">
			<description>If True predicted  visibility  is stored in MODEL_DATA column</description>
			<value>False</value>
		</param>

		<constraints>
		        <when param="nterms">
			      <notequals type="int" value="1">
			           <default param="reffreq"><value type="string"></value></default>
			      </notequals>
			</when>
		</constraints>

	</input>
<returns type="void"/>

<example>

       A source model (souce.model image) or components list is converted into a
       model visibility that is inserted into the MODEL_DATA column.
       Wprojection is included to deal with large FOV.
       If you do not need Wprojection than simply use the task ft.

       TASK created by R. J. van Weeren (adapted from ft-task)




       Keyword arguments:
       vis -- Name of input visibility file
              default: none; example: vis='ngc5921.ms'
       field -- Field name list
               default: '' ==> all
               NOTE: BUT, only one source can be specified in a multi-source vis.
               field = '1328+307'  specifies source '1328+307'
               field = '4' specified field with index 4
       spw -- Spw selection
               default: spw = '' (all spw)
       model -- Name of input model image
               default: '' ==> None;
               example: model='/usr/lib/casapy/data/nrao/VLA/CalModels/3C286_X.im'
               Note: The model visibilities are scaled from the model frequency
                     to the observed frequency of the data.
       nterms -- Number of terms used to model the sky frequency dependence
                 default: 1  ==> one model image is required
                 example : nterms=3  represents a 2nd order Taylor-polynomial in frequency
                           and should be used in conjuction with coefficient model images as
		           model=['xxx.model.tt0','xxx.model.tt1', 'xxx.model.tt2']
             reffreq -- Reference-frequency about which this Taylor-expansion is defined.
	                default: '' ==> reads the reference frequency from the model image
                        example : reffreq = '1.5GHz'
       wprojplanes is the number of pre-computed w-planes used for
                   the W-Projection algorithm.  wprojplanes=1 disables
       complist -- Name of component list
               default: None; ; example: complist='test.cl'
               component lists are difficult to make.
       incremental -- Add model visibility to the existing model visibilties stored in the MS
               default: False; example: incremental=True
       usescratch  -- if True model visibilities will be stored in the scratch column
                            MODEL_DATA; when false the model visibilities will be generated
                            on the fly (this mode may save some disk space equivalent to
			    the volume of the observed data).
                            default: False; example usescratch=True






 </example>
 </task>
 </casaxml>
