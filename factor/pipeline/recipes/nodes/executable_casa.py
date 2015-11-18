#                                                          LOFAR PIPELINE SCRIPT
#
#                                           running an executable with arguments
#                                                         Stefan Froehlich, 2014
#                                                      s.froehlich@fz-juelich.de
# ------------------------------------------------------------------------------

from __future__ import with_statement
from subprocess import CalledProcessError
import os
import shutil
import sys
import errno
import tempfile

from lofarpipe.support.pipelinelogging import CatchLog4CPlus
from lofarpipe.support.pipelinelogging import log_time
from lofarpipe.support.utilities import create_directory
from lofarpipe.support.utilities import catch_segfaults
from lofarpipe.support.lofarnode import LOFARnodeTCP
from lofarpipe.support.parset import Parset


class executable_casa(LOFARnodeTCP):
    """
    Basic script for running an executable with arguments.
    """

    def run(self, infile, executable, args, kwargs, work_dir='/tmp', parsetasfile=False, args_format='', environment=''):
        """
        This function contains all the needed functionality
        """
        # Debugging info
        self.logger.debug("infile            = %s" % infile)
        self.logger.debug("executable        = %s" % executable)
        self.logger.debug("working directory = %s" % work_dir)
        self.logger.debug("arguments         = %s" % args)
        self.logger.debug("arg dictionary    = %s" % kwargs)
        self.logger.debug("environment       = %s" % environment)

        self.environment.update(environment)

        # Time execution of this job
        with log_time(self.logger):
            if infile[0] == '[':
                infiles = [ms.strip(" []\'\"") for ms in infile.split(',')]
                reffile = infiles[0]
            else:
                reffile = infile

            if os.path.exists(reffile):
                self.logger.info("Processing %s" % reffile)
            else:
                self.logger.error("Dataset %s does not exist" % reffile)
                return 1
 
            # Check if executable is present
            if not os.access(executable, os.X_OK):
                self.logger.error("Executable %s not found" % executable)
                return 1

            # hurray! race condition when running with than one process on one filesystem
            if not os.path.isdir(work_dir):
                try:
                    os.mkdir(work_dir, )
                except OSError as exc:  # Python >2.5
                    if exc.errno == errno.EEXIST and os.path.isdir(work_dir):
                        pass
                    else:
                        raise

            #print 'KWARGS: ', kwargs
            if not parsetasfile:
                self.logger.error("Nodescript \"executable_casa.py\" requires \"parsetasfile\" to be True!")
                return 1
            else:
                nodeparset = Parset()
                sublist = []
                for k, v in kwargs.items():
                    nodeparset.add(k, v)
                    if str(k).find('.'):
                        #print 'DOTPOS: ',str(k).find('.')
                        #print 'SPLIT: ', str(k).split('.')[0]
                        #print 'SPLIT: ', str(k).split('.')[1]
                        if not str(k).split('.')[0] in sublist:
                            sublist.append(str(k).split('.')[0])
                #print 'SUBPARSETLIST: ', sublist

                #subpar = Parset()
                #quick hacks below. for proof of concept.
                subparsetlist = []
                casastring = ''
                for sub in sublist:
                    subpar = nodeparset.makeSubset(nodeparset.fullModuleName(sub) + '.')
                    #print 'SUBPAR: ',subpar.keys()
                    casastring = sub + '('
                    for k in subpar.keys():
                        #print 'SUBPARSET: ',k ,' ',subpar[k]
                        #args.append('--' + k + '=' + subpar[k])
                        # Find type for parameter
                        if str(subpar[k]).find('/') == 0:
                            casastring += str(k) + '=' + "'" + str(subpar[k]) + "'" + ','
                        elif str(subpar[k]).find('casastr/') == 0:
                            casastring += str(k) + '=' + "'" + str(subpar[k]).strip('casastr/') + "'" + ','
                        elif str(subpar[k]).lower() == 'false' or str(subpar[k]).lower() == 'true':
                            casastring += str(k) + '=' + str(subpar[k]) + ','
                        else:
                            # Test if int/float or list of int/float
                            try:
                                self.logger.info('value: {}'.format(subpar[k]))
                                test = float(str(subpar[k]))
                                is_int_float = True
                            except:
                                is_int_float = False
                            if is_int_float:
                                casastring += str(k) + '=' + str(subpar[k]) + ','
                            else: 
                                if '[' in str(subpar[k]) or '(' in str(subpar[k]):
                                    # Check if list of int/float or strings
                                    list_vals = [f.strip() for f in str(subpar[k]).strip('[]()').split(',')]
                                    is_int_float = True
                                    for list_val in list_vals:
                                        try:
                                            test = float(list_val)
                                        except:
                                            is_int_float = False
                                            break
                                    if is_int_float:
                                        casastring += str(k) + '=' + str(subpar[k]) + ','
                                    else:
                                        casastring += str(k) + '=' + '[{}]'.format(','.join(["'"+list_val+"'" for list_val in list_vals])) + ','
                                else:
                                    # Simple string
                                    casastring += str(k) + '=' + "'" + str(subpar[k]) + "'" + ','

                    casastring = casastring.rstrip(',')
                    casastring += ')\n'
                #print 'CASASTRING:'
                #print casastring
                # 1) return code of a casapy is not properly recognized by the pipeline
                # wrapping in shellscript works for succesful runs.
                # failed runs seem to hang the pipeline...
                # 2) casapy can not have two instances running from the same directory.
                # create tmp dirs
                casapydir = tempfile.mkdtemp(dir=work_dir)
                if casastring != '':
                    casafilename = os.path.join(work_dir, os.path.basename(reffile) + '.casacommand.py')
                    casacommandfile = open(casafilename, 'w')
                    casacommandfile.write(casastring)
                    casacommandfile.close()
                    args.append(casafilename)
                somename = os.path.join(work_dir, os.path.basename(reffile) + '.casashell.sh')
                commandstring = ''
                commandstring += executable
                for item in args:                    
                    if str(item).find(' ')>-1 or str(item).find('[')>-1:
                        commandstring += ' "' + item + '"'
                    else:
                        commandstring += ' ' + item

                #print 'COMMANDSTRING: ',commandstring
                crap = open(somename, 'w')
                crap.write('#!/bin/bash \n')
                crap.write('echo "Trying CASAPY command" \n')
                crap.write(commandstring + ' >& casa.log\n')
                crap.close()

                import stat
                st = os.stat(somename)
                os.chmod(somename, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

            try:
            # ****************************************************************
            # Run
                cmd = [somename]
                with CatchLog4CPlus(
                    casapydir,
                    self.logger.name + "." + os.path.basename(reffile),
                    os.path.basename(executable),
                ) as logger:
                    # Catch segfaults and retry
                    catch_segfaults(
                        cmd, casapydir, self.environment, logger
                    )
            except CalledProcessError, err:
                # CalledProcessError isn't properly propagated by IPython
                self.logger.error(str(err))
                return 1
            except Exception, err:
                self.logger.error(str(err))
                return 1

        # We need some signal to the master script that the script ran ok.
        self.outputs['ok'] = True
        return 0

if __name__ == "__main__":
    #   If invoked directly, parse command line arguments for logger information
    #                        and pass the rest to the run() method defined above
    # --------------------------------------------------------------------------
    jobid, jobhost, jobport = sys.argv[1:4]
    sys.exit(executable_casa(jobid, jobhost, jobport).run_with_stored_arguments())
