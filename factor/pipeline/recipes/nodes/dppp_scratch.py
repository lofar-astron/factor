#                                                          LOFAR PIPELINE SCRIPT
#
#                                             running DPPP in a temp scratch dir
# ------------------------------------------------------------------------------

from __future__ import with_statement
from subprocess import CalledProcessError
import os
import shutil
import sys
import errno
import subprocess

from lofarpipe.support.pipelinelogging import CatchLog4CPlus
from lofarpipe.support.pipelinelogging import log_time
from lofarpipe.support.utilities import catch_segfaults
from lofarpipe.support.lofarnode import LOFARnodeTCP
from lofarpipe.support.parset import Parset


class dppp_scratch(LOFARnodeTCP):
    """
    Basic script for running DPPP in a scratch directory.
    """

    def run(self, infile, executable, args, kwargs, work_dir='/tmp', parsetasfile=True, args_format='', environment=''):
        """
        This method contains all the needed functionality
        """

        # Debugging info
        self.logger.debug("infile            = %s" % infile)
        self.logger.debug("executable        = %s" % executable)
        self.logger.debug("working directory = %s" % work_dir)
        self.logger.debug("arguments         = %s" % args)
        self.logger.debug("arg dictionary    = %s" % kwargs)
        self.logger.debug("environment       = %s" % environment)

        self.environment.update(environment)

        self.work_dir = work_dir
        self.infile = infile
        self.executable = executable

        if 'msin' in kwargs:
            self.msin_original = kwargs['msin'].rstrip('/')
            kwargs.pop('msin')

        if 'msout' in kwargs:
            self.msout_original = kwargs['msout'].rstrip('/')
            kwargs.pop('msout')
        else:
            self.msout_original = '.'
        if self.msout_original == '.':
            # Input file is modified
            self.msout_destination_dir = os.path.dirname(self.msin_original)
        else:
            # New output is created
            self.msout_destination_dir = os.path.dirname(self.msout_original)

        if 'local_scratch_dir' in kwargs:
            self.scratch_dir = kwargs['local_scratch_dir'].rstrip('/')
            kwargs.pop('local_scratch_dir')

        self.logger.info(self.scratch_dir)

        # Set up scratch paths
        self.msin_scratch = os.path.join(self.scratch_dir, os.path.basename(self.msin_original))
        args.append('msin=' + self.msin_scratch)
        if self.msout_original == '.':
            # Input file is modified
            self.msout_scratch = os.path.join(self.scratch_dir, os.path.basename(self.msin_original))
            args.append('msout=.')
        else:
            # New output is created
            self.msout_scratch = os.path.join(self.scratch_dir, os.path.basename(self.msout_original))
            args.append('msout=' + self.msout_scratch)

        self.copy_to_scratch()

        # Time execution of this job
        with log_time(self.logger):
            #if os.path.exists(infile):
            self.logger.info("Processing %s" % infile)

            # Check if script is present
            if not os.path.isfile(executable):
                self.logger.error("Executable %s not found" % executable)
                return 1

            # hurray! race condition when running with more than one process on one filesystem
            if not os.path.isdir(work_dir):
                try:
                    os.mkdir(work_dir, )
                except OSError as exc:  # Python >2.5
                    if exc.errno == errno.EEXIST and os.path.isdir(work_dir):
                        pass
                    else:
                        raise

            if parsetasfile:
                nodeparset = Parset()
                parsetname = os.path.join(work_dir, os.path.basename(infile) + '.parset')
                for k, v in kwargs.items():
                    nodeparset.add(k, v)
                nodeparset.writeFile(parsetname)
                #args.insert(0, parsetname)
                args.append(parsetname)

            try:
            # ****************************************************************
            #Run
                cmd = [executable] + args
                with CatchLog4CPlus(
                    work_dir,
                    self.logger.name + "." + os.path.basename(infile),
                    os.path.basename(executable),
                ) as logger:
                    # Catch segfaults and retry
                    catch_segfaults(
                        cmd, work_dir, self.environment, logger
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

        # Copy output data back to origin
        self.copy_to_origin()

        # Delete scratch directory
        self.cleanup()
        return 0

    def copy_to_scratch(self):
        self.logger.info("Copying input data to scratch directory")
        if os.path.exists(self.msin_scratch):
            shutil.rmtree(self.msin_scratch)
        args = ['-a', self.msin_original, self.scratch_dir]
        prog = '/usr/bin/rsync'
        self.execute(prog, args)

    def copy_to_origin(self):
        self.logger.info("Copying output data to original directory")
        args = ['-a', self.msout_scratch, self.msout_destination_dir ]
        prog = '/usr/bin/rsync'
        self.execute(prog, args)

    def cleanup(self):
        self.logger.info("Deleting scratch directory")
        shutil.rmtree(self.msin_scratch)
        shutil.rmtree(self.msout_scratch)

    def execute(self, executable, args):
        try:
        # ****************************************************************
        # Run
            cmd = [executable] + args
            with CatchLog4CPlus(
                self.work_dir,
                self.logger.name + "." + os.path.basename(self.infile),
                os.path.basename(self.executable),
            ) as logger:
                # Catch segfaults and retry
                catch_segfaults(
                    cmd, self.work_dir, self.environment, logger
                )
        except CalledProcessError, err:
            # CalledProcessError isn't properly propagated by IPython
            self.logger.error(str(err))
            return 1
        except Exception, err:
            self.logger.error(str(err))
            return 1

if __name__ == "__main__":
    #   If invoked directly, parse command line arguments for logger information
    #                        and pass the rest to the run() method defined above
    # --------------------------------------------------------------------------
    jobid, jobhost, jobport = sys.argv[1:4]
    sys.exit(dppp_scratch(jobid, jobhost, jobport).run_with_stored_arguments())
