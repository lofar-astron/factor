"""
Operation: init_subtract
An example operation, to be used as a template
"""

from factor.lib.operation import operation
import factor.actions as a

class init_subtract(operation):

    def run(self):

        import init_subtract from factor.operation.hardcoded_param as p

        mss = self.parset['mss']

        # this is a typical multi-thread block
        self.log.info('High-res imaging...')
        actions = [ a.imager_mask(ms, niter=p['niterh'], imsize=p['imsizeh'], cell=p['cellh'], uvrange=p['uvrangeh']) for ms in mss ]
        self.s.run_action_parallel( actions )

        self.log.info('Make high-res skymodel...')
        actions = [ a.make_skymodel(ms, ) for ms in mss ]
        self.s.run_action_parallel( actions )

        self.log.info('Subtract high-res skymodel...')
        actions = [ a.subtract(ms, ) for ms in mss ]
        self.s.run_action_parallel( actions )

        self.log.info('Average...')
        actions = [ a.avg(ms, ) for ms in mss ]
        self.s.run_action_parallel( actions )

        self.log.info('Low-res imaging...')
        actions = [ a.imager_mask(ms, niter=p['niterl'], imsize=p['imsizel'], cell=p['celll'], uvrange=p['uvrangel']) for ms in mss ]
        self.s.run_action_parallel( actions )

        self.log.info('Make low-res skymodel...')
        actions = [ a.make_skymodel(ms, ) for ms in mss ]
        self.s.run_action_parallel( actions )

        self.log.info('Subtract low-res skymodel...')
        actions = [ a.subtract(ms, ) for ms in mss ]
        self.s.run_action_parallel( actions )

        self.log.info('Merge skymodels...')
        actions = [ a.merge_skymodel(ms, ) for ms in mss ]
        self.s.run_action_parallel( actions )
