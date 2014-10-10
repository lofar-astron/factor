"""
Operation: init_subtract
An example operation, to be used as a template
"""

from factor.lib.operation import operation
import factor.actions as a

class init_subtract(operation):

    def run(self):

        from factor.operations.hardcoded_param import init_subtract as p

        mss = self.parset['mss']

        # this is a typical multi-thread block
        self.log.info('High-res imaging...')
        actions = [ a.imager_mask.imager_mask(self.name, ms, p['imagerh'], 'init_hires') for ms in mss ]
        self.s.run_action_parallel( actions )

        self.log.info('Make high-res skymodel...')
        actions = [ a.make_skymodel.make_skymodel(self.name, ms, p['xxx'], 'init_hires') for ms in mss ]
        self.s.run_action_parallel( actions )

        self.log.info('Subtract high-res skymodel...')
        actions = [ a.subtract.subtract(self.name, ms, p['xxx'], 'init_hires') for ms in mss ]
        self.s.run_action_parallel( actions )

        self.log.info('Average...')
        actions = [ a.avg.avg(self.name, ms) for ms in mss ]
        self.s.run_action_parallel( actions )

        self.log.info('Low-res imaging...')
        actions = [ a.imager_mask.imager_mask(self.name, ms, p['imagerl'], 'init_lowres') for ms in mss ]
        self.s.run_action_parallel( actions )

        self.log.info('Make low-res skymodel...')
        actions = [ a.make_skymodel.make_skymodel(self.name, ms, p['xxx'], 'init_lowres') for ms in mss ]
        self.s.run_action_parallel( actions )

        self.log.info('Subtract low-res skymodel...')
        actions = [ a.subtract.subtract(self.name, ms, p['xxx'], 'init_lowres') for ms in mss ]
        self.s.run_action_parallel( actions )

        self.log.info('Merge skymodels...')
        actions = [ a.merge_skymodel.merge_skymodel(self.name, ms) for ms in mss ]
        self.s.run_action_parallel( actions )
