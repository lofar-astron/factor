"""
Definition of context managers (with statements) used for operations
"""

class op_timer:
    """
    context manager used to time the operations
    """
    import time

    def __enter__(self):
        self.start = time.clock()

    def __exit__(self, type, value, tb):
        if type is not None:
            pass
        elapsed = (time.clock() - self.start)
        logging.debug('Time for operation: %i sec'.format(elapsed))

class op_init:
    """
    context manager used to initialize an operations
    """
    import factor.operations as op
    operations = {
        'tessellation': op.tessellation
        'init_subtract': op.init_subtract
        'facet_setup': op.facet_setup
        'facet_selfcal': op.facet_selfcal
        'facet_calib': op.facet_calib
        'facet_image': op.facet_image
        'facet_final': op.facet_final
        'facet_comb_image': op.facet_comb_image
        'final_image': op.final_image
        }

    def __enter__(self, name, parset, direction = None):
        self.name = name
        self.parset = parset
        self.direction = direction
        #obj = getattr(operations[name], name)
        obj = getattr(getattr(op, name), name)
        if direction != None: return obj(parset, direction)
        else: return obj(parset)

    def __exit__(self, type, value, tb):
        if type is not None:
            logging.error('Problems running operation: %s'.format(self.name))
            raise(type)
        logging.info('Operation %s, terminated successfully.'.format(self.name))



