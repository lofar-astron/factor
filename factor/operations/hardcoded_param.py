# this file contains some hardcoded parameters
# that can be easily edited by hand. The dict names
# refer to a specific operation. Actions should have all
# the parameters passed, not hardcoded

init_subtract = {
'imagerh' : {'niter': 40000,
        'imsize': 6144,
        'cell': '7.5arcsec',
        'uvrange': "0.08~7.0klambda",
        'threshpix': 4,
        'threshisl': 2.5,
        'nterms' : 1},
'imagerl' : {'niter' : 10000,
        'imsize': 4800,
        'cell': '25arcsec',
        'uvrange': "0.08~2.0klambda",
        'threshisl': 5,
        'threshpix': 5,
        'nterms': 1}
}