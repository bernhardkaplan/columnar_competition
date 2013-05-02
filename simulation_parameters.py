import json
import numpy as np
import numpy.random as rnd
import os
from NeuroTools import parameters as ntp
#import utils

class parameter_storage(object):
    """
    This class contains the simulation parameters in a dictionary called params.
    """

    def __init__(self):

        self.params = {}
        self.set_default_params()


    def set_default_params(self):

        self.params['simulator'] = 'nest' # 'brian' #

        self.params['n_hc'] = 2
        self.params['n_mc_per_hc'] = 2
        self.params['n_inh_per_hc'] = 1

