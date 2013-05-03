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

        self.params['n_hc'] = 3
        self.params['n_mc_per_hc'] = 5
        self.params['n_inh_per_hc'] = 10
        self.params['n_exc_per_mc'] = 10

        ###################
        # CELL PARAMETERS #
        ###################
        self.params['tau_syn_exc'] = 5.0 # 10.
        self.params['tau_syn_inh'] = 10.0 # 20.
        self.params['cell_params_exc'] = {'cm':1.0, 'tau_refrac':1.0, 'v_thresh':-50.0, \
                'tau_syn_E': self.params['tau_syn_exc'], 'tau_syn_I':self.params['tau_syn_inh'], \
                'tau_m' : 10., 'v_reset' : -70., 'v_rest':-70}
        self.params['cell_params_inh'] = {'cm':1.0, 'tau_refrac':1.0, 'v_thresh':-50.0, \
                'tau_syn_E': self.params['tau_syn_exc'], 'tau_syn_I':self.params['tau_syn_inh'], \
                'tau_m' : 10., 'v_reset' : -70., 'v_rest':-70}

        ##########################
        # CONNECTIVITY PARAMETERS 
        ##########################
        # probabilities
        self.params['p_exc_within_mc'] = 0.10

        # w_mean, w_sigma
        self.params['w_exc_within_mc'] = 0.005
        self.params['w_sigma'] = 0.25 # *100 % of the mean value


        #########
        # SEEDS 
        ########
        self.params['seed_conn'] = 0

