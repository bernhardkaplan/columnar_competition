#pierre
import numpy as np
import simulation_parameters
ps = simulation_parameters.parameter_storage()
params = ps.params
exec("from pyNN.%s import *" % params['simulator'])
#import pyNN


class NetworkModel(object):

    def __init__(self, params):
        self.params = params

    def setup(self, load_tuning_prop=False):

        self.projections = {}
        self.projections['ee'] = []
        self.projections['ei'] = []
        self.projections['ie'] = []
        self.projections['ii'] = []



if __name__ == '__main__':

    NM = NetworkModel(ps.params)
    NM.setup()
