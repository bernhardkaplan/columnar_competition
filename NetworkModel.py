import numpy as np
import utils
import simulation_parameters
ps = simulation_parameters.parameter_storage()
params = ps.params
exec("from pyNN.%s import *" % params['simulator'])
#import pyNN


try:
#    I_fail_because_I_do_not_want_to_use_MPI
    from mpi4py import MPI
    USE_MPI = True
    comm = MPI.COMM_WORLD
    pc_id, n_proc = comm.rank, comm.size
    print "USE_MPI:", USE_MPI, 'pc_id, n_proc:', pc_id, n_proc
except:
    USE_MPI = False
    pc_id, n_proc, comm = 0, 1, None
    print "MPI not used"

def get_local_indices(pop, offset=0):
    """
    Returns the list of indices (not IDs) local to the MPI node
    of a population
    """
    list_of_locals = []
    for tgt_id in pop.all():
        tgt = int(tgt_id) - offset - 1 # IDs are 1 aligned
        print tgt_id, tgt, pop.is_local(tgt_id), pc_id
        if pop.is_local(tgt_id) and (tgt < pop.size):
            list_of_locals.append(tgt)
    return list_of_locals


class NetworkModel(object):

    def __init__(self, params):
        self.params = params
        self.rng_conn = NumpyRNG(seed = self.params['seed_conn'], parallel_safe=True) #if True, slower but does not depend on number of nodes

    def setup(self, load_tuning_prop=False):

        self.projections = {}
        self.projections['ee'] = []
        self.projections['ei'] = []
        self.projections['ie'] = []
        self.projections['ii'] = []

    def create(self):
        
        self.exc_populations = []
        self.inh_populations = []
        for hc in xrange(self.params['n_hc']):
            inh_pop = Population(self.params['n_inh_per_hc'], IF_cond_exp, self.params['cell_params_inh'], label="inh_pop")
            self.inh_populations.append(inh_pop)
            list_of_mcs = []
            for mc in xrange(self.params['n_mc_per_hc']):
                exc_pop = Population(self.params['n_exc_per_mc'], IF_cond_exp, self.params['cell_params_exc'], label='exc_cells')
                list_of_mcs.append(exc_pop)
            print 'debug', len(self.exc_populations)
            self.exc_populations.append(list_of_mcs)
            for mc in xrange(self.params['n_mc_per_hc']):
                print 'hc %d mc %d local exc ids' % (hc, mc), len(get_local_indices(self.exc_populations[-1][mc]))
#        self.local_idx_exc = 
#        self.local_idx_inh = utils.get_local_indices(self.inh_pop, offset=self.params['n_exc'])

    def connect(self):
        
        self.connect_inside_mc()


    def connect_inside_mc(self):

        w_dist = RandomDistribution('normal',
                (self.params['w_exc_within_mc'], self.params['w_sigma'] * self.params['w_exc_within_mc']),
                rng=self.rng_conn,
                constrain='redraw',
                boundaries=(0, w_mean * 10.))

        conn = FixedProbabilityConnector(self.params['p_exc_within_mc'], allow_self_connections=False, weights=w_dist)
        for hc in xrange(self.params['n_hc']):
            for mc in self.exc_populations[hc]:
                prj = Projection(mc, mc, conn)
                self.projections['ee'].append(prj)



if __name__ == '__main__':

    NM = NetworkModel(ps.params)
    NM.setup()
    NM.create()
