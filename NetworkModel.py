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


class NetworkModel(object):

    def __init__(self, params):
        self.params = params
        self.conn_rng = NumpyRNG(seed = self.params['seed_conn'], parallel_safe=True) #if True, slower but does not depend on number of nodes

    def setup(self, load_tuning_prop=False):

        setup(timestep=0.1, min_delay=self.params['delay_range'][0], max_delay=self.params['delay_range'][1], rng_seeds_seed=self.params['rng_seeds_seed'])
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
            self.exc_populations.append(list_of_mcs)

    def connect_input(self):

        poisson_parameters = {'duration' : self.params['stim_duration'], 'start' : self.params['stim_start'], 'rate' : self.params['stim_rate']}
        rng = NumpyRNG(seed=self.params['stim_seed'], parallel_safe=True)
        for hc in xrange(self.params['n_hc']):
            for i_, mc in enumerate(self.exc_populations[hc]):
                n_spikes = int(self.params['n_exc_per_mc'] * self.params['stim_stop'] * self.params['stim_rate']/1000.0)
                spike_times = numpy.add.accumulate(rng.next(n_spikes, 'exponential',
                                                                [1000.0/self.params['stim_rate']], mask_local=False))
                input_population  = Population(self.params['n_exc_per_mc'], SpikeSourceArray, {'spike_times': spike_times }, label="input")

                # create a spike source poisson 



    def connect_populations(self):
        
        self.connect_inside_mc()
        self.connect_mcs_to_inh_pop()
        self.connect_inh_pop_to_mcs()


    def connect_mcs_to_inh_pop(self):

        w_dist = RandomDistribution('normal',
                (self.params['w_ei'], self.params['w_sigma'] * self.params['w_ei']),
                rng=self.conn_rng,
                constrain='redraw',
                boundaries=(0, self.params['w_ei'] * 10.))

        delay_dist = RandomDistribution('normal',
                (self.params['standard_delay'], self.params['standard_delay_sigma']),
                rng=self.conn_rng,
                constrain='redraw',
                boundaries=(self.params['delay_range'][0], self.params['delay_range'][1]))

        conn = FixedProbabilityConnector(self.params['p_ei'], weights=w_dist, delays=delay_dist)

        for hc in xrange(self.params['n_hc']):
            inh_pop = self.inh_populations[hc]
            for i_, mc in enumerate(self.exc_populations[hc]):
                mc_cnt = hc * self.params['n_mc_per_hc'] + i_
                prj = Projection(mc, inh_pop, conn, label='mc-inh-prj-%d' % mc_cnt)
                self.projections['ei'].append(prj)


    def connect_inh_pop_to_mcs(self):

        w_dist = RandomDistribution('normal',
                (self.params['w_ie'], self.params['w_sigma'] * self.params['w_ie']),
                rng=self.conn_rng,
                constrain='redraw',
                boundaries=(0, self.params['w_ie'] * 10.))

        delay_dist = RandomDistribution('normal',
                (self.params['standard_delay'], self.params['standard_delay_sigma']),
                rng=self.conn_rng,
                constrain='redraw',
                boundaries=(self.params['delay_range'][0], self.params['delay_range'][1]))

        conn = FixedProbabilityConnector(self.params['p_ie'], weights=w_dist, delays=delay_dist)

        for hc in xrange(self.params['n_hc']):
            inh_pop = self.inh_populations[hc]
            for i_, mc in enumerate(self.exc_populations[hc]):
                mc_cnt = hc * self.params['n_mc_per_hc'] + i_
                prj = Projection(inh_pop, mc, conn, label='inh-mc-prj-%d' % mc_cnt)
                self.projections['ei'].append(prj)


    def connect_inside_mc(self):

        w_dist = RandomDistribution('normal',
                (self.params['w_exc_within_mc'], self.params['w_sigma'] * self.params['w_exc_within_mc']),
                rng=self.conn_rng,
                constrain='redraw',
                boundaries=(0, self.params['w_exc_within_mc'] * 10.))

        delay_dist = RandomDistribution('normal',
                (self.params['standard_delay'], self.params['standard_delay_sigma']),
                rng=self.conn_rng,
                constrain='redraw',
                boundaries=(self.params['delay_range'][0], self.params['delay_range'][1]))

        conn = FixedProbabilityConnector(self.params['p_exc_within_mc'], allow_self_connections=False, \
                weights=w_dist, delays=delay_dist)
        
        for hc in xrange(self.params['n_hc']):
            for i_, mc in enumerate(self.exc_populations[hc]):
                mc_cnt = hc * self.params['n_mc_per_hc'] + i_
                prj = Projection(mc, mc, conn, label='mc-mc-prj-%d' % mc_cnt)
                self.projections['ee'].append(prj)




if __name__ == '__main__':

    NM = NetworkModel(ps.params)
    NM.setup()
    NM.create()
    NM.connect_populations()
    NM.connect_input()
