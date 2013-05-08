import numpy as np
import utils
import simulation_parameters
ps = simulation_parameters.parameter_storage()
params = ps.params
exec("from pyNN.%s import *" % params['simulator'])
#import pyNN
import time
times = {}
t0 = time.time()


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

    def __init__(self, params, comm=None):
        """
        params : parameter dictionary
        comm : MPI communicator
        """
        self.params = params
        self.conn_rng = NumpyRNG(seed = self.params['seed_conn'], parallel_safe=True) #if True, slower but does not depend on number of nodes
        self.comm = comm
        if self.comm != None:
            self.pc_id, self.n_proc = self.comm.rank, self.comm.size
        else:
            self.pc_id, self.n_proc = 0, 1
        if self.comm != None:
            self.comm.Barrier()
        from pyNN.utility import Timer
        self.timer = Timer()
        self.timer.start()
        self.times = {} # dictionary to record connection times, simulation time, etc. ...

    def setup(self, load_tuning_prop=False):

        setup(timestep=self.params['sim_dt'], min_delay=self.params['delay_range'][0], max_delay=self.params['delay_range'][1], rng_seeds_seed=self.params['rng_seeds_seed'])
        self.projections = {}
        self.projections['ee'] = []
        self.projections['ei'] = []
        self.projections['ie'] = []
        self.projections['ii'] = []
        if self.comm != None:
            self.comm.Barrier()

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
        if self.comm != None:
            self.comm.Barrier()

    def connect_input(self, save_input_files=False):

        self.projections['stim'] = []
        poisson_parameters = {'duration' : self.params['stim_duration'], 'start' : self.params['stim_start'], 'rate' : self.params['stim_rate']}
        rng = NumpyRNG(seed=self.params['stim_seed'], parallel_safe=True)
        n_spikes = 1 + int(self.params['n_exc_per_mc'] * (self.params['stim_stop'] - self.params['stim_start']) * self.params['stim_rate']/1000.0)



        w_dist = RandomDistribution('normal',
                (self.params['w_stim'], self.params['w_sigma'] * self.params['w_stim']),
                rng=self.conn_rng,
                constrain='redraw',
                boundaries=(0, self.params['w_stim'] * 10.))

        delay_dist = RandomDistribution('normal',
                (self.params['stim_delay_sigma'], self.params['stim_delay_sigma']),
                rng=self.conn_rng,
                constrain='redraw',
                boundaries=(self.params['delay_range'][0], self.params['delay_range'][1]))


        # + 1 is to discard the last one which would always be at stim_stop for all populations
        for hc in xrange(self.params['n_hc']):
            for i_, mc in enumerate(self.exc_populations[hc]):
                mc_cnt = hc * self.params['n_mc_per_hc'] + i_
                # create a population of poisson distribution intput spikes
                spike_times = numpy.add.accumulate(rng.next(n_spikes, 'exponential',
                                                                [1000.0/self.params['stim_rate']], mask_local=False))
                # bring the n_spikes spiketimes into the right range
                spike_times += self.params['stim_start'] 
                spike_times /= spike_times.max()
                spike_times *= self.params['stim_stop']
                spike_times = spike_times[:-1]
                input_population  = Population(self.params['n_exc_per_mc'], SpikeSourceArray, {'spike_times': spike_times }, label="input")
#                conn = OneToOneConnector(self.params['w_stim'], delays=None)
#                conn = AllToAllConnector(weights=w_dist, delays=delay_dist)
                conn = FixedProbabilityConnector(self.params['p_stim_exc'], weights=w_dist, delays=delay_dist)
                prj = Projection(input_population, mc, conn, label='input-mc-prj-%d' % mc_cnt)

                if save_input_files:
                    fn = self.params['input_st_fn_base']
                    np.save(fn, spike_times)
        if self.comm != None:
            self.comm.Barrier()



    def connect_populations(self):
        
        t0 = time.time()
        self.connect_inside_mc()
        self.connect_mcs_to_inh_pop()
        self.connect_inh_pop_to_mcs()
        t1 = time.time()
        self.times['t_connect'] = t1 - t0


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
        if self.comm != None:
            self.comm.Barrier()
#        self.times['t_connect_mcs_to_inh_pop'] = self.timer.diff()



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
        if self.comm != None:
            self.comm.Barrier()
#        self.times['t_connect_inh_pop_to_mcs'] = self.timer.diff()


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
        if self.comm != None:
            self.comm.Barrier()
#        self.times['t_connect_inside_mc'] = self.timer.diff()

    def record(self, record_v=False):
        
        self.recorded_cells_exc = []
        self.recorded_cells_inh = []
        for hc in xrange(self.params['n_hc']):
            if record_v:
                if self.params['n_inh_to_record_per_hc'] < self.params['n_inh_per_hc']:
                    gids_to_record = np.random.randint(0, self.params['n_inh_to_record_per_hc'], self.params['n_inh_per_hc']),
                else:
                    gids_to_record = np.arange(0, self.params['n_inh_per_hc'], dtype=np.int)
                inh_pop_view = PopulationView(self.inh_populations[hc], gids_to_record, label='inh-sample-%d' % hc)
                inh_pop_view.record_v()
                self.recorded_cells_inh.append(inh_pop_view)
            for i_, mc in enumerate(self.exc_populations[hc]):
                mc_cnt = hc * self.params['n_mc_per_hc'] + i_
                mc.record() # record spikes
                if record_v:
                    if self.params['n_exc_to_record_per_mc'] < self.params['n_exc_per_mc']:
                        gids_to_record = np.random.randint(0, self.params['n_exc_per_mc'], self.params['n_exc_to_record_per_mc']),
                    else:
                        gids_to_record = np.arange(0, self.params['n_exc_per_mc'], dtype=np.int)

                    exc_pop_view = PopulationView(mc, gids_to_record, label='exc-sample-%d' % mc_cnt)
                    exc_pop_view.record_v()
                    self.recorded_cells_exc.append(exc_pop_view)
        
        if self.comm != None:
            self.comm.Barrier()
        self.times['t_record'] = self.timer.diff()

    def print_results(self, record_v=False):
        print 'print_results...'

        if record_v:
            for i_, recorded_mc in enumerate(self.recorded_cells_exc):
                fn = self.params['file_names']['exc_volt_fn_base'] + '%d.v' % i_
                recorded_mc.print_v(fn, compatible_output=False)
            for i_, recorded_inh_cells in enumerate(self.recorded_cells_inh):
                fn = self.params['file_names']['inh_volt_fn_base'] + '%d.v' % i_
                recorded_inh_cells.print_v(fn, compatible_output=False)


        for hc in xrange(self.params['n_hc']):
            for i_, mc in enumerate(self.exc_populations[hc]):
                mc_cnt = hc * self.params['n_mc_per_hc'] + i_
                mc.printSpikes(self.params['file_names']['exc_spiketimes_fn_base'] + '%d.ras' % mc_cnt)
        if self.comm != None:
            self.comm.Barrier()

        if self.pc_id == 0:
            self.times['t_all'] = 0.
            for k in self.times.keys():
                self.times['t_all'] += self.times[k]

            self.n_cells = {}
            self.n_cells['n_exc'] = self.params['n_exc']
            self.n_cells['n_inh'] = self.params['n_inh']
            self.n_cells['n_cells'] = self.params['n_cells']
            output = {'times' : self.times, 'n_cells_proc' : self.n_cells}
            print "Proc %d Simulation time: %d sec or %.1f min for %d cells (%d exc %d inh)" % (self.pc_id, self.times['t_sim'], (self.times['t_sim'])/60., self.params['n_cells'], self.params['n_exc'], self.params['n_inh'])
            print "Proc %d Full pyNN run time: %d sec or %.1f min for %d cells (%d exc %d inh)" % (self.pc_id, self.times['t_all'], (self.times['t_all'])/60., self.params['n_cells'], self.params['n_exc'], self.params['n_inh'])
        self.times['t_print'] = self.timer.diff()

    def run_sim(self):

        if self.pc_id == 0:
            print "Running simulation ... "
        run(self.params['sim_stop'])
        if self.comm != None:
            self.comm.Barrier()
        self.times['t_sim'] = self.timer.diff()

if __name__ == '__main__':

    ps.write_parameters_to_file()
    NM = NetworkModel(ps.params)
    NM.setup()
    NM.create()
    NM.connect_populations()
    NM.connect_input()
    NM.record(record_v=True)
    NM.run_sim()
    NM.print_results(record_v=True)
