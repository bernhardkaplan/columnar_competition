#!/usr/bin/python
"""
A single IF neuron with exponential, conductance-based synapses, fed by two
spike sources.

Run as::

   $ python IF_cond_exp.py <simulator>

where <simulator> is 'neuron', 'nest', etc

Andrew Davison, UNIC, CNRS
May 2006

Expected results in
  .. figure::  ./examples/results/IF_cond_exp.png


$Id: IF_cond_exp.py 917 2011-01-31 15:23:34Z apdavison $
"""

#!/usr/bin/python
#simulator_name = 'spiNNaker'
simulator_name = 'brian'
#simulator_name = 'nest'

from pyNN.utility import get_script_args
from pyNN.errors import RecordingError
from pyNN.random import NumpyRNG, RandomDistribution
import numpy as np

exec("import pyNN.%s as p" % simulator_name)



def cortical_cluster(Ne = 5, Ni = 1, Nc = 10):

    p.setup(timestep=1.0,min_delay=1.0,max_delay=10.0, db_name='if_cond.sqlite')
    
    uniformDistr = RandomDistribution('uniform', [0,1])
    uniformDistr_inh = RandomDistribution('uniform',[-1,0]) 
    #cell_params = {     'i_offset' : .1,    'tau_refrac' : 3.0, 'v_rest' : -65.0,
    #                    'v_thresh' : -51.0,  'tau_syn_E'  : 2.0,
     #                   'tau_syn_I': 5.0,    'v_reset'    : -70.0,
     #                   'e_rev_E'  : 0.,     'e_rev_I'    : -80.}

    cell_params = {     'i_offset' : .1,    'tau_refrac' : 3.0, 'v_rest' : -65.0,
                        'v_thresh' : -51.0,  'tau_syn_E'  : 2.0,
                        'tau_syn_I': 5.0,    'v_reset'    : -70.0}
                        
    source_times = {'spike_times':[i for i in range(5,105,10)]}
    exc_pop = []
    inh_pop = []
    spike_source  = []
    proj_input = []
    proj_exc_to_inh = []
    proj_inh_to_exc = []
    print p.SpikeSourceArray
    
    for i in range (Nc):
        exc_pop.append(p.Population(Ne, p.IF_curr_exp, cell_params, label='exc_pop'))

        inh_pop.append(p.Population(Ni, p.IF_curr_exp, cell_params, label='inh_pop'))

        spike_source.append(p.Population(Ne, p.SpikeSourceArray, source_times, label='spike_sourceE'))

        #spike_sourceI = p.Population(1, p.SpikeSourceArray, {'spike_times': [[i for i in range(155,255,10)],]}, label='spike_sourceI')

        proj_input.append(p.Projection(spike_source[i], exc_pop[i], p.AllToAllConnector(weights=uniformDistr, delays=1), target='excitatory'))
        #proj_exc_to_exc[i] = p.Projection(exc_pop[i], exc_pop[i], p.AllToAllConnector(weights=uniformDistr, delays=1), target='excitatory')
        proj_exc_to_inh.append(p.Projection(exc_pop[i], inh_pop[i], p.AllToAllConnector(weights=uniformDistr, delays=1), target='excitatory'))
        proj_inh_to_exc.append(p.Projection(inh_pop[i], exc_pop[i], p.AllToAllConnector(weights=uniformDistr_inh, delays=1), target='inhibitory'))
        #connE = p.Projection(spike_sourceE, ifcell, p.OneToOneConnector(weights=0.006, delays=2), target='excitatory')
        #connI = p.Projection(spike_sourceI, ifcell, p.OneToOneConnector(weights=0.02, delays=4), target='inhibitory')
            
        exc_pop[i].record_v()
        exc_pop[i].record()

        inh_pop[i].record_v()
        inh_pop[i].record()
        
    proj_column = [[0 for i in range(Nc)] for j in range(Nc)]
    for i in range (Nc):
        for j in range (Nc):
            proj_column[i][j] = p.Projection(exc_pop[i], exc_pop[j], p.AllToAllConnector(weights=uniformDistr, delays=1), target='excitatory')


    p.run(200.0)
    #ifcell = exc_pop[1].all_cells[0]
    #print ifcell
    #ifcell.print_v('results/IF_cond_exp_%s.v' % simulator_name)
    #recorded_v =  ifcell.get_v()

    #import pylab
    #f = pylab.figure()
    #f.add_subplot(211)
    #pylab.plot([ i[2] for i in recorded_v ])

    #f.add_subplot(212)
    #pylab.plot([ i[2] for i in recorded_gsyn ], color='green')

    #pylab.show()
    #p.end()

