import json
import numpy as np
import numpy.random as rnd
import os
from NeuroTools import parameters as ntp
import utils

class parameter_storage(object):
    """
    This class contains the simulation parameters in a dictionary called params.
    """

    def __init__(self):

        self.params = {}
        self.set_default_params()
        self.set_filenames()


    def set_default_params(self):

        self.params['simulator'] = 'nest' # 'brian' #

        self.params['n_hc'] = 3
        self.params['n_mc_per_hc'] = 5
        self.params['n_inh_per_hc'] = 10
        self.params['n_exc_per_mc'] = 10

        self.params['n_exc'] = self.params['n_hc'] * self.params['n_mc_per_hc'] * self.params['n_exc_per_mc']
        self.params['n_inh'] = self.params['n_hc'] * self.params['n_inh_per_hc']
        self.params['n_cells'] = self.params['n_exc'] + self.params['n_inh']


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
        self.params['p_ei'] = 0.10
        self.params['p_ie'] = 0.10

        # w_mean, w_sigma
        self.params['w_exc_within_mc'] = 0.005
        self.params['w_sigma'] = 0.33 # *100 % of the mean value

        self.params['w_ei'] = 0.005 #
        self.params['w_ie'] = 0.005 #


        # delays
        self.params['standard_delay'], self.params['standard_delay_sigma'] = 1, .1
        self.params['delay_range'] = (.1, 2.)

        ##################
        # INPUT PARAMS
        ##################
        self.params['stim_rate'] = 100.
        self.params['stim_start'] = 100.
        self.params['stim_duration'] = 100.
        self.params['stim_stop'] = self.params['stim_start'] + self.params['stim_duration']
        self.params['w_stim'] = 0.002
        self.params['stim_delay_sigma'] = 5. # [ms] how accurate / precise one minicolumn sees the input stimulus
        self.params['p_stim_exc'] = .25


        ########################
        # SIMULATION PARAMETERS
        ########################
        self.params['sim_stop'] = 1000.  # [ms] 
        self.params['sim_dt'] = 0.1  # [ms] 


        #########
        # SEEDS 
        ########
        self.params['rng_seeds_seed'] = 1234
        self.params['seed_conn'] = 0
        self.params['stim_seed'] = 4312


        ##############
        # RECORDING
        ##############
        self.params['n_exc_to_record_per_mc'] = self.params['n_exc_per_mc']
#        self.params['n_exc_to_record_per_mc'] = 5
        self.params['n_inh_to_record_per_hc'] = self.params['n_inh_per_hc']
#        self.params['n_inh_to_record_per_hc'] = 5
        assert (self.params['n_exc_to_record_per_mc'] <= self.params['n_exc_per_mc']), 'Can not record more cells than exist in a minicolumn'
        assert (self.params['n_inh_to_record_per_hc'] <= self.params['n_inh_per_hc']), 'Can not record more cells than exist in a minicolumn'


    def set_folder_name(self, folder_name):

        if folder_name == None:
            self.params['folder_name'] = 'Results/'
        else:
            self.params['folder_name'] = folder_name



    def set_filenames(self, folder_name=None):

        self.set_folder_name(folder_name) # make sure that self.params['folder_name'] is set

        self.params['input_folder'] = "%sInputSpikeTrains/"   % self.params['folder_name']# folder containing the input spike trains for the network generated from a certain stimulus
        self.params['spiketimes_folder'] = "%sSpikes/" % self.params['folder_name']
        self.params['volt_folder'] = "%sVoltageTraces/" % self.params['folder_name']
        self.params['parameters_folder'] = "%sParameters/" % self.params['folder_name']
        self.params['connections_folder'] = "%sConnections/" % self.params['folder_name']
        self.params['figures_folder'] = "%sFigures/" % self.params['folder_name']
        self.params['tmp_folder'] = "%stmp/" % self.params['folder_name']
        self.params['data_folder'] = '%sData/' % (self.params['folder_name']) # for storage of analysis results etc

        self.params['folder_names'] = [self.params['folder_name'], \
                            self.params['spiketimes_folder'], \
                            self.params['volt_folder'], \
                            self.params['parameters_folder'], \
                            self.params['connections_folder'], \
                            self.params['figures_folder'], \
                            self.params['tmp_folder'], \
                            self.params['data_folder'], \
#                            self.params['training_input_folder'], \
                            self.params['input_folder']] # to be created if not yet existing


        self.params['file_names'] = {}
        self.params['params_fn'] = '%ssimulation_parameters.info' % (self.params['parameters_folder'])
        self.params['params_fn_json'] = '%ssimulation_parameters.json' % (self.params['parameters_folder'])
        self.params['file_names']['input_st_fn_base'] = "%sstim_spike_train_" % self.params['input_folder']# input spike trains filename base
        self.params['file_names']['exc_spiketimes_fn_base'] = '%sexc_spikes_' % self.params['spiketimes_folder']
        self.params['file_names']['exc_nspikes_fn_merged'] = '%sexc_nspikes' % self.params['spiketimes_folder']
        self.params['file_names']['exc_volt_fn_base'] = '%sexc_volt' % self.params['volt_folder']
        self.params['file_names']['inh_volt_fn_base'] = '%sinh_volt' % self.params['volt_folder']

        #####################
        # FIGURE FILENAMES
        #####################
        self.params['rasterplot_exc_fig'] = '%srasterplot_exc.png' % (self.params['figures_folder'])
        self.params['rasterplot_inh_fig'] = '%srasterplot_inh.png' % (self.params['figures_folder'])


    def write_parameters_to_file(self, fn=None):
        """
        The NeuroTools ParameterSet save is not necessary, and could be removed but is 
        currently kept for readability issues (json files don't have linebreaks).
        """
        if not (os.path.isdir(self.params['folder_name'])):
            print 'Creating folder:\n\t%s' % self.params['folder_name']
            self.create_folders()

        if fn == None:
            fn = self.params['params_fn']
        print 'Writing parameters to: %s' % (fn)

        self.ParamSet = ntp.ParameterSet(self.params)
        fn = utils.convert_to_url(fn)
        self.ParamSet.save(fn)

        output_file = file(self.params['params_fn_json'], 'w')
        d = json.dump(self.params, output_file)


    def create_folders(self):
        """
        Create folders if not already existant
        """

        for f in self.params['folder_names']:
            if not os.path.exists(f):
                print 'Creating folder:\t%s' % f
                os.system("mkdir %s" % (f))

