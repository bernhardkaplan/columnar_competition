import numpy as np
import numpy.random as rnd
import os

def extract_trace(d, gid):
    """
    d : voltage trace from a saved with compatible_output=False
    gid : cell_gid
    """
    mask = gid * np.ones(d[:, 0].size)
    indices = mask == d[:, 0]
    time_axis, volt = d[indices, 1], d[indices, 2]
    return time_axis, volt


def get_local_indices(pop, offset=0):
    """
    Returns the list of indices (not IDs) local to the MPI node
    of a population
    ATTENTION: the ids returned by pop.all() assigned to the cell, depends _naturally_ on the 
    time of creation of the respective population --> important for offset
    """
    list_of_locals = []
    for tgt_id in pop.all():
        tgt = int(tgt_id) - offset - 1 # IDs are 1 aligned
        print tgt_id, tgt, pop.is_local(tgt_id)
        if pop.is_local(tgt_id) and (tgt < pop.size):
            list_of_locals.append(tgt)
    return list_of_locals


def convert_to_url(fn):
    """
    Needed for NeuroTools ParameterSet class
    """
    p = os.path.realpath('.')
    s = 'file://%s/%s' % (p, fn)
    return s
