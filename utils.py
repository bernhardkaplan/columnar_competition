
def get_local_indices(pop, offset=0):
    """
    Returns the list of indices (not IDs) local to the MPI node
    of a population
    """
    list_of_locals = []
    for tgt_id in pop.all():
        tgt = int(tgt_id) - offset - 1 # IDs are 1 aligned
        print tgt_id, tgt, pop.is_local(tgt_id)
        if pop.is_local(tgt_id) and (tgt < pop.size):
            list_of_locals.append(tgt)
    return list_of_locals

