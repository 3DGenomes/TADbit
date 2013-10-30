"""
30 Oct 2013


"""

def center_of_mass(model):
    """
    returns the center of mass of a given model
    """
    r_x = sum(model['x'])/len(model['x'])
    r_y = sum(model['y'])/len(model['y'])
    r_z = sum(model['z'])/len(model['z'])
    return dict(('x', r_x), ('y', r_y), ('z', r_z))


def radius_of_gyration(model):
    """
    returns the radius of gyration for the components of the tensor
    """
    com = center_of_mass(model)
    dist = sum([((model['x'][i] - com['x']) *
                 (model['y'][i] - com['y']) *
                 (model['z'][i] - com['z']))
                for i in xrange(len(model['x']))]) / len(model['x'])
    return dist
