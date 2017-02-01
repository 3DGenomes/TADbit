"""
25 Oct 2016


"""
from pytadbit.modelling.structuralmodel import StructuralModel

class LAMMPSmodel(StructuralModel):
    """
    A container for the LAMMPS modelling results. The container is a dictionary
    with the following keys:
    
    - rand_init: Random number generator feed (needed for model reproducibility)
    - x, y, z: 3D coordinates of each particles. Each represented as a list

    """
    def __str__(self):
        try:
            return ('LAMMPS model ranked %s (%s particles) with: \n' +
                    ' - random initial value: %s\n' +
                    ' - first coordinates:\n'+
                    '        X      Y      Z\n'+
                    '  %7s%7s%7s\n'+
                    '  %7s%7s%7s\n'+
                    '  %7s%7s%7s\n') % (
                self['index'] + 1,
                len(self['x']), self['rand_init'],
                int(self['x'][0]), int(self['y'][0]), int(self['z'][0]),
                int(self['x'][1]), int(self['y'][1]), int(self['z'][1]),
                int(self['x'][2]), int(self['y'][2]), int(self['z'][2]))
        except IndexError:
            return ('LAMMPS model of %s particles with: \n' +
                    ' - random initial value: %s\n' +
                    ' - first coordinates:\n'+
                    '      X    Y    Z\n'+
                    '  %5s%5s%5s\n') % (
                len(self['x']), self['rand_init'],
                self['x'][0], self['y'][0], self['z'][0])
