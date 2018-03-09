#!/usr/bin/env python
# -*- coding:utf-8 -*-

from subprocess import call
import os

class Lammmps_setup():
    ''' 
    A Lammps setup that can be used to generate a Lammps file to be run on 
    Hebbe using mpi.
    '''

    def __init__(self):

        # parameters are everything that is not a variable right now, was lazy 
        # as usual.
        self.parameters = {'dimension': '3',
                           'boundary': 'p p p',
                           'units': 'metal',
                           'atom_style': 'atomic',
                           'neighbor': '0.3 bin',
                           'neigh_modify': 'every 10 delay 0 check no',
                           'lattice': 'fcc ${a}',
                           'region': 'box block 0 ${nx} 0 ${ny} 0 ${nz}',
                           'create_box': '1 box',
                           'create_atoms': '1 box',
                           'velocity': 'all create ${T2} 1337',
                           'mass': '1 26.9815',
                           'pair_style': 'eam/alloy',
                           'pair_coeff': '* * my_al_potential.alloy Al',
                           'dump': 'id all atom 1000 atomsdump',
                           'timestep': '0.001',
                           'thermo_style': ('custom step temp etotal pe press'
                                            'pxx pyy pzz lx ly lz'),
                            'thermo': '100',
                            'fix': ('1 all npt temp ${T} ${T}'
                                    '1 iso 0 0 3 drag 1.0'),
                            'run': '${eqsteps}'}

        
        # Variables are treated specialy
        self.variables = {'eqsteps': '100000',
                          'T': '600.0',
                          'Tmelt': '2000.0',
                          'nx': '6',
                          'ny': '6',
                          'nz': '6',
                          'a': '4.05',
                          'T2': '2*${T}'}

        # The order that the file is constructed, basically which lines to  
        # write in which order.
        self.order = ['dimension', 
                      'boundary', 
                      'units', 
                      'atom_style', 
                      'neighbor', 
                      'neigh_modify', 
                      'eqsteps', 
                      'T', 'Tmelt', 
                      'nx', 'ny', 'nz', 'a', 
                      'lattice', 
                      'region', 
                      'create_box', 
                      'create_atoms', 
                      'mass', 
                      'pair_style', 
                      'pair_coeff', 
                      'T2', 
                      'velocity', 
                      'dump', 
                      'timestep', 
                      'thermo_style', 
                      'thermo', 
                      'fix', 
                      'run']

    def set_parameter(self, key, val):
        '''
        Sets a parameter with the string name 'key' to have the value 'val'
        '''
        if not isinstance(key, str):
            raise ValueError('Dict key \'{0}\' is not a string!'.format(key))

        self.parameters[key] = val

    def set_variable(self, key, val):
        '''
        Sets a variable with the string name 'key' to have the value 'val'
        '''
        if not isinstance(key, str):
            raise ValueError('Dict key \'{0}\' is not a string!'.format(key))

        self.variables[key] = val

    def set_order(self, order):
        '''
        Sets the order in which the file will be written
        '''
        self.order = order

    def param_str(self, key):
        '''
        Creates the string representation of a parameter that is a correct 
        command in Lammps.
        '''
        s = '{0} \t{1}\n'.format(key, self.parameters[key])
        return s

    def var_str(self, key):
        '''
        Creates the string representation of a variable that is a correct 
        command in Lammps.
        '''
        s = 'variable \t{0} equal {1}\n'.format(key, self.variables[key])
        return s

    def gen_lmp_file(self, filename):
        '''
        Creates a Lammps file from the setup variables.
        '''
        p = self.parameters
        v = self.variables
        with open(filename, 'w') as lmp_file:
            for action in self.order:
                if action in p:
                    lmp_file.write(self.param_str(action))
                elif action in v:
                    lmp_file.write(self.var_str(action))
                else:
                    err_msg = ('The action \'{}\' is not defined!'
                                .format(action))
                    raise ValueError(err_msg)

def run_lammps(lmp_file, actually_run=True):
    '''
    Runs a Lammps script using mpi.
    lmp_file: file path for the Lammps file.
    actually_run: debug flag, if set to false Lammps will not actually be run 
    and a string will be printed to stdout instead.
    '''
    if not os.path.isfile(lmp_file):
        raise FileNotFoundError('File \'{0}\' could not be found.'
                                .format(lmp_file))
    # Remove the comment below to actually run the program, this part has not 
    # been tested yet though...
    if actually_run:
        call(['mpirun', 'lmp_mpi', '-in', lmp_file])
    else:
        print('In debug mode, not running Lammps.')
        

def read_log(log_file):
    '''
    Returns a dictionary with the numerical data from the output log generated 
    by Lammps.
    ''' 
    data_dict = {}

    # I'm lazy, so I will just split this string to see how many vaiables we 
    # have.
    for word in 'Step Temp TotEng PotEng Press Pxx Pyy Pzz Lx Ly Lz'.split():
        data_dict[word.lower()] = []

    in_data = False # Are we in the data we want to extract?
    total_lines = 0
    data_lines = 0

    with open(log_file, 'r') as log:
        for line in log:
            total_lines += 1

            # This line preceeds the data
            if (not in_data and 
                'Step Temp TotEng PotEng Press Pxx Pyy Pzz Lx Ly Lz' in line):
                in_data = True
            # We have just reached the first line on none data after the data 
            elif (in_data and 'Loop time of ' in line):
                  in_data = False
            # If we are in the data, extract the values and put them into the 
            # array.
            elif in_data:
                data_lines += 1
                values = line.split()
                for val, key in zip(values, data_dict):
                    data_dict[key].append(float(val))
                
    print('Log was {0} lines long.'.format(total_lines))
    print('There was {0} lines containing data in the log.'.format(data_lines))

    return data_dict

    