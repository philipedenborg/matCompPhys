#!/usr/bin/env python
# -*- coding:utf-8 -*-

from subprocess import call
import os

class Lammps():
    ''' 
    A Lammps setup that can be used to generate a Lammps file to be run on 
    Hebbe using mpi.
    '''
    def __init__(self):


        self.dimension = None
        self.boundary = None
        self.units = None

        self.atom_style = None
        self.neighbor = None
        self.neigh_modify = None
        
        self.lattice = None

        self.region = None
        
        self.mass = None
        
        self.pair_style = None
        self.pair_coeff = None
        
        self.velocity = None
        
        self.timestep = None
        self.thermo_style = None
        self.thermo = None
        
        self.groups = []

        self.commands = []


    def set_boundary(self, bc):
        self.boundary = bc
        cmd_str = 'boundary \t{0} {1} {2}'.format(bc[0], bc[1], bc[2])
        self.commands.append(cmd_str)

    def set_dimension(self, dim):
        self.dimension = dim
        cmd_str = 'dimension \t{}'.format(dim)
        self.commands.append(cmd_str)

    def set_units(self, u):
        self.units = u
        cmd_str = 'units \t{0}'.format(u)
        self.commands.append(cmd_str)

    def set_atom_style(self, style):
        self.atom_style = style
        cmd_str = 'atom_style \t{0}'.format(style)
        self.commands.append(cmd_str)

    def set_neighbor(self, skin, style):
        self.neighbor = (skin, style)
        cmd_str = 'neighbor \t{0} {1}'.format(skin, style)
        self.commands.append(cmd_str)

    def set_neigh_modify(self, keyval):
        self.neigh_modify = keyval
        cmd_str = 'neigh_modify \t' + ' '.join(['{0} {1}'.format(key, val) 
                                           for key, val in keyval.items()])
        self.commands.append(cmd_str)

    def set_lattice(self, style, scale):
        self.lattice = (style, scale)
        cmd_str = 'lattice \t{} {}'.format(style, scale)
        self.commands.append(cmd_str)

    def set_mass(self, atom_type, m):
        self.mass = (atom_type, m)
        cmd_str = 'mass \t{0} {1}'.format(atom_type, m)
        self.commands.append(cmd_str)

    def set_velocity(self, group_id, style, args='', keywords='', 
                     value=''):
        self.velocity = (group_id, style, args, keywords, value)
        args_str = ' '.join([str(a) for a in args])
        cmd_str = 'velocity \t{0} {1} {2} {3} {4}'.format(group_id, style, args_str, keywords, value)
        self.commands.append(cmd_str)

    def set_timestep(self, ts):
        self.timestep = ts
        cmd_str = 'timestep \t{0}'.format(ts)
        self.commands.append(cmd_str)

    def set_thermo_style(self, 
            style='custom step temp etotal pe press pxx pyy pzz lx ly lz'):
        self.thermo_style = style
        cmd_str = 'thermo_style \t{0}'.format(style)
        self.commands.append(cmd_str)

    def set_thermo(self, th):
        self.thermo = th
        cmd_str = 'thermo \t{0}'.format(th)
        self.commands.append(cmd_str)

    def set_region(self, region_id, style='block', 
                   dim=[(0,3), (0,3), (0,3)]):
        self.region = [region_id, style, dim]
        dims = dim[0] + dim[1] + dim[2]
        dim_str = ' '.join([str(a) for a in dims])
        cmd_str = 'region \t{0} {1} {2}'.format(region_id, style, dim_str)
        self.commands.append(cmd_str)

    def create_box(self, N, box_id):
        if not self.region:
            raise ValueError('Needs region!!!!!!')
            
        cmd_str = 'create_box \t{} {}'.format(N, box_id)
        self.commands.append(cmd_str)

    def create_atoms(self, atom_type, style):
        if not self.region:
            raise ValueError('Needs region!!!!!!')
            
        cmd_str = 'create_atoms \t{} {}'.format(atom_type, style)
        self.commands.append(cmd_str)

    def set_pair_style(self, ps):
        self.pair_style = ps
        cmd_str = 'pair_style \t{}'.format(ps)
        self.commands.append(cmd_str)

    def set_pair_coeff(self, pc):
        self.pair_coeff = pc
        cmd_str = 'pair_coeff \t * * {}'.format(pc)
        self.commands.append(cmd_str)

    def log(self, log_file='log.lammps'):
        self.log_name(log_file)
        cmd_str = 'log \t{}'.format(log_file)
        self.commands.append(cmd_str)

    def dump(self, dump_str='id all atom 1000 atomsdump'):
        cmd_str = 'dump \t{0}'.format(dump_str)
        self.commands.append(cmd_str)

    def fix(self, fix_id, group_id, style, args):
        args_str = ' '.join([str(a) for a in args])
        cmd_str = 'fix \t{0} {1} {2} {3}'.format(fix_id, group_id, style, args_str)
        self.commands.append(cmd_str)
    
    def unfix(self, fix_id):
        cmd_str = 'unfix \t{0}'.format(fix_id)
        self.commands.append(cmd_str)

    def run(self, steps):
        cmd_str = 'run \t{0}'.format(steps)
        self.commands.append(cmd_str)

    def group(self, group_id, style, args):

        self.groups.append((group_id, style, args))
        if isinstance(args, (list, tuple)):
            args_str = ' '.join([str(a) for a in args])
        elif isinstance(args, str):
            args_str = args
        cmd_str = 'group \t{0} {1} {2}'.format(group_id, style, args_str)
        self.commands.append(cmd_str)

    def gen_from_lmp_file(self, lmp_file):
        file_commands = []

        with open(self, lmp_file) as lmp:
            for line in lmp:
                cmds = line.split()
                if '#' in cmds[0]:
                    continue
                else:
                    cmd_str = cmds[0]
                    for cmd in cmds[1:]:
                        if '#' in cmd:
                            break
                        else:
                            cmd_str += cmd

                file_commands.append(cmd_str)

        self.commands = file_commands

    def gen_lmp_file(self, filename, align=True):
        '''
        Creates a Lammps file from the setup variables.
        '''

        # Align all rows neatly
        left_str = []
        left_align = 0
        if align:
            for line in self.commands:
                left_str.append(line.split()[0])
            left_align = max([len(s) for s in left_str])

        with open(filename, 'w') as lmp_file:
            for cmd in self.commands:
                if align:
                    cmd_split = cmd.split()
                    left_cmd = cmd_split[0]
                    padding = left_align - len(left_cmd)
                    left_cmd += ' '*padding
                    cmd_split[0] = left_cmd
                    cmd_str = ' '.join(cmd_split)
                    lmp_file.write(cmd_str + '\n')
                else:
                    lmp_file.write(cmd + '\n')

def run_lammps(lmp_file, actually_run=True):
    '''
    Runs a Lammps script using mpi.
    lmp_file: file path for the Lammps file.
    actually_run: debug flag, if set to false Lammps will not actually be run 
    and a string will be printed to stdout instead.
    '''
    if not os.path.isfile(lmp_file):
        raise FileError('File \'{0}\' could not be found.'
                                .format(lmp_file))
    # Remove the comment below to actually run the program, this part has not 
    # been tested yet though...
    if actually_run:
        call(['mpirun', './lmp_mpi', '-in', lmp_file])
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
    data_order = []
    for word in 'Step Temp TotEng PotEng Press Pxx Pyy Pzz Lx Ly Lz'.split():
        data_dict[word.lower()] = []
        data_order.append(word.lower())

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
                for val, key in zip(values, data_order):
                    data_dict[key].append(float(val))
                
    print('Log was {0} lines long.'.format(total_lines))
    print('There was {0} lines containing data in the log.'.format(data_lines))

    return data_dict

    
