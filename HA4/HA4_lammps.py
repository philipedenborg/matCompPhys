#!/usr/bin/env python
# -*- coding:utf-8 -*-

import lammps_wrapper as lmp
from pprint import pprint

print('Starting...')

setup = lmp.Lammmps_setup()
setup.gen_lmp_file('test.in')

lmp.run_lammps('test.in')

data = lmp.read_log('log.lammps')
pprint(data)