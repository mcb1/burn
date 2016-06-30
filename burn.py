#!/usr/bin/env python

"""

burn.py
A model of burning materials

"""

import sys
import os
from time import clock

# add path to scripts
sys.path.insert(0, "/home/fenics/shared/scripts")

import meshing
import model
import solution

# read input file module from command line
input_file  = os.path.splitext(sys.argv[1])[0]
data        = __import__(input_file)

# create directory in "./output_name/" for storing results
out_path    = r"./output_" + data.name
if not os.path.exists(out_path):
    os.makedirs(out_path)

print "Burning " + data.name

# build mesh and boundary
mesh_objects = meshing.define_objects( data )

# define transport problem
[ problem, z, z0 ] = model.define_problem( data, mesh_objects )

start = clock()

# solve coupled problem in time
solution.solve( data, out_path, mesh_objects, problem, z, z0 )

elapsed = clock() - start
print "Finished. Total wall-clock time = ", elapsed, " s"
