"""

copper_slab.py
Input module for heat conduction across a copper slab

"""

from dolfin import *

# Using material properties from Mill for coppert at 300 K
def k( z ):
    # W/m-K
    return 401.

def c_p( z ):
    # J/kg-K
    return 385.

# constants
sigma   = 5.67e-8     # Stefan-Boltzman constant, W/m^2-K^4

name    = "copper_slab"  # problem name

# problem descriptors
mesh_type   = "Interval"

# geometry and mesh parameters
bound_type  = "Interval"

N_m     = 64                    # Number of elements
L       = 0.04                  # length of slab in x direction, m

# simulation parameters
dt      = 1e-2                  # time step, s
dt_save = 5e-1                    # time interval of data saves
t_f     = 5.                  # final time, s

# material model
N     = 1                       # number of components
N_v   = 0                       # number of volatile components
N_r   = 0                       # number of reactions

# scenario parameters (boundary conditions)
T_h     = 398.15                # hot surface temperature, K
T_inf   = 298.15                 # initial temperature, K
z_i     = (T_inf,)                # initial state vector

# natural boundary conditions for transport: (equation, surface, function)
bcs_tn  = []

# essential boundary conditons for transport: (equation, surface, value)
bcs_te  = [ (0, "Left", T_inf),
            (0, "Right", T_h) ]

# point values to record
pt_vals = [ (0, 0.01),
            (0, 0.02),
            (0, 0.03) ]

# reaction model
nu      = []
reactant_idx = []

# material properties
rho     = 8933.         # density, kg/m^3
A       = []    # reaction pre-exponentials, 1/s
E       = []     # reaction activation energies, J/mol
dh      = []             # heats of reaction, J/kg
eps     = 0.0                  # emissivity
D_v     = []             # volatile diffusivity, m^2/s
