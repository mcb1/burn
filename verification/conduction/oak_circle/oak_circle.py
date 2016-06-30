"""

oak_circle.py
Input module for heat conduction into an oak circle

"""

from dolfin import *

# Using material properties from Mill for coppert at 300 K
def k( z ):
    # W/m-K
    return 0.21

def c_p( z ):
    # J/kg-K
    return 2400.

def q_surf( z, normal ):
    # W/m^2
    return -q_e

name    = "oak_circle"  # problem name

# geometry and mesh parameters
mesh_type   = "Circle"          # circular mesh geometry
bound_type  = "Circle"          # boundary type for circular geometry
N_m         = 32                # number of elements
R           = 0.01              # radius of circle, m

# simulation parameters
ss          = False             # flag for steady-state problems
dt          = 1e-1              # time step, s
dt_save     = 1e0               # time interval of data saves
t_f         = 120.              # final time, s

# material model
N     = 1                       # number of components
N_v   = 0                       # number of volatile components
N_r   = 0                       # number of reactions

# scenario parameters (boundary conditions)
q_e     = 6e3                   # external heat flux, W/m^2
T_inf   = 298.15                 # initial temperature, K
z_i     = (T_inf,)                # initial state vector

# natural boundary conditions for transport: (equation, surface, function)
bcs_tn  = []
bcs_tn  = [ (0, "Side", q_surf) ]

# essential boundary conditons for transport: (equation, surface, value)
bcs_te  = [ ]

# point values to record
pt_vals = [ (0, (0.,0.25*R)),
            (0, (0.,0.50*R)),
            (0, (0.,0.75*R)) ]

# reaction model
nu      = []
reactant_idx = []

# material properties
rho     = 820.         # density, kg/m^3
A       = []    # reaction pre-exponentials, 1/s
E       = []     # reaction activation energies, J/mol
dh      = []             # heats of reaction, J/kg
eps     = 0.0                  # emissivity
D_v     = []             # volatile diffusivity, m^2/s
