"""

plate.py
Input module for heat conduction across a 2d plate

"""

from dolfin import *

# Using material properties from Mill for coppert at 300 K
def k( z ):
    # W/m-K
    return 1.

def c_p( z ):
    # J/kg-K
    return 1.

name        = "plate"  # problem name

# problem descriptors
mesh_type   = "Rectangular"

# geometry and mesh parameters
bound_type  = "Rectangular"

N_y        = 64                    # Number of elements in y direction
N_x        = 2*N_y                    # Number of elements in x direction
L_x         = 0.01                  # length of slab in x direction, m
L_y         = 0.005                  # length of slab in y direction, m

# simulation parameters
ss          = True

# material model
N           = 1                     # number of components
N_v         = 0                     # number of volatile components
N_r         = 0                     # number of reactions

# scenario parameters (boundary conditions)
T_p     = 398.15                # peak hot surface temperature, K
T_inf   = 298.15                 # initial temperature, K
z_i     = (T_inf,)                # initial state vector

# natural boundary conditions for transport: (equation, surface, function)
bcs_tn  = []

# essential boundary conditons for transport: (equation, surface, value)
T_h = Expression('T_inf + (T_p-T_inf)*sin(pi*x[0]/L_x)', T_inf=T_inf, T_p=T_p, L_x=L_x)

bcs_te  = [ (0, "Bottom", T_h),
            (0, "Left", Constant(T_inf)),
            (0, "Top", Constant(T_inf)),
            (0, "Right", Constant(T_inf)) ]

# exact solutions for comparison: (equation, function for exact solution)
T_e = Expression('T_inf + (T_p-T_inf)*sin(pi*x[0]/L_x)*sinh(pi*(L_y-x[1])/L_x)/sinh(pi*L_y/L_x)',
                 T_inf=T_inf, T_p=T_p, L_x=L_x, L_y=L_y)
soln_e  = [ (0, T_e) ]

# reaction model
nu      = []
reactant_idx = []

# material properties
rho     = 1.         # density, kg/m^3
A       = []    # reaction pre-exponentials, 1/s
E       = []     # reaction activation energies, J/mol
dh      = []             # heats of reaction, J/kg
eps     = 1.                  # emissivity
D_v     = []             # volatile diffusivity, m^2/s
