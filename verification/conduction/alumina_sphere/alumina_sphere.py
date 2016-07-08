"""

alumina_sphere.py
Input module for heat conduction in an internally heated alumina sphere

"""

from dolfin import *

# Using material properties from Mill for coppert at 300 K
def k( z ):
    # W/m-K
    return 36.

def c_p( z ):
    # J/kg-K
    return 765.

name        = "alumina_sphere"  # problem name

# geometry and mesh parameters
mesh_type   = "File"
mesh_file   = "sphere.xml"
bound_type  = "Total"

# geometry and mesh parameters
R           = 0.05              # radius of circle, m

# simulation parameters
ss          = False             # steady state problem
dt          = 1e-1              # time step, s
dt_save     = 5e-1              # time interval of data saves
t_f         = 60.               # final time, s
pt_vals     = [ (0, (0.,0.,0.25*R)),    # point values to record
                (0, (0.,0.,0.50*R)),
                (0, (0.,0.,0.75*R)) ]

# material model
N           = 1                 # number of components
N_v         = 0                 # number of volatile components
N_r         = 0                 # number of reactions

# scenario parameters (boundary conditions)
q_int       = 1e7               # internal heat generation, W/m^3
T_inf       = 298.15            # initial temperature, K
z_i         = (T_inf,)          # initial state vector

# natural boundary conditions for transport: (equation, surface, function)
bcs_tn      = []

# essential boundary conditons for transport: (equation, surface, value)
bcs_te      = [ (0, "Total", T_inf) ]

# reaction model
nu          = []                # stoichiometric coefficients
reactant_idx= []                # index of reactants

# material properties
rho         = 3970.             # density, kg/m^3
A           = []                # reaction pre-exponentials, 1/s
E           = []                # reaction activation energies, J/mol
dh          = []                # heats of reaction, J/kg
D_v         = []                # volatile diffusivity, m^2/s
