"""

hips_30.py
Input module for the gasification of a 1d slab of HIPS at 30 kW/m^2

"""

from dolfin import *

# using material properties from Jing Li and Stas Stoliarov
def k( z ):
    # W/m-K
    return 0.10 - 0.0001*z[0]

def c_p( z ):
    # J/kg-K
    return 590. + 3.4*z[0]

def q_top( z, normal ):
    # W/m^2
    q_conv  = h*(z[0] - T_inf)
    q_rad   = eps*sigma*(z[0]**4 - T_inf**4)
    q_e     = dot( Constant( q_ext ), normal )

    # only include positive values of q_e
    q_e_neg = q_e/2. - abs(q_e)/2.

    return q_conv + q_rad + q_e_neg

def q_bottom( z, normal ):
    # W/m^2
    q_conv  = h*(z[0] - T_inf)
    q_rad   = eps*sigma*(z[0]**4 - T_inf**4)

    return q_conv + q_rad

name        = "hips_30"         # problem name

# constants
sigma       = 5.67e-8           # Stefan-Boltzman constant, W/m^2-K^4

# geometry and mesh parameters
mesh_type   = "Interval"        # specify a 1d mesh
bound_type  = "Interval"
N_m         = 256               # characteristic number of elements
L           = 6.0e-3            # height of slab, m

# simulation parameters
ss          = False             # steady state problem flag
dt          = 5e-2              # time step, s
dt_save     = 1.                # time interval of data saves
t_f         = 1800.             # final time, s
pt_vals     = []                # point values: (equation, coordinates)

# material model
N           = 2                 # number of components
N_v         = 1                 # number of volatile components
N_r         = 1                 # number of reactions

# scenario parameters (boundary conditions)
q_ext       = (-30e3,)          # externally applied heat flux vector, W/m^2
h           = 10.               # convection coefficient, W/m^2-K
T_inf       = 298.              # ambient temperature, K
q_int       = 0.                # internal energy generation, W/m^3
z_i         = (T_inf, 0.0)      # initial state vector

# natural boundary conditions for transport: (equation, surface, function)
bcs_tn      = [ (0, "Right", q_top),
                (0, "Left", q_bottom)]

# essential boundary conditons for transport: (equation, surface, value)
bcs_te      = [ (1, "Right", 0.0),
                (1, "Left", 0.0) ]

# reaction model
nu          = [ [1,], [-1,] ]   # stoichiometric coefficients
reactant_idx= [ 2, ]            # indices of reactant components

# material properties
rho         = 1060.             # density, kg/m^3
A           = [1.70e20,]        # reaction pre-exponentials, 1/s
E           = [3.01e5,]         # reaction activation energies, J/mol
dh          = [6.89e5,]         # heats of reaction, J/kg
eps         = 0.95              # emissivity
D_v         = [1.0e-1,]         # volatile diffusivity, m^2/s
