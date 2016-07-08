"""

pmma_test_gasification.py
Input module for gasification of test specimen of PMMA

"""

from dolfin import *

# constants
sigma   = 5.67e-8     # Stefan-Boltzman constant, W/m^2-K^4

# using material properties from Jing Li and Stas Stoliarov
def k( z ):
    # W/m-K (assumed to be average of two piecewise linear parts)
    return 0.36 - 0.00031*z[0]

def c_p( z ):
    # J/kg-K
    return 600. + 3.6*z[0]

def q_surf( z, normal ):
    # W/m^2
    q_conv  = h*(z[0] - T_inf)
    q_rad   = eps*sigma*(z[0]**4 - T_inf**4)
    q_e     = dot( Constant( q_ext ), normal )

    # only include positive values of q_e
    q_e_neg = q_e/2. - abs(q_e)/2.

    return q_conv + q_rad + q_e_neg

name        = "pmma_test_gasification"

# geometry and mesh parameters
mesh_type   = "File"
mesh_file   = "test.xml"
bound_type  = "Total"
N_m         = 50                # Number of elements
R           = 0.01              # radius of circle, m

# simulation parameters
ss          = False
dt          = 2e-2              # time step, s
dt_save     = 2e0               # time interval of data saves
t_f         = 200.              # final time, s
pt_vals     = []

# material model
N           = 2                 # number of components
N_v         = 1                 # number of volatile components
N_r         = 1                 # number of reactions

# scenario parameters (boundary conditions)
q_ext       = (0.,-70e3)        # externally applied heat flux vector, W/m^2
h           = 10.               # convective heat transfer coefficient, W/m^2-K
T_inf       = 298.              # ambient temperature, K
z_i         = (T_inf, 0.0)      # initial state vector
q_int       = 0.                # internal heat generation, W/m^3

# natural boundary conditions for transport: (equation, surface, function)
bcs_tn      = [ (0, "Total", q_surf) ]

# essential boundary conditons for transport: (equation, surface, value)
bcs_te      = [ (1, "Total", 0.0) ]
            
# reaction model
nu          = [ [1,], [-1,] ]
reactant_idx= [ 2, ]

# material properties
rho         = 1160.             # density, kg/m^3
A           = [8.60e12,]        # reaction pre-exponentials, 1/s
E           = [1.88e5,]         # reaction activation energies, J/mol
dh          = [8.46e5,]         # heats of reaction, J/kg
eps         = 0.95              # emissivity
D_v         = [1.0e-1,]         # volatile diffusivity, m^2/s
