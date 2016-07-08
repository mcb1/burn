"""

pom_50.py
Input module for the gasification of a 1d slab of POM at 70 kW/m^2.

"""

from dolfin import *

# using material properties from Jing Li and Stas Stoliarov

def k( z ):
    # W/m-K
    k_2 = 0.19 - 0.00006*z[0]
    k_3 = 0.21 + 0.00001*z[0]
    k_4 = 0.25 + 0.00002*z[0]
    Y_4 = 1. - z[2] - z[3]
    return k_2*z[2] + k_3*z[3] + k_4*Y_4

def c_p( z ):
    # J/kg-K
    c_3 = 1650. + 1.2*z[0]
    c_4 = -1860 + 9.9*z[0]
    Y_4 = 1. - z[2] - z[3]
    return c_3*(z[2] + z[3]) + c_4*Y_4

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

name        = "pom_50"          # problem name

# constants
sigma       = 5.67e-8           # Stefan-Boltzman constant, W/m^2-K^4

# geometry and mesh parameters
mesh_type   = "Interval"
bound_type  = "Interval"
N_m         = 256               # characteristic number of elements
L           = 6.0e-3            # height of slab, m

# simulation parameters
ss          = False             # steady state problem flag
dt          = 5e-2              # time step, s
dt_save     = 1.                # time interval of data saves
t_f         = 900.              # final time, s
pt_vals     = []                # point values: (equation, coordinates)

# material model
N           = 4                 # number of components
N_v         = 1                 # number of volatile components
N_r         = 3                 # number of reactions

# scenario parameters (boundary conditions)
q_ext       = (-50e3,)          # externally applied heat flux vector, W/m^2
h           = 10.               # convection coefficient, W/m^2-K
T_inf       = 298.              # ambient temperature, K
q_int       = 0.                # internal energy generation, W/m^3
z_i         = (T_inf, 0.0, 0.0, 0.0)          # initial state vector

# natural boundary conditions for transport: (equation, surface, function)
bcs_tn      = [ (0, "Right", q_top),
                (0, "Left", q_bottom) ]

# essential boundary conditons for transport: (equation, surface, value)
bcs_te      = [ (1, "Right", 0.0),
                (1, "Left", 0.0) ]

# reaction model
nu          = [ [ 0.0, 0.6, 1.0],
                [ 0.0, 0.4,-1.0],
                [ 1.0,-1.0, 0.0],
                [-1.0, 0.0, 0.0] ]  # stoichiometric coefficients
reactant_idx= [ 4, 3, 2 ]           # indices of reactant components

# material properties
rho         = 1424.             # density, kg/m^3

# # ...original kinetics from Stoliarov et al.
# A           = [2.69e42, 3.84e14, 4.76e44]    # reaction pre-exponentials, 1/s
# E           = [3.82e5, 2.00e5, 5.90e5]     # reaction activation energies, J/mol

# ..."smoothed" kinetics
A           = [1.09e27, 3.84e14, 3.55e21]   # reaction pre-exponentials, 1/s
E           = [2.50e5, 2.00e5, 3.00e5]      # reaction activation energies, J/mol

dh          = [1.92e5, 11.92e5, 13.52e5]    # heats of reaction, J/kg
eps         = 0.95                          # emissivity
D_v         = [1.0e-1,]                     # volatile diffusivity, m^2/s
