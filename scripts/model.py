"""

model.py
Module to define physical problem

"""

from dolfin import *

# constants
R = 8.314           # gas constant, J/K-mol

# Class for interfacing with the FEniCS Newton solver
class NonlinearEquation(NonlinearProblem):
    def __init__(self, a, L, bcs):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a
        self.bcs = bcs
    def F(self, b, x):
        assemble(self.L, tensor=b)
        [bc.apply(b, x) for bc in self.bcs]
    def J(self, A, x):
        assemble(self.a, tensor=A)
        [bc.apply(A) for bc in self.bcs]

def arrhenius_rate( A, E, T ):

    # compute temperature dependent rate, 1/s
    return A*exp(-E/(R*T))

def rhs_weak_form( data, mesh_objects, z, w ):

    # get mesh objects
    mesh        = mesh_objects[0]
    boundary    = mesh_objects[1]
    bound_dict  = mesh_objects[2]
    num_dim     = mesh.topology().dim()

    # compute mesh characteristics
    dss     = ds(domain=mesh, subdomain_data=boundary)
    normal  = FacetNormal(mesh)

    # copy problem data to local parameters
    N = data.N
    N_v = data.N_v
    N_r = data.N_r
    rho = data.rho
    D_v = data.D_v
    reactant_idx = data.reactant_idx
    A = data.A
    E = data.E
    dh = data.dh
    nu = data.nu

    # compute additional mass fractions
    Y_N     = 1. - sum(z[1:N])          # mass fraction of largest species
    Y_s     = 1. - sum(z[1:(N_v + 1)])  # mass fraction of stable species

    # compute volatile diffusion mass flux
    # ...only if N > 1
    if N > 1:

        j_v = -rho*D_v[0]*grad( z[1] )

        for i in range(1, N_v):
            i_z = i + 1     # index of equaiton in array of funciton spaces
            j_i = -rho*D_v[i]*grad( z[i_z] )
            j_v = j_v + j_i
    else:

        j_v = Constant( (0.0,)*num_dim )

    # compute source terms
    s   = [0.]*(N - 1)              # volumetric mass generation rate, kg/m^3-s
    q_c = Constant(-data.q_int)      # volumetric energy generation rate, W/m^3

    for j in range(0, N_r):

        if reactant_idx[j] == N:
            Y_r = Y_N
        else:
            Y_r = z[ reactant_idx[j] ]

        rxn_rate    = rho*Y_r*arrhenius_rate( A[j], E[j], z[0] )
        q_c         = q_c + dh[j]*rxn_rate

        for i in range(0, N - 1):
            s[i] = s[i] + nu[i][j]*rxn_rate

    # initialize list of forms for RHS functions
    F = []

    # ...(1) energy equation
    F_T = -data.k( z )*inner( grad( z[0] ), grad( w[0] ) )*dx               \
            - q_c*w[0]*dx                                                   \
            - rho*data.c_p( z )*inner( grad( z[-1] ), grad( z[0] ) )*w[0]*dx
    F.append( F_T )

    # ...(2) volatile component transport equations
    for i in range(0, N_v):
        i_z = i + 1     # index of equaiton in array of funciton spaces
        F_i = -rho*D_v[i]*inner( grad( z[i_z] ), grad( w[i_z] ))*dx         \
                   + s[i]*w[i_z]*dx                                         \
                   - rho*inner( grad( z[-1] ), grad( z[i_z] ) )*w[i_z]*dx
        F.append( F_i )

    # ...(3) stable componenet transport equations
    for i in range(N_v, N - 1):
        i_z = i + 1     # index of equaiton in array of funciton spaces
        F_i = -(z[i_z]/Y_s)*inner( j_v, grad( w[i_z] ) )*dx                 \
                   + s[i]*w[i_z]*dx                                         \
                   - rho*inner( grad( z[-1] ), grad( z[i_z] ) )*w[i_z]*dx
        F.append( F_i )

    # ...(4) potential equation
    alpha_m = 1000.*(0.2e-6)
    F_phi = -alpha_m*inner( grad( z[-1] ), grad( w[-1] ) )*dx               \
                - inner( grad( z[-1] ), grad( z[-1] ))*w[-1]*dx             \
                + (w[-1]*alpha_m/(rho*Y_s))*dot( j_v, normal )*dss
    F.append( F_phi )

    # ...apply natural boundary conditions
    for bc in data.bcs_tn:
        idx = bc[0]
        surf = bound_dict[bc[1]]
        F[idx] = F[idx] - bc[2]( z, normal )*w[idx]*dss(surf)

    return F

def define_problem( data, mesh_objects ):

    # get mesh objects
    mesh        = mesh_objects[0]
    boundary    = mesh_objects[1]
    bound_dict  = mesh_objects[2]

    # define transport problem function space
    Z_1 = FunctionSpace( mesh, "CG", 1 )    # temperature and component space
    Z_2 = FunctionSpace( mesh, "CG", 2 )    # potential function space
    Z   = MixedFunctionSpace( [Z_1]*data.N + [Z_2,] )

    # essential BCs
    bcs     = []
    for bc in data.bcs_te:
        bcs.append( DirichletBC( Z.sub(bc[0]),
                                 bc[2],
                                 boundary, bound_dict[bc[1]] ) )

    # define trial and test functions
    dz      = TrialFunction( Z )
    w_lst   = TestFunctions( Z )

    # define functions
    z       = Function( Z )
    z0      = Function( Z )

    # split mixed functions
    z_lst   = split( z )
    z0_lst  = split( z0 )

    # weak statement of equations

    # ...conservation equation RHSs
    F   = rhs_weak_form( data, mesh_objects, z_lst, w_lst )

    if data.ss == True:
       
        R = F[0]
        
        for i in range(1, data.N + 1):
            R = R + F[i]

    else:
    
        # ...build weak form residuals
    
        # ...(1) add energy equation residual
        R   = data.rho*data.c_p( z_lst )*(z_lst[0] - z0_lst[0])*w_lst[0]*dx  \
                - data.dt*F[0]
    
        # ...(2) add component scalar transport equation residuals
        for i in range(1, data.N):
            R_i = data.rho*(z_lst[i] - z0_lst[i])*w_lst[i]*dx - data.dt*F[i]
            R = R + R_i
    
        # ...(3) add potential equation residual
        R_phi = (z_lst[-1] - z0_lst[-1])*w_lst[-1]*dx - data.dt*F[-1]
        R   = R + R_phi
    
    # compute directional derivative
    a   = derivative(R, z, dz)
    
    # assemble problem matrices
    problem = NonlinearEquation(a, R, bcs)

    return (problem, z, z0)
