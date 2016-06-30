"""

solution.py
Module to define solution functions

"""

from dolfin import *
import numpy as np

def solve( data, out_path, mesh_objects, problem, z, z0 ):

    # get mesh objects
    mesh        = mesh_objects[0]
    boundary    = mesh_objects[1]
    bound_dict  = mesh_objects[2]
    num_dim     = mesh.topology().dim()

    # set up solver

    # ...form compiler options
    parameters["form_compiler"]["optimize"]         = True
    parameters["form_compiler"]["cpp_optimize"]     = True
    parameters["form_compiler"]["representation"]   = "quadrature"
#    parameters["form_compiler"]["quadrature_degree"]   = 4

    # ...create Newton solver
    solver = NewtonSolver()

    # ...set solver parameters
    prm = solver.parameters
    prm["linear_solver"]              = "lu"
    #prm["linear_solver"]              = "gmres"
    #prm["preconditioner"]             = "petsc_amg"
    prm["convergence_criterion"]      = "incremental"
    prm["relative_tolerance"]         = 1e-6
    prm["maximum_iterations"]         = 20
    prm["relaxation_parameter"]       = 1.0
    prm["absolute_tolerance"]         = 1e-9
    #prm["report"]                     = True

    # initial conditions
    z0.interpolate(Constant(data.z_i + (10.0,)))
    z.interpolate(Constant(data.z_i + (10.0,)))

    # set up output files

    # ...Paraview output
    pvd_files = [ File(out_path + "/paraview/T_data.pvd", "compressed") ]
    pvd_files += [ File(out_path + "/paraview/Y_" + str(i+1) + "_data.pvd",
                       "compressed") for i in range(0, data.N - 1)]
    pvd_files.append( File(out_path + "/paraview/phi_data.pvd", "compressed") )

    if data.ss == True:

        # solve system step
        solver.solve( problem, z.vector() )
        z_lst = z.split( True )

        # save paraview data
        for i in range(0, data.N + 1):
            pvd_files[i] << z_lst[i]

        # compute L_2 norm with exact solutions
        for soln in data.soln_e:
            E = errornorm(soln[1], z_lst[soln[0]], norm_type='L2', degree_rise=0)
            print "L_2 norm = ", E
    
    else:

        # ...name mass output *.csv output file
        mass_file    = out_path + "/mass.csv"
        pt_vals_file = out_path + "/pt_vals.csv"
    
        # ...compute initial mass data
        z_lst    = z.split( True )
        m0       = assemble( data.rho*dx( mesh ) )
    
        # ...initiate mass output *.csv file
        f_handle = open(mass_file, 'w')
        np.savetxt( f_handle, np.array([0., m0, 0.]).reshape(1,3), fmt='%.5e',
                    delimiter=",",
                    header="Time (s), Mass (kg/m^" + str(3 - num_dim)+ "), "
                            "MLR (kg/m^" + str(3 - num_dim) + "-s)",
                    comments="")
        f_handle.close()
    
        # ...initiate pt_vals output file *.csv file
        header_str = "Time (s)"
        pv_array = np.array([0.])
    
        for pv in data.pt_vals:
    
            # add initial value for pv
            pv_array = np.append( pv_array, data.z_i[pv[0]] )
    
            if pv[0] == 0:
                header_str += ", Temperature (K) at " + str(pv[1]).replace(",", "")
            else:
                header_str += ", Mass Fraction at " + str(pv[1]).replace(",", "")
    
        f_handle = open(pt_vals_file, 'w')
        np.savetxt( f_handle, pv_array.reshape(1, len(pv_array)),
                    fmt='%.5e', delimiter=",", header=header_str, comments="")
        f_handle.close()
    
        # initial state
        t = 0.0
        time_since_save = 0.0
    
        # save initial fields to Paraview files
        for i in range(0, data.N + 1):
            pvd_files[i] << (z_lst[i], t)
    
        while (t < data.t_f):
    
            # print status to terminal
            print "Time = ", t, "s"
    
            # solve system at next time step
            solver.solve( problem, z.vector() )
    
            # update time
            t += data.dt
    
            # clip mass fractions
            z.vector()[ z.vector() < 0.0 ] = 0.0
            z.vector()[ (z.vector() < 2.0)*(z.vector() > 1.0) ] = 1.0
    
            # move mesh
            #z_lst   = split( z )
            z_lst   = z.split( True )
           
            dx_s    = project(-grad( z_lst[-1] )*data.dt,
                                VectorFunctionSpace(mesh, "CG", 2),
                                solver_type="cg", preconditioner_type="ilu")
    
            #print dx_s.vector()[:].min(), dx_s.vector()[:].max()
            
            # temporary fix to prevent moving mesh for conduction problems
            if data.N > 1:
                mesh.move( dx_s )
    
            # compute current mass on current mesh
            m = assemble( data.rho*dx( mesh ) )
    
            # save data if necessary
            if (time_since_save <= (data.dt_save - 2.*data.dt)):
                time_since_save += data.dt
            else:
    
                # direct method for MLR
                mlr  = -(m - m0)/data.dt
    
                # save paraview data
                for i in range(0, data.N + 1):
                    pvd_files[i] << (z_lst[i], t)
    
                # save data to mass file
                f_handle = open(mass_file, 'a')
                np.savetxt( f_handle,
                            np.array([t, m, mlr]).reshape(1,3),
                            fmt='%.5e', delimiter=',', comments='')
                f_handle.close()
    
                # save point value data
                pv_array = np.array([t])
                for pv in data.pt_vals:
                    pv_array = np.append( pv_array, z_lst[pv[0]](pv[1]) )
    
                f_handle = open(pt_vals_file, 'a')
                np.savetxt( f_handle,
                            pv_array.reshape(1, len(pv_array)),
                            fmt='%.5e', delimiter=',', comments='')
                f_handle.close()
    
                # reset time since save
                time_since_save = 0.0
    
            # update solution
            z0.vector()[:] = z.vector()
            m0 = m
