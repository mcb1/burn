Burn
====

Burn is a model for simulating the burning of materials. It uses tools from
the [FEniCS project](http://fenicsproject.org) to automate mesh generation,
discretization, and solution of the problem. The goal of developing Burn is to
enable predictions of burning rate and flame spread across common flammable
items. Such a capability will enable materials scientists and engineers to
design safer products.

Overview
--------

The physical model underlying Burn is expressed as a system of partial
differential equations describing conservation of energy and component mass
fraction. Diffusional heat and mass transfer are allowed, and state dependent
thermophysical properties may be specified by the user. The current
version of Burn assumes constant density, and so it is most appropriate for
non-charring and non-intumescent materials. Both Dirichlet and Neumann-type
boundary conditions may be specified by the user. One of the unique features of
Burn as compared to similar models of burning materials is that 3d material
deformation is modeled. This deformation is achieved through an arbitrary
Lagrangian-Eulerian (ALE) formulation of the governing equations.

The governing equations of Burn are discretized using linear finite elements,
and time-stepping is by the implicit Euler method. Since the discretized
equations are nonlinear, Newton's method is used to update the solution at each
time step.

Burn is implemented as a set of Python modules. These modules are found in the
`scripts/` directory. A simulation is executed using the `burn.py` script along
with an input file specified as a command line argument. Input files for Burn
take the form of Python modules, but require little to no knowledge of Python to
create or modify. For examples of input file modules look in `demos/`,
`validation/`, or `verification/`.  Several meshes generated using
[Gmsh](http://gmesh.info) may be found in `meshes/`. Gmsh `*.msh` files may be
converted to XML files as used by Burn using the `dolfin-convert.py` script found
[here](https://people.sc.fsu.edu/~jburkardt/py_src/dolfin-convert/dolfin-convert.py).

Setup and Usage
---------------

Burn requires an installation of FEniCS. The easiest way to use FEniCS is
through the FEniCS [Docker](https://www.docker.com) images. Installers for
Docker may be found
[here](https://www.docker.com/products/overview#/install_the_platform). Note
that FEniCS requires that Windows users install the [Docker
Toolbox](https://www.docker.com/products/docker-toolbox), and not
Docker for Windows. A nice discussion of
how to use FEniCS with Docker may be found
[here](http://fenics.readthedocs.io/projects/containers/en/latest/).

Once you have Docker installed, it is straightforward to start running Burn. The
process is outlined as follows:

1. Start a Docker Quickstart terminal. Your installation should create an
executable for launching such a terminal.
2. Change directories into your Burn repository. The Burn repository may be
obtained by clicking the green "Clone or download" button above.
3. Start a FEniCS session with access to your working directory by executing the
following command:  
`$ docker run -ti -v $(pwd):/home/fenics/shared quay.io/fenicsproject/stable`
4. Change directories to the location of your input file (e.g., `input.py`). An
example of this is given below.
5. Run Burn:  `$ ~/shared/burn.py input.py`

The output of the simulation will be stored locally in an automatically
generated folder. Burn writes temperature and mass fraction
fields to `vtk` files, which may be visualized using
[ParaView](http://www.paraview.org).

Example Simulation
------------------

To run the plate demo, begin by starting a Docker Quickstart Terminal. Change
directories to your Burn repository. Start a FEniCS session (see Step 3 above).
Change directories to the plate verification case:
```
$ cd verification/conduction/plate/
```
Then run the simulation using
```
$ ~/shared/burn.py plate.py
```
This example is steady state, and so it won't take long to finish. Once the
simulation is complete, start ParaView. From the File menu, click on open and
find and open
`.../burn/verification/conduction/plate/output_plate/paraview/T_data.pvd`.
Click the **Apply** button on the left of the ParaView window, and a colored
plot of the temperature field will appear.

All of the other input files in the `demos/`, `validation/`, and `verification/`
folders can be run similarly. Be aware that some of the 3d simulations can take a
while to run.
