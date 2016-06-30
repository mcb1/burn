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

Burn is implemented as a set of Python modules. These modules are found in the
`scripts/` directory. A simulation is executed using the `burn.py` script along
with an input file. Input files for Burn take the form of Python modules, but
require minimal knowledge of Python to create or modify. For examples of input
file modules look in `demos/`, `validation/`, or `verification/`. 
Several meshes generate using
[Gmsh](gmesh.info) may be found in `meshes/`.

Setup and Usage
---------------

Burn requires an installation of FEniCS. The easiest way to use FEniCS is
through the FEniCS [Docker](https://www.docker.com) images. A nice discussion of
how to use FEniCS with Docker may be found
[here](http://fenics.readthedocs.io/projects/containers/en/latest/).

Once you have Docker and the FEniCS Docker script installed. It is
straightforward to start running Burn. The process is described as follows:

1. Start a Docker terminal.
2. Change directories into your Burn repository.
3. Start a FEniCS session:
```
$ fenicsproject run
```
All Burn files may be found in `~/shared/`.
4. Change directories to wherever your input file is.
5. Run Burn:
```
$ ~/shared/burn.py input.py
```

The output of the simulation will be stored locally in an automatically
generated folder. Burn writes temperature and mass fraction
fields to `vtk` files, which may be visualized using
[ParaView](www.paraview.org).

Example Usage
-------------

To run the plate demo, change directories to
```
$ cd verification/conduction/plate/
```
and then launch the simulation using
```
$ ~/shared/burn.py plate.py
```

