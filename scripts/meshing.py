"""

meshing.py
Tools for building simple meshes using MSHR.

"""

import sys
from dolfin import *
from mshr import *

# boundary class for total boundary
class TotalBoundary( SubDomain ):
    def inside( self, x, on_boundary ):
        return on_boundary

# boundary classes for rectilinear problems
class LeftBoundary( SubDomain ):
    def inside( self, x, on_boundary ):
        return on_boundary and abs( x[0] ) < DOLFIN_EPS

class BottomBoundary( SubDomain ):
    def inside( self, x, on_boundary ):
        return on_boundary and abs( x[1] ) < DOLFIN_EPS

class RightBoundary( SubDomain):
    def __init__( self, length ):
        self.length = length
        SubDomain.__init__( self )
    def inside( self, x, on_boundary ):
        return on_boundary and abs( x[0] - self.length ) < DOLFIN_EPS

class TopBoundary( SubDomain ):
    def __init__( self, height ):
        self.height = height
        SubDomain.__init__( self )
    def inside( self, x, on_boundary ):
        return on_boundary and abs( x[1] - self.height ) < DOLFIN_EPS

# generate mesh objects
def define_objects( data ):

    # build mesh
    if data.mesh_type == "Rectangle":
        mesh = rectangle_mesh( data )
    elif data.mesh_type == "Interval":
        mesh = interval_mesh( data )
    elif data.mesh_type == "Circle":
        mesh = circle_mesh( data )
    elif data.mesh_type == "Sphere":
        mesh = sphere_mesh( data )
    elif data.mesh_type == "File":
        mesh = Mesh( data.mesh_file )
    else:
        sys.exit( "Error: Mesh type is not recognized." )

    # build and mark boundaries and boundary dictionary
    if data.mesh_type == "Rectangle" or data.bound_type == "Rectangle":
        boundary, bound_dict = rectangle_boundaries( data, mesh )
    elif data.mesh_type == "Interval":
        boundary, bound_dict = interval_boundaries( data, mesh )
    elif data.bound_type == "Total":
        boundary, bound_dict = total_boundaries( data, mesh )
    else:
        sys.exit( "Error: Boundary type is not recognized." )

    mesh_objects = (mesh, boundary, bound_dict)

    return mesh_objects

# functions for creating meshes
def rectangle_mesh( data ):

    # build mesh
    mesh = RectangleMesh( Point(0., 0.),
                          Point(data.L_x, data.L_y),
                          data.N_x, data.N_y, "crossed")

    # smooth mesh
    mesh.smooth( 50 )

    return mesh

def interval_mesh( data ):

    # build mesh
    mesh = IntervalMesh( data.N_m, 0., data.L )

    return mesh

def circle_mesh( data ):

    # build mesh
    domain = Circle( Point(0., 0.), data.R, 3*data.N_m )
    mesh = generate_mesh( domain, data.N_m )
    
    # smooth mesh
    mesh.smooth( 50 )

    return mesh

def sphere_mesh( data ):

    # build mesh
    domain = Sphere( Point(0., 0., 0.), data.R, 3*data.N_m )
    mesh = generate_mesh( domain, data.N_m )
    
    # smooth mesh
    mesh.smooth( 50 )

    return mesh

# functions for creating boundaries
def rectangle_boundaries( data, mesh ):

    # create boundary
    boundary = MeshFunction( "size_t", mesh, mesh.topology().dim() - 1 )

    # set all facets to 999 to avoid ambiguous initializations
    boundary.set_all( 999 )

    # create boundary parts
    total_boundary  = TotalBoundary()
    left_boundary   = LeftBoundary()
    bottom_boundary = BottomBoundary()
    right_boundary  = RightBoundary( data.L_x )
    top_boundary    = TopBoundary( data.L_y )

    # mark boundary parts
    left_boundary.mark( boundary, 1 )
    bottom_boundary.mark( boundary, 2 )
    right_boundary.mark( boundary, 3 )
    top_boundary.mark( boundary, 4 )

    # create dictionary of mesh boundary labels
    bound_dict = { "Left" : 1,
                   "Bottom" : 2,
                   "Right" : 3,
                   "Top" : 4 }

    return boundary, bound_dict

def interval_boundaries( data, mesh ):

    # create boundary
    boundary = MeshFunction( "size_t", mesh, mesh.topology().dim() - 1 )

    # set all facets to 999 to avoid ambiguous initializations
    boundary.set_all( 999 )

    # create boundary parts
    left_boundary   = LeftBoundary()
    right_boundary  = RightBoundary( data.L )

    # mark boundary parts
    left_boundary.mark( boundary, 1 )
    right_boundary.mark( boundary, 2 )

    # create dictionary of mesh boundary labels
    bound_dict = { "Left" : 1,
                   "Right" : 2 }

    return boundary, bound_dict

def total_boundaries( data, mesh ):

    # create boundary
    boundary = MeshFunction( "size_t", mesh, mesh.topology().dim() - 1 )

    # set all facets to 999 to avoid ambiguous initializations
    boundary.set_all( 999 )

    # create boundary parts
    total_boundary  = TotalBoundary()

    # mark boundary parts
    total_boundary.mark( boundary, 1 )

    # create dictionary of mesh boundary labels
    bound_dict = { "Total" : 1 }

    return boundary, bound_dict

