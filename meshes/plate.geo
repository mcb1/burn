// Gmsh project created on Thu Oct  8 16:42:19 2015
Point(1) = {0, 0, 0, 0.001};
Point(2) = {0.04, 0, 0, 0.001};
Point(3) = {0.04, 0.0085, 0, 0.0001};
Point(4) = {0, 0.0085, 0, 0.0001};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
