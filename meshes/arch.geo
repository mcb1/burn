r = 0.005;
R = 0.01;
h = 0.001;

Point(1) = {R,0,0,h};
Point(2) = {R+r,0,0,h};
Point(3) = {R,r,0,h};
Circle(1) = {2,1,3};
Point(4) = {R-r,0,0,h};
Circle(2) = {3,1,4};
Point(5) = {R,-r,0,h};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};

Extrude{{0,1,0},{0,0,0}, Pi/2} {
    Surface{6}; Layers{15};
}
Line Loop(29) = {10, 11, 8, 9};
Plane Surface(30) = {29};
