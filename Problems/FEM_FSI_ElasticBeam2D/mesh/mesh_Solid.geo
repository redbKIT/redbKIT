//
Merge "meshFS.msh";


//Physical Line(1) = {4};// Inlet Fluid
//Physical Line(2) = {2};// Outlet Fluid
//Physical Line(3) = {1, 3};// lateral Fluid
//Physical Line(4) = {8, 9, 5, 6, 7};// Obstacle
Physical Line(5) = {13, 11, 10, 12}; // Beam
Physical Line(6) = {14}; // Obstacle-beam interface

Physical Surface(2) = {19};// Solid
//Physical Surface(1) = {17,27};// Fluid

Physical Point(1) = {9, 10}; // rings
