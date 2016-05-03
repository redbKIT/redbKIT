cl1 = 0.05;

Lx = 2.0; // x
Ly = 2.0; // y
Lz = 0.3; // z

Point(1) = {-Lx/2, -Ly/2, -Lz/2, cl1};
Point(2) = {-Lx/2,  Ly/2, -Lz/2, cl1};
Point(3) = {Lx/2,   Ly/2, -Lz/2, cl1};
Point(4) = {Lx/2,  -Ly/2, -Lz/2, cl1};
Point(5) = {Lx/2,  -Ly/2,  Lz/2, cl1};
Point(6) = {-Lx/2, -Ly/2,  Lz/2, cl1};
Point(7) = {-Lx/2,  Ly/2,  Lz/2, cl1};
Point(8) = {Lx/2,   Ly/2,  Lz/2, cl1};

Line(1) = {6, 1};
Line(2) = {1, 4};
Line(3) = {4, 5};
Line(4) = {5, 6};
Line(5) = {7, 2};
Line(6) = {2, 3};
Line(7) = {3, 8};
Line(8) = {8, 7};
Line(9) = {6, 7};
Line(10) = {1, 2};
Line(11) = {3, 4};
Line(12) = {5, 8};

Line Loop(13) = {2, 3, 4, 1};
Plane Surface(14) = {13};
Line Loop(15) = {11, -2, 10, 6};
Plane Surface(16) = {15};
Line Loop(17) = {5, 6, 7, 8};
Plane Surface(18) = {17};
Line Loop(19) = {10, -5, -9, 1};
Plane Surface(20) = {19};
Line Loop(21) = {4, 9, -8, -12};
Plane Surface(22) = {21};
Line Loop(23) = {3, 12, -7, 11};
Plane Surface(24) = {23};
Surface Loop(25) = {22, 14, 16, 24, 18, 20};

Volume(26) = {25};

Physical Surface(1) = {14};// y = -Ly/2
Physical Surface(2) = {18};// y = Ly/2
Physical Surface(3) = {16};// z = - Lz/2
Physical Surface(4) = {20};// x = - Lx/2 -- Dirichlet
Physical Surface(5) = {24};// x = Lx/2
Physical Surface(6) = {22};// z = Lz/2 


Physical Volume(1) = {26};
