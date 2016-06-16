
// medium
//cl1  = 0.03;
//cl2  = 0.007;

// fine
//cl1  = 0.02;
//cl2  = 0.005;

// Coarse
cl1  = 0.05;
cl2  = 0.01;

// fine2
cl1  = 0.015;
cl2  = 0.003;

R     = 1;  
a     = 0.05;
b     = 0.3;
l     = 0.5;
d     = 0.5;
thick = 0.02;
  
Point(1) = {-R, 0, 0, cl2};
Point(2) = { R, 0, 0, cl1};
Point(3) = { 0, R, 0, cl1};

Point(4) = { -R, a, 0, cl2};
Point(5) = { -R, a+thick, 0, cl2};

Point(6) = { -R+d, a, 0, cl2};
Point(7) = { -R+d+l, a+b, 0, cl2};
Point(8) = { -R+d+l- 0.5145*thick, a+b+ 0.8575*thick, 0, cl2};

Point(9) = { -R+d, a+thick, 0, cl2};

Point(10) = { 0, 0, 0, cl1};


Line(1) = {1, 2};
Line(2) = {1, 4};
Line(3) = {4, 6};
Line(4) = {6, 7};
Line(5) = {7, 8};
Line(6) = {8, 9};
Line(7) = {9, 5};
Circle(8) = {2, 10, 3};
Circle(9) = {3, 10, 5};
Line Loop(10) = {9, -7, -6, -5, -4, -3, -2, 1, 8};
Plane Surface(11) = {10};


Physical Line(1) = {1};// symmetry
Physical Line(2) = {2};// in
Physical Line(3) = {9, 8};// out 
Physical Line(4) = {7, 6, 5, 4, 3}; // channel

Physical Surface(1) = {11};
