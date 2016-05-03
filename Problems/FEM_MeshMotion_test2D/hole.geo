// Fine
mesh_s1 = 0.015;
mesh_s2 = 0.025;
mesh_s3 = 0.005;

//Fine3
mesh_s1 = 0.01;
mesh_s2 = 0.02;
mesh_s3 = 0.001;

// Fine2
//mesh_s1 = 0.012;
//mesh_s2 = 0.022;
//mesh_s3 = 0.002;

// Coarse
mesh_s1 = 0.025;
mesh_s2 = 0.025;
mesh_s3 = 0.008;

h1 = 0.5;
h2 = 0.5;
L  = 1.0;

D  = 0.5;
r  = 0.15;


Point(1) = {0, 0,    0, mesh_s1};
Point(2) = {L, 0,    0, mesh_s1};


Point(3) = {L, h1+h2,    0, mesh_s1};
Point(4) = {0, h1+h2,    0, mesh_s1};



Point(5) = {D, D, 0, 1.0};
Point(6) = {D+r, D, 0, mesh_s3};
Point(7) = {D, D+r, 0, mesh_s3};
Point(8) = {D-r, D, 0, mesh_s3};
Point(9) = {D, D-r, 0, mesh_s3};




Circle(1) = {6, 5, 7};
Circle(2) = {7, 5, 8};
Circle(3) = {8, 5, 9};
Circle(4) = {9, 5, 6};


Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};

Line Loop(9) = {7, 8, 5, 6};
Line Loop(10) = {1, 2, 3, 4};

Plane Surface(11) = {9, 10};

Physical Line(3) = {8}; // right
Physical Line(1) = {6}; // left
Physical Line(4) = {7}; // top
Physical Line(2) = {5}; // bottom
Physical Line(5) = {1, 2, 3, 4}; // obstacle

Physical Surface(1) = {11};
