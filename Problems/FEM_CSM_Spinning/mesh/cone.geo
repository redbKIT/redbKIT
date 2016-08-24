lc = 0.01;
lc_apex = 0.0004;

// Coarse
lc = 0.012;
lc_apex = 0.001;

// Coarse2
lc = 0.008;
lc_apex = 0.008;

// Super-Coarse
lc = 0.02;
lc_apex = 0.001;


H  = 0.1;
R  = 0.05;

Point(1) = {0,0,H,lc};
Point(2) = {R,0,H,lc};
Point(6) = {0,0,0,lc_apex};
Line(1) = {1,2};
Line(2) = {2,6};
Line(3) = {6,1};
Line Loop(4) = {2,3,1};
Plane Surface(5) = {4};
Extrude {{0,0,1}, {0,0,0}, Pi/2} { Surface{5}; }
Extrude {{0,0,1}, {0,0,0}, Pi/2} { Surface{17}; }
Extrude {{0,0,1}, {0,0,0}, Pi/2} { Surface{29}; }
Extrude {{0,0,1}, {0,0,0}, Pi/2} { Surface{41}; }

Physical Point(1) = {6};
Physical Surface(1) = {12, 16, 48, 52, 24, 40, 28, 36, 29, 5, 41, 17};
Physical Volume(2) = {2, 1, 4, 3};
