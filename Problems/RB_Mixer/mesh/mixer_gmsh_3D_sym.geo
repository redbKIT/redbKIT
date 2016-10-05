//  This file is part of redbKIT.
//  Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
//  Author: Federico Negri <federico.negri at epfl.ch> 

lc1      = 0.05;
lc2      = 0.025;

mu1      = 0.1;
mu2      = 0.2;

Point(1) = {-0.5, -0.2, 0, lc1};
Point(2) = {0, -0.2, 0, lc1};
Point(3) = {0, -1, 0, lc1};
Point(4) = {2-mu1/2, -1, 0, lc1};
Point(5) = {2-mu1/2, mu2, 0, lc2};
Point(6) = {2+mu1/2, mu2, 0, lc2};
Point(7) = {2+mu1/2, -1, 0, lc1};
Point(8) = {4, -1, 0, lc1};
Point(9) = {4, -0.2, 0, lc2};
Point(10) = {4.5, -0.2, 0, lc1};
Point(11) = {4.5, 0.2, 0, lc1};
Point(12) = {4, 0.2, 0, lc1};
Point(13) = {4, 1, 0, lc1};
Point(14) = {3+mu1/2, 1, 0, lc1};
Point(15) = {3+mu1/2, -mu2, 0, lc2};
Point(16) = {3-mu1/2, -mu2, 0, lc2};
Point(17) = {3-mu1/2, 1, 0, lc1};
Point(18) = {1+mu1/2, 1, 0, lc1};
Point(19) = {1+mu1/2, -mu2, 0, lc2};
Point(20) = {1-mu1/2, -mu2, 0, lc2};
Point(21) = {1-mu1/2, 1, 0, lc1};
Point(22) = { 0,1, 0, lc1};
Point(23) = {0, 0.2, 0, lc1};
Point(24) = {-0.5,0.2, 0, lc1};


Line(1) = {24, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 9};
Line(10) = {9, 10};
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 13};
Line(14) = {13, 14};
Line(15) = {14, 15};
Line(16) = {15, 16};
Line(17) = {16, 17};
Line(18) = {17, 18};
Line(19) = {18, 19};
Line(20) = {19, 20};
Line(21) = {20, 21};
Line(22) = {21, 22};
Line(23) = {22, 23};
Line(24) = {23, 24};
Line Loop(25) = {23, 24, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22};
Plane Surface(26) = {25};


Extrude {0, 0, 0.3} {
  Surface{26};
}
Surface Loop(149) = {87, 26, 55, 59, 63, 67, 71, 75, 79, 83, 148, 91, 95, 99, 103, 107, 111, 115, 119, 123, 127, 131, 135, 139, 143, 147};
Volume(150) = {149};


Physical Surface(1) = {63};// inlet
Physical Surface(2) = {59, 55, 147, 143, 135, 139, 131, 119, 127, 123, 115, 111, 107, 99, 95, 91, 87, 83, 79, 26, 75, 71, 67};// wall
Physical Surface(4) = {148};//simmetry
Physical Surface(3) = {103};//outflow
Physical Volume(1)  = {1};
