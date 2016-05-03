
// Coarse
cl1 = 0.085;
cl2 = 0.025;
cl3 = cl1;

// Medium
//cl1 = 0.06;
//cl2 = 0.0075;
//cl3 = cl1;

// Fine
//cl1 = 0.05;
//cl2 = 0.005;
//cl3 = cl1;

// Fin2
//cl1 = 0.05;
//cl2 = 0.005;
//cl3 = cl1;


z1 = 0.4;
z2 = -0.4;
z3 = 0.1;

// Middle surface
Point(1) = {-0.5, 0.4, 0, cl1};
Point(2) = {-0.5, 0, 0, cl1};

Point(4) = {0.8, 0, 0, cl1};
Point(5) = {0, 0.4, 0, cl1};
Point(6) = {0.8, 0.4, 0, cl1};
Point(7) = {0.15, 0, 0, cl2};
Point(8) = {0.15, 0.1, 0, cl2};
Point(9) = {0.15, 0.4, 0, cl1};

Point(11) = {0, 0.1, 0, cl2};

Point(13) = {-0.1, 0, 0, cl2};
Point(14) = {-0, -0, 0, cl2};






 
Line(1) = {1, 2};
Line(2) = {2, 13};
Line(3) = {13, 14};
Line(4) = {14, 7};
Line(5) = {7, 4};
Line(6) = {4, 6};
Line(7) = {6, 9};
Line(8) = {9, 5};
Line(9) = {5, 1};
Line(10) = {7, 8};
Line(11) = {8, 11};

Circle(17) = {11, 14, 13};

Line(18) = {11, 14};
Line(19) = {11, 5};
Line(20) = {9, 8};






Line Loop(21) = {1, 2, -17, 19, 9};
Plane Surface(22) = {21};
Line Loop(23) = {17, 3, -18};
Plane Surface(24) = {23};
Line Loop(25) = {18, 4, 10, 11};
Plane Surface(26) = {25};
Line Loop(27) = {19, -8, 20, 11};
Plane Surface(28) = {27};
Line Loop(29) = {7, 20, -10, 5, 6};
Plane Surface(30) = {29};
Extrude {0, 0, z3} {
  Surface{24, 26, 22, 28, 30};
}

Line Loop(261) = {71, 72, 33, 34, 74, 75};
Plane Surface(262) = {261};
Line Loop(263) = {34, 74, 99, 100, -51, -50};
Plane Surface(264) = {263};
Line Loop(265) = {120, 100, -51, 123, 124};
Plane Surface(266) = {265};


Extrude {0, 0, z1} {
  Surface{266, 262, 264};
}
Delete {
  Surface{42, 60};
}
Delete {
  Surface{60, 42};
}
Delete {
  Volume{1, 2};
}
Delete {
  Surface{26, 24};
}
Delete {
  Surface{60, 42};
}
Delete {
  Surface{46};
}


Characteristic Length {86, 72, 82} = cl3;

Delete {
  Surface{69};
}
Delete {
  Surface{47};
}


Delete {
  Surface{264};
}
Delete {
  Surface{118};
}
Delete {
  Surface{118};
}
Delete {
  Surface{264};
}
Delete {
  Surface{264};
}
Line Loop(358) = {32, 33, 34};
Plane Surface(359) = {358};
Line Loop(360) = {34, -52, -51, -50};
Plane Surface(361) = {360};
Delete {
  Surface{118};
}
Delete {
  Surface{264};
}
Delete {
  Surface{264};
}
Delete {
  Surface{264};
}
Delete {
  Volume{4, 8};
}
Delete {
  Surface{118};
}
Delete {
  Surface{264};
}
Delete {
  Surface{361};
}
Delete {
  Surface{359};
}
Plane Surface(361) = {358};
Plane Surface(362) = {360};
Surface Loop(363) = {91, 28, 109, 344, 357, 356, 362, 68, 320, 316, 284, 280, 113};
Volume(364) = {363};


Delete {
  Surface{96};
}
Delete {
  Surface{96};
}
Delete {
  Volume{3, 7};
}
Delete {
  Surface{96};
}
Delete {
  Surface{262};
}
Surface Loop(365) = {79, 22, 83, 38, 361, 312, 325, 304, 308, 324, 95, 316, 320, 91};
Volume(366) = {365};


Physical Surface(1) = {79, 304};
Physical Surface(2) = {292, 144};
Physical Surface(3) = {324, 95, 109, 344, 276, 128};
Physical Surface(4) = {308, 83, 140, 288, 356, 312};

Physical Surface(5) = {38, 361, 362, 64};
Physical Surface(6) = {68};
Physical Surface(7) = {325, 293, 357};
Physical Surface(8) = {30, 28, 22};


Physical Volume(1) = {366};
Physical Volume(2) = {364};
Physical Volume(3) = {5, 6};
Characteristic Length {90, 68} = cl1;
