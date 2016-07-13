// medium
mesh_s1 = 0.75;
mesh_s2 = 0.2;
mesh_s3 = 0.015;

// coarse
mesh_s1 = 0.8;
mesh_s2 = 0.3;
mesh_s3 = 0.02;


H  = 12.0;
L  = 19.5;
l  = 4; 
h = 0.06;
d  = 4.5;

Cx  = 5;
Cy  = 6;
r   = 0.5;

Ccx = 8.5;
Ccy = 6;
dc  = 5;
hc  = 2.5;


Point(1) = {0, 0,    0, mesh_s1};
Point(2) = {L, 0,    0, mesh_s1};


Point(3) = {L, H,    0, mesh_s1};
Point(4) = {0, H,    0, mesh_s1};



Point(5) = {Cx-r, Cy-r, 0, mesh_s3};
Point(6) = {Cx+r, Cy-r, 0, mesh_s3};
Point(7) = {Cx+r, Cy+r, 0, mesh_s3};
Point(8) = {Cx-r, Cy+r, 0, mesh_s3};
Point(9) = {Cx+r, Cy+h/2, 0, mesh_s3};
Point(10) = {Cx+r, Cy-h/2, 0, mesh_s3};


Point(11) = {Cx+r+l, Cy-h/2, 0, mesh_s3};
Point(12) = {Cx+r+l, Cy+h/2, 0, mesh_s3};
Point(13) = {Cx+r+l, Cy, 0, mesh_s3};


Point(21) = {Ccx-dc, Ccy-hc,    0, mesh_s2};
Point(22) = {Ccx+dc, Ccy-hc,    0, mesh_s2};
Point(23) = {Ccx+dc, Ccy+hc,    0, mesh_s2};
Point(24) = {Ccx-dc, Ccy+hc,    0, mesh_s2};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 10};
Line(7) = {9, 7};
Line(8) = {7, 8};
Line(9) = {8, 5};


Line(10) = {10, 11};
Line(11) = {11, 13};
Line(12) = {13, 12};
Line(13) = {12, 9};
Line(14) = {9, 10};
Line Loop(15) = {4, 1, 2, 3};
Line Loop(16) = {8, 9, 5, 6, 10, 11, 12, 13, 7};


Line(21) = {21, 22};
Line(22) = {22, 23};
Line(23) = {23, 24};
Line(24) = {24, 21};

Line Loop(25) = {24, 21, 22, 23};


Plane Surface(17) = {15, 25};
Plane Surface(27) = {25, 16};

Line Loop(18) = {13, 14, 10, 11, 12};
Plane Surface(19) = {18};
