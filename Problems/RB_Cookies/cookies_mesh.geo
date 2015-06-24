cl  = 0.04;

x0 = -1;
x1 = 1;
y0 = -1;
y1 = 1;
   
Point(1) = {x0, y0, 0, cl};
Point(2) = {x1, y0, 0, cl};
Point(3) = {x1, y1, 0, cl};
Point(4) = {x0, y1, 0, cl};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};


rad   = 0.2; 
f     = 0.5;

xc    = -f;
yc    = -f;
k = 4;

Point(k+1)  = {xc+rad, yc, 0, cl};
Point(k+2)  = {xc-rad, yc, 0, cl};
Point(k+3)  = {xc, yc+rad, 0, cl};
Point(k+4)  = {xc, yc-rad, 0, cl};
Point(k+5)  = {xc, yc, 0, cl};

xc    =  f;
yc    = -f;
k = k + 5;

Point(k+1)  = {xc+rad, yc, 0, cl};
Point(k+2)  = {xc-rad, yc, 0, cl};
Point(k+3)  = {xc, yc+rad, 0, cl};
Point(k+4)  = {xc, yc-rad, 0, cl};
Point(k+5)  = {xc, yc, 0, cl};


xc    = f;
yc    = f;
k = k + 5;

Point(k+1)  = {xc+rad, yc, 0, cl};
Point(k+2)  = {xc-rad, yc, 0, cl};
Point(k+3)  = {xc, yc+rad, 0, cl};
Point(k+4)  = {xc, yc-rad, 0, cl};
Point(k+5)  = {xc, yc, 0, cl};


xc    = -f;
yc    = f;
k = k + 5;

Point(k+1)  = {xc+rad, yc, 0, cl};
Point(k+2)  = {xc-rad, yc, 0, cl};
Point(k+3)  = {xc, yc+rad, 0, cl};
Point(k+4)  = {xc, yc-rad, 0, cl};
Point(k+5)  = {xc, yc, 0, cl};



Circle(5) = {7, 9, 8};
Circle(6) = {8, 9, 7};

Circle(7) = {12, 14, 13};
Circle(8) = {13, 14, 12};

Circle(9) = {17, 19, 18};
Circle(10) = {18, 19, 17};

Circle(11) = {22, 24, 23};
Circle(12) = {23, 24, 22};
Line Loop(13) = {11, 12};
Plane Surface(14) = {13};
Line Loop(15) = {9, 10};
Plane Surface(16) = {15};
Line Loop(17) = {7, 8};
Plane Surface(18) = {17};
Line Loop(19) = {6, 5};
Plane Surface(20) = {19};
Line Loop(21) = {4, 1, 2, 3};
Plane Surface(22) = {13, 15, 17, 19, 21};

Physical Line(1) = {4};
Physical Line(2) = {1};
Physical Line(3) = {2};
Physical Line(4) = {3};
Physical Surface(1) = {20};
Physical Surface(2) = {18};
Physical Surface(3) = {16};
Physical Surface(4) = {14};
Physical Surface(5) = {22};
