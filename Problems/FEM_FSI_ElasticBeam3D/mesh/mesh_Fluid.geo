Merge "meshFS.msh";


Physical Surface(1) = {40};// inflow
Physical Surface(2) = {48};// outflow
Physical Surface(3) = {69, 136, 17, 27};// laterali
Physical Surface(4) = {52, 44};// top / bottom

Physical Surface(5) = {135, 131, 127, 123, 103};// obstacle
Physical Surface(6) = {107, 119, 111, 115};// FSI interface
Physical Volume(1) = {2, 1};// Fluid volume


Physical Line(12) = {102,118};// Rings

//Physical Surface(7) = {39};// beam clamped
//Physical Surface(8) = {163, 119};// beam

//Physical Volume(2) = {1};
