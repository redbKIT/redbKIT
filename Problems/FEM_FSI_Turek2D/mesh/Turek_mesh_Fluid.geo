//
Merge "meshFS.msh";


//rings
Physical Point(1) = {6, 10};
Physical Point(2) = {13};//A point
Physical Point(3) = {8};//B point


Physical Line(2) = {4}; // inlet fluid
Physical Line(3) = {2}; // outlet fluid
Physical Line(4) = {1, 3}; // top/bottom fluid
Physical Line(5) = {6, 5, 7, 8}; // obstacle fixed


Physical Line(7) = {9, 11}; // top/bottom solid
Physical Line(8) = {10, 13}; // right solid

Physical Surface(1) = {15};


