//-*- C++ -*-

//!@file

//!This file generates a structured cubic mesh made of tetrahedra.
//!The flags that are set on the mesh are for pure traction tests
//!with symmetryc boundary conditions.
//!
//!The #Elements can be doubled easily and also the dicretization 
//!of the edges.
//!The scripts incorporates the transformation of the file *.msh into
//!the proper *.mesh file

//Mesh.CharacteristicLengthMin=0.1; useless in this case
//Mesh.CharacteristicLengthFromPoints=1; useless in this case

lc = 0.05;
R= 1;
L = 1.0;

// ############# POINTS ##############
//Base
Point(0) = {-1, -R, R, lc};
Point(1) = {-1, R, R, lc};
Point(2) = {-1, R, -R, lc};
Point(3) = {-1, -R, -R, lc};


//############### LINES ##############
//Base`
//Constrained face
Line (10) = {0,1};
Line (20) = {1,2};
Line (30) = {2,3};
Line (40) = {3,0};




// ############# LINE LOOPS ############
Line Loop(1000) = {10,20,30,40}; //fixed face

// ############# LINE LOOPS ############
Plane Surface(700) = {1000}; //on line loops

//h = 0.5;
nn = 10;
// ############# GENERATE A STRUCTURED SURFACE MESH ############
Transfinite Line{10,20,30,40}=nn;
Transfinite Surface { 700 } = {0,1,2,3} ; //on points

// ############# EXTRUSION ############
nb_layers= nn-1 ;
out[] = Extrude{2, 0, 0}{Surface{700}; Layers{nb_layers}; }; 
IndexOutlet = out[0];

// ############# VOLUME ############


// ############# EXTRUDED SURFACE ############

// ############# LATERAL SURFACE ############

//physical groups
//Physical Line(10) = {10};
//Physical Line(20) = {20};
//Physical Line(30) = {1008};
Physical Surface(6) = {-1009};//6
Physical Surface(4) = {-1013};//4
Physical Surface(5) = {-1017};//5
Physical Surface(3) = {-1021};//3
Physical Surface(2) = {-1022};//2
Physical Surface(1) = {700};//1
Physical Volume(0)={700, out[0], out[1]};


