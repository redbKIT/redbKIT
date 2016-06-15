// See this website for the geometry: 
// http://www.featflow.de/en/benchmarks/cfdbenchmarking/flow/dfg_flow3d/dfg_flow3d_configuration.html
// Parameters for the mesh-size

Mesh.Algorithm3D = 4;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

//L1
lc     = 0.05; 
lc_cyl = 0.016;

//L2
lc     = 0.035; 
lc_cyl = 0.0125;

//L3
lc     = 0.025;
lc_cyl = 0.008;

//L4
lc     = 0.0125; 
lc_cyl = 0.004;

//L0
lc     = 0.085; 
lc_cyl = 0.025;

// Parameters to define the geometry

H  = 0.41;
L  = 2.5;
R  = 0.05;
XC = 0.45 + R;
YC = 0.15 + R;

// Points at z = 0

    // Lines rectangular domain
	
	Point(1) = {0,0,0,lc};
	Point(2) = {L,0,0,lc};
	Point(3) = {L,H,0,lc};
	Point(4) = {0,H,0,lc};
	
	Line(1) = {1,2};
	Line(2) = {2,3};
	Line(3) = {3,4};
	Line(4) = {4,1};

    // Lines Circle

	Point(5) = {XC,YC,0,lc_cyl};
	Point(6) = {XC,YC+R,0,lc_cyl};
	Point(7) = {XC+R,YC,0,lc_cyl};
	Point(8) = {XC,YC-R,0,lc_cyl};
	Point(9) = {XC-R,YC,0,lc_cyl};

	Circle(5) = {6,5,7};
	Circle(6) = {7,5,8};
	Circle(7) = {8,5,9};
	Circle(8) = {9,5,6};

    // Line Loops

	Line Loop(1) = {1,2,3,4,5,6,7,8};
	Plane Surface(1) = {1};

    // Extrude the domain in direction z of a quantity H

	out[] = Extrude{0, 0, H}{Surface{1}; }; 
 
    // Physical groups for the flags

	Physical Surface(2)  = {33}; 		// Inflow (v = v_imposed)
	Physical Surface(3)  = {25}; 		// Outflow (do nothing)
	Physical Surface(10) = {-1,21,50,29}; 	// Lateral walls (v = 0)
	Physical Surface(6)  = {37,41,45,49}; 	// Cylinder (v = 0)
	
	Physical Volume(1)   = {1};
	


