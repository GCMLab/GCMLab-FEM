//+
W = 3.2;
H = 1.23;
s = 1;


SetFactory("OpenCASCADE");
Point(1) = {0, 0, 0, s};
Point(2) = {W, 0, 0, s};
Point(3) = {W, H, 0, s};
Point(4) = {0, H, 0, s};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+//+
Transfinite Curve {2} = 3 Using Progression 1;
//+
Transfinite Curve {3} = 4 Using Progression 1;
