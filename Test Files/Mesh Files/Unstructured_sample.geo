//+
W = 1;
H = 1;
s = 0.5;


SetFactory("OpenCASCADE");
Point(1) = {0, 0, 0, s};
Point(2) = {1, 0, 0, s};
Point(3) = {1, 1, 0, s};
Point(4) = {0, 1, 0, s};

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
//+