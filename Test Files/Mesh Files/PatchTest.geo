//+
W = 1.45;
H = 2.30;
s = 0.25;


SetFactory("OpenCASCADE");
Point(1) = {0, 0, 0, s*1.5};
Point(2) = {W, 0, 0, s*1};
Point(3) = {W, H, 0, s*0.5};
Point(4) = {0, H, 0, s*2};

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
Transfinite Curve {3, 2} = 6 Using Progression 1;
