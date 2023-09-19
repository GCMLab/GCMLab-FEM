//+
W = 2;
H = 2;
R = 0.5;
s = 0.5;


SetFactory("OpenCASCADE");
Point(1) = {0, 0, 0, s};
Point(2) = {R, 0, 0, s};
Point(3) = {W, 0, 0, s};
Point(4) = {W, H, 0, s};
Point(5) = {0, H, 0, s};
Point(6) = {0, R, 0, s};
//+
Circle(1) = {2, 1, 6};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Curve Loop(1) = {5, -1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Recursive Delete {
  Point{1}; 
}
//+
Transfinite Curve {4, 3, 1} = 3 Using Progression 1;
//+
Transfinite Curve {5, 2} = 3 Using Progression 1;
