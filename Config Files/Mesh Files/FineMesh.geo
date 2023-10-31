// Gmsh project created on Wed Jan 18 20:31:00 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.05};
//+
Point(2) = {5, 0, 0, 0.01};
//+
Point(3) = {5, 5, 0, 0.01};
//+
Point(4) = {4, 5, 0, 0.01};
//+
Point(5) = {4, 1, 0, 0.01};
//+
Point(6) = {0, 1, 0, 0.01};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Curve Loop(1) = {5, 6, 1, 2, 3, 4};
//+
Plane Surface(1) = {1};
