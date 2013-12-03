
n = 10;
h = 1e10;

Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {2, 0, 0, h};
Point(4) = {0, 1, 0, h};
Point(5) = {1, 1, 0, h};
Point(6) = {2, 1, 0, h};
Point(7) = {0.9, 0.5, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {1, 4};
//Line(4) = {2, 5};
Spline(4) = {2, 7, 5};
Line(5) = {3, 6};
Line(6) = {4, 5};
Line(7) = {5, 6};

Transfinite Line {1, 2, 3, 4, 5, 6, 7}  = n+1;

Line Loop(1) = {1, 4, -6, -3};
Plane Surface(1) = {1};

Line Loop(2) = {2, 5, -7, -4};
Plane Surface(2) = {2};

Transfinite Surface {2} = {2, 3, 6, 5};
Transfinite Surface {1} = {1, 2, 5, 4};

Recombine Surface {1, 2};

Physical Surface(11) = {1};
Physical Surface(12) = {2};
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {5};
Physical Line(4) = {7};
Physical Line(5) = {6};
Physical Line(6) = {3};
Physical Line(7) = {4};

