h = 0.25;

Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {2, 0, 0, h};
Point(4) = {2, 1, 0, h};
Point(5) = {1, 1, 0, h};
Point(6) = {0, 1, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 6};
Line(4) = {6, 1};
Line(5) = {5, 4};
Line(6) = {4, 3};
Line(7) = {3, 2};

Line Loop(8) = {1, 2, 3, 4};
Plane Surface(9) = {8};

Line Loop(10) = {7, 2, 5, 6};
Plane Surface(11) = {10};

Physical Surface(11) = {9};
Physical Surface(12) = {11};
Physical Line(1) = {1, 7};
Physical Line(2) = {6};
Physical Line(3) = {5, 3};
Physical Line(4) = {4};

Recombine Surface {9, 11};
