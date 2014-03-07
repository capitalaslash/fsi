h = 1e20;
Lx = 1.0;
Ly = 2.0;
Vy = 0.1;
Vx = 0.1;

Point(11) = { 0, 0, 0, h};
Point(12) = {Vx, 0, 0, h};
Point(13) = {Lx, 0, 0, h};
Point(21) = { 0, 0.5*(Ly-Vy), 0, h};
Point(22) = {Vx, 0.5*(Ly-Vy), 0, h};
Point(23) = {Lx, 0.5*(Ly-Vy), 0, h};
Point(31) = { 0, 0.5*(Ly+Vy), 0, h};
Point(32) = {Vx, 0.5*(Ly+Vy), 0, h};
Point(33) = {Lx, 0.5*(Ly+Vy), 0, h};
Point(41) = { 0, Ly, 0, h};
Point(42) = {Vx, Ly, 0, h};
Point(43) = {Lx, Ly, 0, h};

Line(11) = {11,12};
Line(12) = {21,22};
Line(13) = {31,32};
Line(14) = {41,42};

Line(21) = {12,13};
Line(22) = {22,23};
Line(23) = {32,33};
Line(24) = {42,43};

Line(31) = {11,21};
Line(32) = {12,22};
Line(33) = {13,23};

Line(41) = {21,31};
Line(42) = {22,32};
Line(43) = {23,33};

Line(51) = {31,41};
Line(52) = {32,42};
Line(53) = {33,43};

Line Loop(1) = {11,32,-12,-31};
Plane Surface(1) = {1};
Line Loop(2) = {21,33,-22,-32};
Plane Surface(2) = {2};
Line Loop(3) = {12,42,-13,-41};
Plane Surface(3) = {3};
Line Loop(4) = {22,43,-23,-42};
Plane Surface(4) = {4};
Line Loop(5) = {13,52,-14,-51};
Plane Surface(5) = {5};
Line Loop(6) = {23,53,-24,-52};
Plane Surface(6) = {6};

n_Vx = 4;
n_Vy = 4;
n_Lx = 8;
n_Ly = 10;
f = 0.8;

Transfinite Line {11,12,13,14} = n_Vx+1;
Transfinite Line {21,22,23,24} = n_Lx+1;
Transfinite Line {31,32,33} = n_Ly+1 Using Progression f;
Transfinite Line {41,42,43} = n_Vy+1;
Transfinite Line {51,52,53} = n_Ly+1 Using Progression 1./f;

Transfinite Surface "*";
Recombine Surface "*";

// INLET
Physical Line(1) = {11,21};

// WALL
Physical Line(2) = {33,43,53};

// OUTLET
Physical Line(3) = {14,24};

// SYMMETRY
Physical Line(4) = {31,41,51};

// INTERFACE
Physical Line(10) = {22,23,42};

// FLUID
Physical Surface(11) = {1,2,3,5,6};

// SOLID
Physical Surface(12) = {4};
