
#ifndef parametric_descriptions_stuff
#define parametric_descriptions_stuff

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int f1(double u, double xy[2]);
// describes a circle

int f2(double u, double xy[2]);
// describes rounded square x^4 + y^4 = 1

int f3(double u, double xy[2]);
// describes a square: |x| + |y| = 1

int f4(double u, double xy[2]);
// describes an asteroid: sqrt(|x|) + sqrt(|y|) = 1

int f5(double u, double xy[2]);
// describes a hyberbola: x^2 - y^2 = 1

int f6(double u, double xy[2]);
// describes a parabola: y = x^2

int f7(double u, double xy[2]);
// describes a lemon: x^2 - (1-y^2)^3 = 0


int f8(double u, double v, double xyz[3]);
// describes a sphere

int f9(double u, double v, double xyz[3]);
// describes a hyperboloid

int sphere(double u, double v, double xyz[3]);
// wrapper for f8

int hyperboloid(double u, double v, double xyz[3]);
// wrapper for f9

int torus(double u, double v, double xyz[3]); //might probs just give full donut, delete first two params
// x = cos(u)
// y = cos(v) * ( sin(u) + 4) <-- 4 is arbitrary, this is the central r
// z = sin(v) * ( sin(u) + 4)

#endif
