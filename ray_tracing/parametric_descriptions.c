#include <parametric_descriptions.h>


double sign(double v)
{
  if (v > 0) return  1.0  ;
  if (v < 0) return -1.0 ;
  return 0.0 ;
}

//-----------------------------------------------------------2d functions
// circle :   x^2 + y^2 = 1
int f1 (double u, double xy[2])
{
  xy[0] = cos(u) ;  
  xy[1] = sin(u) ;
  return 1;
}


// sum4 :   x^4 + y^4 = 1
int f2 (double u, double xy[2])
{
  double c,s ;
  c = cos(u) ; s = sin(u) ;
  xy[0] = sign(c)*sqrt(fabs(c)) ;
  xy[1] = sign(s)*sqrt(fabs(s)) ;
  return 1;
}


// square :   |x| + |y| = 1
int f3 (double u, double xy[2])
{
  double c,s ;
  c = cos(u) ; s = sin(u) ;
  xy[0] = sign(c)*c*c ;
  xy[1] = sign(s)*s*s ;
  return 1;
}


// astroid :   sqrt(|x|) + sqrt(|y|) = 1
int f4 (double u, double xy[2])
{
  double c,s ;
  c = cos(u) ; s = sin(u) ;
  xy[0] = sign(c)*c*c*c*c ;
  xy[1] = sign(s)*s*s*s*s ;
  return 1;
}


// hyperbola :   x^2 - y^2 = 1 
int f5 (double u, double xy[2])
// right branch :
{
  xy[0] = cosh(u) ;
  xy[1] = sinh(u) ;
  return 1;
}


// parabola :   y = x^2
int f6 (double u, double xy[2])
{
  xy[0] = u ;
  xy[1] = u*u ;
  return 1;
}


// lemon :   x^2 - (1 - y^2)^3 = 0
int f7 (double u, double xy[2])
{
  double c ;
  c = cos(u) ;
  xy[0] = c*c*c ;
  xy[1] = sin(u) ;
  return 1;
}

int f8 (double u, double v, double xyz[3]) {
  double sv = sin(v); double su = sin(u);
  double cv = cos(v); double cu = cos(u);

  xyz[0] = cv*cu;
  xyz[1] = sv;
  xyz[2] = cv * su;
  return 1;
}
int f9 (double u, double v, double xyz[3]) {
  double hcu = cosh(u); double hsu = sinh(u);
  double cv = cos(v); double sv = sin(v);

  xyz[0] = hcu*cv;
  xyz[1] = hcu*sv;
  xyz[2] = hsu;
  return 1;
}

int sphere(double u, double v, double xyz[3]) {
  return f8(u,v,xyz);
}

int hyperboloid(double u, double v, double xyz[3]) {
  return f9(u,v,xyz);
}
