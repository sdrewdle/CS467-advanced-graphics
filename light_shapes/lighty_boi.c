#include <FPT.h>
#include <D3d_matrix.h>
double half_angle = 30 *M_PI/180;
double hither = .1, yon = 10, view_window[6][4];
double ambient = 0.2, diffuse_max = 0.5, spec_pow = 30;
double zbuff[800][800];
double eye[3], ls[3];

void clear_zbuff() {
  int i,j;
  for(int i=0;i<800;i++){for(int j=0;j<800;j++) {zbuff[i][j]=pow(10,50);}}
}

//-----------------------------------------------------------universal tools
void normalize(double *v, int l) {
  int i;
  double length = 0;
  for(int i=0; i<l;i++) {
    length += v[i]*v[i];
  } length = sqrt(length);
  for(int i=0;i<l;i++) {
    v[i] /= length;
  }
}

void scale(double *v, double scale, int l) {
  int i;
  for(i=0;i<l;i++) {
    v[i] *= scale;
  }
}

double sgn(double v)
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
}


// sum4 :   x^4 + y^4 = 1
int f2 (double u, double xy[2])
{
  double c,s ;
  c = cos(u) ; s = sin(u) ;
  xy[0] = sgn(c)*sqrt(fabs(c)) ;
  xy[1] = sgn(s)*sqrt(fabs(s)) ;
}


// square :   |x| + |y| = 1
int f3 (double u, double xy[2])
{
  double c,s ;
  c = cos(u) ; s = sin(u) ;
  xy[0] = sgn(c)*c*c ;
  xy[1] = sgn(s)*s*s ;
}


// astroid :   sqrt(|x|) + sqrt(|y|) = 1
int f4 (double u, double xy[2])
{
  double c,s ;
  c = cos(u) ; s = sin(u) ;
  xy[0] = sgn(c)*c*c*c*c ;
  xy[1] = sgn(s)*s*s*s*s ;
}


// hyperbola :   x^2 - y^2 = 1 
int f5 (double u, double xy[2])
// right branch :
{
  xy[0] = cosh(u) ;
  xy[1] = sinh(u) ;
}


// parabola :   y = x^2
int f6 (double u, double xy[2])
{
  xy[0] = u ;
  xy[1] = u*u ;
}


// lemon :   x^2 - (1 - y^2)^3 = 0
int f7 (double u, double xy[2])
{
  double c ;
  c = cos(u) ;
  xy[0] = c*c*c ;
  xy[1] = sin(u) ;
}

//-----------------------------------------------------------3d tools
void vector_to_pts(double x1, double y1, double z1,
	       double x2, double y2, double z2,
	       double v[3]){
  //create a vector from p2 to p1
  v[2] = z2-z1;
  v[1] = y2-y1;
  v[0] = x2-x1;
}
double dot(double a[3], double b[3]) {
  int i;
  double prod = 0;
  for(i=0;i<3;i++) {
    prod += a[i] * b[i];
  }
  return prod;
}
void vector_to(double a[3], double b[3], double v[3]) {
  int i;
  for(i=0;i<3;i++) {v[i] = a[i]-b[i];}
}
void orthogonal(double u[3], double v[3], double orthog[3]){
  double t[3];
  D3d_x_product(t,u,v);

  orthog[0] = 0 - t[0];
  orthog[1] = 0 - t[1];
  orthog[2] = 0 - t[2];
}
void normal_pts(double x1, double y1, double z1,
	    double x2, double y2, double z2,
	    double x3, double y3, double z3, double n[3]) {
  double u[3], v[3];
  vector_to_pts(x1,y1,z1, x2,y2,z2, u);
  vector_to_pts(x3,y3,z3, x2,y2,z2, v);
  orthogonal(u,v, n);
  normalize(n,3);
}
void normal(double a[3], double b[3], double c[3], double n[3]) {
  double u[3], v[3];
  vector_to(a,b,u);
  vector_to(c,b,v);
  orthogonal(u,v,n);
  normalize(n,3);
}
void reflection(double l[3], double n[3], double r[3]) {
  int i;
  for(i=0;i<3;i++) {
    r[i] = 2*dot(n,l)*n[i]/dot(n,n)-l[i];
  }
  normalize(r,3);
}
void build_view_window() {
  int i;
  double x=0,y,z, norm[3];
  //hither plane
  view_window[0][0] = 0; view_window[0][1] = 0;
  view_window[0][2] = hither+1;
  view_window[0][3] = -view_window[0][2] * hither;

  //yonder plane
  view_window[1][0] = 0; view_window[1][1] = 0;
  view_window[1][2] = yon+1;
  view_window[1][3] = -view_window[1][2] * yon;

  y=tan(half_angle)/hither; z = hither;
  normal_pts(1,0,0, 0,0,0, x,y,z, view_window[2]);
  view_window[2][3] = 
    - view_window[2][0] * x
    - view_window[2][1] * y
    - view_window[2][2] * z;
  
  y = -y;
  normal_pts(1,0,0, 0,0,0, x,y,z, view_window[3]);
  view_window[3][3] = 
    - view_window[3][0] * x
    - view_window[3][1] * y
    - view_window[3][2] * z;

  x = -y; y=0;
  normal_pts(0,1,0, 0,0,0, x,y,z, view_window[4]);
  view_window[4][3] = 
    - view_window[4][0] * x
    - view_window[4][1] * y
    - view_window[4][2] * z;

  x = -x;
  normal_pts(0,1,0, 0,0,0, x,y,z, view_window[5]);
  view_window[5][3] = 
    - view_window[5][0] * x
    - view_window[5][1] * y
    - view_window[5][2] * z;
}

double add(double xyz[3], int i){
  return view_window[i][0] * xyz[0] +
    view_window[i][1] * xyz[1] +
    view_window[i][2] * xyz[2] +
    view_window[i][3];
}
int in_view(double xyz[3]) {
  double checkp[3] = {0,0,hither+yon/2};
  
  double p,check;
  int i;
  for (i=0; i<6; i++) { //check edge planes
    p = add(xyz,i);
    check = add(checkp,i);
    if(sgn(p) != sgn (check)) {return 0;}
    }
  return 1; //point is in view
}

//-----------------------------------------------------------3d functions
// sphere :  x^2 + y^2 + z^2 = 0
int f8 (double u, double v, double xyz[3]) {
  double sv = sin(v); double su = sin(u);
  double cv = cos(v); double cu = cos(u);

  xyz[0] = cv*cu;
  xyz[1] = sv;
  xyz[2] = cv * su;
}

void shade(double rgb[3], double sum, double intensity) {
  double ratio;
  double r = rgb[0], g = rgb[1], b = rgb[2];
  
  if (intensity <= sum) {
    ratio = intensity/sum;
    G_rgb(r * ratio,
	  g * ratio,
	  b * ratio);
  } else {
    ratio = (intensity-sum)/(1-sum);
    G_rgb((1-ratio)*r + ratio,
	  (1-ratio)*g + ratio,
	  (1-ratio)*b + ratio);
  }
}

/* sets the rgb value to the appropriate color
 * u,v are the current angles
 * rgb is an array of the base color
 * func is the same function used to get the point;
 *** used to approx the norm
 */
void light(double u, double v, double mat[4][4],
	   double rgb[3],
	   int (*func)(double uu, double vv, double xyz[3]))
{
  int i;

  //-------------------------------------------------define the normal vector
  double norm[3];
  double a[3], b[3], c[3];
  func(u,v,b);         D3d_mat_mult_pt(b,mat,b); //point on vertex
  func(u+0.01,v,a);   D3d_mat_mult_pt(a,mat,a);
  func(u,v+0.01,c);   D3d_mat_mult_pt(c,mat,c);
  normal(a,b,c,norm);

  //------------------------define the eye vector and the light source vector
  double e[3], l[3];
  vector_to(eye,b,e); vector_to(ls,b,l);
  normalize(e,3);
  normalize(l,3);
  //------if the ls and eye vectors are obtuse to the normal, flip the normal
  if(dot(e,norm)<=0 && dot(l,norm)<=0) {
    for (i=0;i<3;i++) {norm[i] = -norm[i];}
  }
  //--------------------------------------------------calculate the intensity
  //constants above: diffuse_max, spec_pow, ambient
  double specular, diffuse;
  
  double r[3];
  reflection(l,norm,r);
  
  double edr = dot(e,r);
  double ndl = dot(norm,l);
  double nde = dot(norm,e);

  if(ndl*nde <= 0) {edr=0; ndl=0;}

  if(edr <= 0) {
    specular = 0;
  } else {
    specular = (1-ambient-diffuse_max) * pow(edr, spec_pow);
  }

  if(ndl <= 0) {
    diffuse = 0;
  } else {
    diffuse = diffuse_max * ndl;
  }

  //calculate correct rgb value and call G_rgb to set
  shade(rgb,ambient+diffuse_max, ambient + diffuse + specular);
}

int plot_3d (double ulo, double uhi, double vlo, double vhi,
	     int (*func)(double u, double v, double xyz[3]),
	     double mat[4][4], double rgb[3]
	 )
{
  double u,v ;
  int p[2];
  double xyz[3] = {0,0,0};
  for (u = ulo ; u <= uhi ; u += 0.005) {
    for(v = vlo; v <= vhi ; v += 0.005) {
      func(u, v, xyz) ;
      D3d_mat_mult_pt(xyz,mat,xyz);
      p[0] = (400/tan(half_angle)) * (xyz[0]/xyz[2]) + 400;
      p[1] = (400/tan(half_angle)) * (xyz[1]/xyz[2]) + 400;

      light(u,v,mat,rgb,func);
      if( (p[0] < 800 && p[0] >-1) &&
          (p[1] < 800 && p[1] >-1) &&
          (zbuff[p[0]][p[1]] > xyz[2]))
        {
          zbuff[p[0]][p[1]] = xyz[2];
          G_point(p[0],p[1]) ;
        }
    }
  }
  return 1;
}

int main()
{
  eye[0] = 0; eye[1] = 0; eye[2] = 0;
  ls[0] = 3; ls[1] = 10; ls[2] = .1;
  clear_zbuff();
  int i, Tn, Ttypelist[100] ;
  double Tvlist[100] ;
  double mat[4][4],imat[4][4] ;

  G_init_graphics(800,800) ;
  G_rgb(0,0,0) ;
  G_clear() ;

  //build view window to draw in 3d
  //build_view_window();

  //---------------------------------------------------------
  // spheres
  double rgb[3];

  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  0 ; Tn++ ; //2.5
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  0.0 ; Tn++ ; //.5
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  10 ; Tn++ ; //7
  //Ttypelist[Tn] = SZ ; Tvlist[Tn] =    3 ; Tn ++;
  D3d_make_movement_sequence_matrix (mat,imat,
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  rgb[0] = .1; rgb[1] = .6; rgb[2] = .3;
  plot_3d(0.0, 2*M_PI, 0,2*M_PI, f8, mat, rgb) ;

  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  1 ; Tn++ ; //2.5
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  0 ; Tn++ ; //.5
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  6 ; Tn++ ; //7
  //Ttypelist[Tn] = SZ ; Tvlist[Tn] =    3 ; Tn ++;
  D3d_make_movement_sequence_matrix (mat,imat,
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  rgb[0] = 1; rgb[1] = 1; rgb[2] = 0;
  //plot_3d(0.0, 2*M_PI, 0,2*M_PI, f8, mat, rgb) ;
  
  G_wait_key() ;
  

}



