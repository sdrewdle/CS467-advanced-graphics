#include <FPT.h>
#include <D3d_matrix.h>
#include <parametric_descriptions.h>
// general display variables:
double half_angle = 25 *M_PI/180;
double hither = .1, yon = 10, view_window[6][4];   //not implemented yet
#define WW  800
#define WH  800

// lighting and viewpoint variables
#define AMBIENT 0.2
#define DIFFUSE_MAX 0.5
#define SPEC_POW 30

double zbuff[WW][WH];
double eye[3], ls[3];
double coi[3], up[3];
double VIEW[4][4];

void clear_zbuff() {
  int i,j;
  for(int i=0;i<800;i++){for(int j=0;j<800;j++) {zbuff[i][j]=pow(10,50);}}
}

void init() {
  clear_zbuff();
  eye[0] = 0; eye[1] =  2.5; eye[2] = -10;
  ls [0] = 2;  ls[1] = 2;  ls[2] = -4;
  coi[0] = 0; coi[1] =  0; coi[2] =  0;

  up[0] = eye[0];
  up[1] = eye[1] + 1;
  up[2] = eye[2];

  D3d_make_identity(VIEW);
  double _[4][4];
  D3d_view(VIEW,_,eye,coi,up);

}
//-----------------------------------------------------------universal tools


void scale(double *v, double scale, int l) {
  int i;
  for(i=0;i<l;i++) {
    v[i] *= scale;
  }
}

//-----------------------------------------------------------3d tools

void reflection(double l[3], double n[3], double r[3]) {
  int i;
  for(i=0;i<3;i++) {
    r[i] = 2*dot(n,l)*n[i]/dot(n,n)-l[i];
  }
  normalize(r,3);
}



//-----------------------------------------------------------3d functions
// sphere :  x^2 + y^2 + z^2 = 0


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
    specular = (1-AMBIENT-DIFFUSE_MAX) * pow(edr, SPEC_POW);
  }

  if(ndl <= 0) {
    diffuse = 0;
  } else {
    diffuse = DIFFUSE_MAX * ndl;
  }

  //calculate correct rgb value and call G_rgb to set
  shade(rgb,AMBIENT+DIFFUSE_MAX, AMBIENT + diffuse + specular);
}

int plot_3d_with_inc (double ulo, double uhi, double vlo, double vhi,
	     int (*func)(double u, double v, double xyz[3]),
                      double mat[4][4], double rgb[3],
                      double incU, double incV
	 )
{
  double u,v ;
  int p[2];
  double xyz[3] = {0,0,0};
  for (u = ulo ; u <= uhi ; u += incU) {
    for(v = vlo; v <= vhi ; v += incV) {
      func(u, v, xyz) ;
      D3d_mat_mult_pt(xyz,mat,xyz);
      D3d_mat_mult_pt(xyz,VIEW,xyz);
      p[0] = (400/tan(half_angle)) * (xyz[0]/xyz[2]) + WW/2;
      p[1] = (400/tan(half_angle)) * (xyz[1]/xyz[2]) + WH/2;

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
int plot_3d (double ulo, double uhi, double vlo, double vhi,
             int (*func)(double u, double v, double xyz[3]),
             double mat[4][4], double rgb[3]
             )
{
  return plot_3d_with_inc (ulo,uhi,vlo,vhi,
                           func, mat, rgb,
                           0.01, 0.01);
}

int main()
{
  int i, Tn, Ttypelist[100] ;
  double Tvlist[100] ;
  double mat[4][4],imat[4][4] ;

  init();

  G_init_graphics(WW,WH) ;
  G_rgb(0,0,0) ;
  G_clear() ;

  //---------------------------------------------------------
  // sphere1
  double rgb[3];

  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  3 ; Tn++ ;
 
  D3d_make_movement_sequence_matrix (mat,imat,
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  //D3d_view(mat,imat,eye,coi,up);
  rgb[0] = .1; rgb[1] = .6; rgb[2] = .3;
  plot_3d(0.0, 2*M_PI, 0,2*M_PI, f8, mat, rgb);

  
  //---------------------------------------------------------
  // sphere2
  D3d_make_identity(mat); D3d_make_identity(imat);
  Tn = 0 ;

  Ttypelist[Tn] = TY ; Tvlist[Tn] =  -3 ; Tn++ ; 
  D3d_make_movement_sequence_matrix (mat,imat,
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  //D3d_view(mat,imat,eye,coi,up);
  rgb[0] = .5; rgb[1] = .5; rgb[2] = .8;
  plot_3d(0.0, 2*M_PI, 0,2*M_PI, f8, mat, rgb) ;

  G_rgb(1,1,0);
  G_fill_circle(WW/2,WH/2,4);

  //--------------------------------------------------------
  // hyperbaloid
  D3d_make_identity(mat); D3d_make_identity(imat);

  Tn=0;
  Ttypelist[Tn] = RX; Tvlist[Tn] = 90; Tn++;
  Ttypelist[Tn] = SY; Tvlist[Tn] = 2;  Tn++;
  Ttypelist[Tn] = SX; Tvlist[Tn] = 0.7; Tn++;
  Ttypelist[Tn] = SZ; Tvlist[Tn] = 0.7; Tn++;
  D3d_make_movement_sequence_matrix (mat, imat,
                                     Tn,
                                     Ttypelist,
                                     Tvlist);
  rgb[0] = 0; rgb[1] = 1; rgb[2] = 1;
  double v = M_PI/2 - 0.65;
  plot_3d_with_inc(-v, v, 0, 2*M_PI, f9, mat, rgb,0.002,0.002);

  G_wait_key();
  return 1;
}



/**

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


double add(double xyz[3], int i){
return view_window[i][0] * xyz[0] +
view_window[i][1] * xyz[1] +
view_window[i][2] * xyz[2] +
view_window[i][3];
}
**/
