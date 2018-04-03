#include <FPT.h>
#include <D3d_matrix.h>
#include <parametric_descriptions.h>
#include <xwd_tools.h>
//--------------------------------------------general display variables:
double half_angle = 25 *M_PI/180;
#define WW  800
#define WH  800

#define BG 0.18
#define X 0
#define Y 1
#define Z 2

#define MMN 20 //maximum number of shapes
int MN = 0;
//--------------------------------------lighting and viewpoint variables
#define AMBIENT 0.2
#define DIFFUSE_MAX 0.5
#define SPEC_POW 30

double zbuff[WW][WH];
double eye[3], ls[3];
double coi[3], up[3];
double VIEW[4][4], iVIEW[4][4];

// scale: used for translating 3d --> 2d screen
void scale(double *v, double scale, int l) {
  int i;
  for(i=0;i<l;i++) {
    v[i] *= scale;
  }
}

// reflection over a normal vector
// guaranteed |r| = 1
void reflection(double l[3], double n[3], double r[3]) {
  int i;
  for(i=0;i<3;i++) {
    r[i] = 2*dot(n,l)*n[i]/dot(n,n)-l[i];
  }
  normalize(r,3);
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

// it's light, but with a checkerboard pattern
void checker(double rgb1[3], double rgb2[3],
             double sum, double intensity, double u, double v,
             double rgb[3])
{
  v += M_PI/2;
  double ratio, r, g, b;
  if((int)(u*100/30) % 2 == 0) {
    if((int)(v*100/30) % 2 == 0){
      r = rgb1[0]; g = rgb1[1]; b = rgb1[2];
    } else {
      r = rgb2[0]; g = rgb2[1]; b = rgb2[2];
    }
  } else {
    if((int)(v*100/20) % 2 == 1){
      r = rgb1[0]; g = rgb1[1]; b = rgb1[2];
    } else {
      r = rgb2[0]; g = rgb2[1]; b = rgb2[2];
    }
  }

  if (intensity <= sum) {
    ratio = intensity/sum;
    rgb[0] = r * ratio;
    rgb[1] = g * ratio;
    rgb[2] = b * ratio;

  } else {
    ratio = (intensity-sum)/(1-sum);
    rgb[0] = (1-ratio)*r + ratio;
    rgb[1] = (1-ratio)*g + ratio;
    rgb[2] = (1-ratio)*b + ratio;
  }
}

/* sets the rgb value to the appropriate color
 * u,v are the current angles
 * rgb is an array of the base color
 * func is the same function used to get the point;
 *** used to approx the norm
 */
void light(double u, double v, double mat[4][4],
           double RGB[3],
           int (*func)(double uu, double vv, double xyz[3]),
           double rgb[3])
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
  double white[3] = {1,1,1};
  checker(RGB,white,AMBIENT+DIFFUSE_MAX, AMBIENT + diffuse + specular, u, v, rgb);
}

void compute_best_intersection(double mat[MMN][4][4], double imat[MMN][4][4],
                               int (*f[MMN])(double u, double v, double xyz[3]),
                               double p1[3], double p2[3],
                               double uvi[3]) {
  double a, b, c; // variables of quadratic equation
  double pos[3], neg[3], ltemp;
  double pa[3]; // point a, in simplified space
  double pb[3]; // point b, in simplified space
  double dir[3];// direction vector of equation
  double tn, tp;

  uvi[2] = -1;
  double shortest_length = 10000;
  int i;

  for(i=0;i<MN; i++) {
    D3d_mat_mult_pt(pa,imat[i],p1);
    D3d_mat_mult_pt(pb,imat[i],p2);

    // find direction vector
    dir[X] = pb[X] - pa[X];
    dir[Y] = pb[Y] - pa[Y];
    dir[Z] = pb[Z] - pa[Z];
    if(dir[X] == 0 && dir[Y] == 0 && dir[Z] == 0) {continue;}

    // finding parameters of quadratic equation:
    // at^2 + bt + c = 0
    a = dir[X]*dir[X] + dir[Y]*dir[Y] + dir[Z]*dir[Z];
    b = 2 * (pa[X]*dir[X] + pa[Y]*dir[Y] + pa[Z]*dir[Z]);
    c =  pa[X]*pa[X] + pa[Y]*pa[Y] + pa[Z]*pa[Z] - 1;

    // solving quadratic equation
    tp = (-b + sqrt(b*b - 4*a*c)) / (2 * a);
    tn = (-b - sqrt(b*b - 4*a*c)) / (2 * a);

    //check for nan; means they don't intersect
    if(isnan(tp) && isnan(tn)){continue;}

    pos[X] = pa[X] + tp*dir[X];
    pos[Y] = pa[Y] + tp*dir[Y];
    pos[Z] = pa[Z] + tp*dir[Z];

    neg[X] = pa[X] + tn*dir[X];
    neg[Y] = pa[Y] + tn*dir[Y];
    neg[Z] = pa[Z] + tn*dir[Z];


    // back in normal view 
    D3d_mat_mult_pt(pos,mat[i],pos);
    D3d_mat_mult_pt(neg,mat[i],neg);

    //test for pos being the shortest
    ltemp = sqrt(pow(pos[X]-p2[X],2) + pow(pos[Y]-p2[Y],2) + pow(pos[Z]-p2[Z],2) );
    if(ltemp < shortest_length && tp >= 0){
      D3d_mat_mult_pt(pos,imat[i],pos);
      uvi[0] = atan2(pos[Y],pos[X]);
      uvi[1] = atan2(pos[Z],pos[X]);
      uvi[2] = i;
      shortest_length = ltemp;

      D3d_mat_mult_pt(pos,mat[i],pos);
    }

    //test for neg being the shortest
    ltemp = sqrt(pow(neg[X]-p2[X],2) + pow(neg[Y]-p2[Y],2) + pow(neg[Z]-p2[Z],2) );
    if (ltemp < shortest_length && tn >= 0){
      D3d_mat_mult_pt(neg,imat[i],neg);
      uvi[0] = atan2(neg[Y],neg[X]);
      uvi[1] = atan2(neg[Z],neg[X]);
      uvi[2] = i;
      shortest_length = ltemp;

      D3d_mat_mult_pt(neg,mat[i],neg);
    }

  }
}

int plot_with_inc (int map, int (*func[MMN])(double u, double v, double xyz[3]),
                   double mat[MMN][4][4], double imat[MMN][4][4]
	 )
{
  int p[2];
  double screen[3] = {0,0,1};
  double xyz[3] = {0,0,0};
  double th = tan(half_angle);
  double uvi[3] = {0,0,0}; // grr give me python tupling
  double origin[3] = {0,0,0};

  int i,j;
  for (i=0; i < WW; i++) {
    for(j=0; j < WH; j++) {
      screen[X] = -th + i * (2*th) / WW;
      screen[Y] = -th + j * (2*th) / WH;
      screen[Z] = 1;

      compute_best_intersection(mat, imat, func, origin, screen, uvi);
      if(uvi[2] == -1) {
        set_xwd_map_color(map,i,j,BG,BG,BG); // set pixel to background color
        continue;
      }
      (func[ (int)(uvi[2]) ])(uvi[0], uvi[1], xyz);

      //light(uvi[0],uvi[1],mat[uvi[i]],RGB,func, rgb);
      set_xwd_map_color(map, i, j, 1, 1, 0.5);
    }
  }
  return 1;
}

int main()
{
  int map;
  char filename[25];
  sprintf(filename,"pic0000.xwd");

  map = create_new_xwd_map(WW,WH);
  xwd_map_to_named_xwd_file(map,filename);


  double mat[MMN][4][4],imat[MMN][4][4] ;
  int (*f[MMN])(double u, double v, double xyz[3]);
  //---------------------------------------------------------
  // create VIEW matrix
  eye[0] = 0;
  eye[1] = 0;
  eye[2] = -10;
  up[0] = eye[0];
  up[1] = eye[1] + 1;
  up[2] = eye[2];
  coi[0] = 0; coi[1] = 0; coi[2] = 0;
  D3d_make_identity(VIEW); D3d_make_identity(iVIEW);
  D3d_view(VIEW,iVIEW,eye,coi,up);

  //---------------------------------------------------------

  // keep track of matrix to define single shape
  int Tn, Ttypelist[100] ;
  double Tvlist[100] ;

  //---------------------------------------------------------
  // sphere1 (inc by 0.005)
  Tn = 0 ;
  //Ttypelist[Tn] = TX ; Tvlist[Tn] =  4 ; Tn++ ;

  D3d_make_movement_sequence_matrix (mat[MN],imat[MN],
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  f[MN] = sphere; MN++;


  int i;
  for(i = 0; i < MN; i++){
    D3d_mat_mult(mat[i], mat[i],VIEW);
    D3d_mat_mult(imat[i], imat[i],iVIEW);
  }

  plot_with_inc(map,f,mat,imat);
  xwd_map_to_named_xwd_file(map, filename);

  return 1;
}
