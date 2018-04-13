#include <FPT.h>
#include <D3d_matrix.h>
#include <xwd_tools.h>
//--------------------------------------------general display variables:
double half_angle = 25 *M_PI/180;
#define WW  800
#define WH  800

#define BG 0.18
#define X 0
#define Y 1
#define Z 2

#define SPHERE 0
#define HYPERB 1

#define MMN 20 //maximum number of shapes
int MN = 0;
//--------------------------------------lighting and viewpoint variables
#define AMBIENT 0.2
#define DIFFUSE_MAX 0.5
#define SPEC_POW 30

double eye[3], ls[3];
double coi[3], up[3];
double VIEW[4][4], iVIEW[4][4];
double light_in_eye_space[3];

// scale: used for translating 3d --> 2d screen
void scale(double *v, double scale, int l) {
  int i;
  for(i=0;i<l;i++) {
    v[i] *= scale;
  }
}


int Light_Model (double irgb[3],
                 double s[3],
                 double p[3],
                 double n[3],
                 double argb[3])
// s,p,n in eyespace

// irgb == inherent color of object (input to this function)
// s = location of start of ray (probably the eye)
// p = point on object (input to this function)
// n = normal to the object at p (input to this function)
// argb == actual color of object (output of this function)
// globals : AMBIENT, MAX_DIFFUSE, SPECPOW, light_in_eye_space[3]

// return 1 if successful, 0 if error
{

  double len ;
  double N[3] ; 
  len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;
  if (len == 0) return 0 ;
  N[0] = n[0]/len ;  N[1] = n[1]/len ;  N[2] = n[2]/len ;

  double E[3] ;
  E[0] = s[0] - p[0] ; 
  E[1] = s[1] - p[1] ; 
  E[2] = s[2] - p[2] ; 
  len = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]) ;
  if (len == 0) return 0 ;
  E[0] /= len ;  E[1] /= len ;  E[2] /= len ;
  double NdotE = N[0]*E[0] + N[1]*E[1] + N[2]*E[2] ;

  double L[3] ;
  L[0] = light_in_eye_space[0] - p[0] ; 
  L[1] = light_in_eye_space[1] - p[1] ; 
  L[2] = light_in_eye_space[2] - p[2] ; 
  len = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) ;
  if (len == 0) return 0 ;
  L[0] /= len ;  L[1] /= len ;  L[2] /= len ;
  double NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2] ;





  double max_ambient_and_diffuse = AMBIENT + DIFFUSE_MAX;
     // this needs to occur BEFORE you possibly jump to LLL below




  double intensity ;
  if (NdotL*NdotE < 0) {
    // eye and light are on opposite sides of polygon
    intensity = AMBIENT ; 
    goto LLL ;
  } else if ((NdotL < 0) && (NdotE < 0)) {
    // eye and light on same side but normal pointing "wrong" way
    N[0] *= (-1.0) ;    N[1] *= (-1.0) ;    N[2] *= (-1.0) ; 
    NdotL *= (-1.0) ;
    NdotE *= (-1.0) ;   // don't use NdotE below, probably should eliminate this
  }


  // ignore Blinn's variant
  double R[3] ; // Reflection vector of incoming light
  R[0] = 2*NdotL*N[0] - L[0] ;
  R[1] = 2*NdotL*N[1] - L[1] ;
  R[2] = 2*NdotL*N[2] - L[2] ;

  double EdotR = E[0]*R[0] + E[1]*R[1] + E[2]*R[2] ;

  double diffuse ;
  if (NdotL <= 0.0) { diffuse = 0.0 ; }
  else { diffuse = DIFFUSE_MAX*NdotL ; }

  double specular ;
  if (EdotR <= 0.0) { specular = 0.0 ; }
  else { specular = (1.0 - max_ambient_and_diffuse)*pow(EdotR,SPEC_POW) ;}

  // printf("%lf %lf\n",diffuse,specular) ;
  intensity = AMBIENT + diffuse + specular ;



 LLL : ;

  double f,g ;
  if (intensity <= max_ambient_and_diffuse) {
    f = intensity / max_ambient_and_diffuse ;
    argb[0] = f * irgb[0] ;
    argb[1] = f * irgb[1] ;
    argb[2] = f * irgb[2] ;
  } else {
    f = (intensity - max_ambient_and_diffuse) / 
                           (1.0 - max_ambient_and_diffuse) ;
    g = 1.0 - f ;
    argb[0] = g * irgb[0] + f ;
    argb[1] = g * irgb[1] + f ;
    argb[2] = g * irgb[2] + f ;
  }

  return 1 ;
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

void shade(double rgb[3], double sum, double intensity,
           double fin[3]) {
  double ratio;
  int i;
  if (intensity <= sum) {
    ratio = intensity/sum;
    for(i=0;i<3;i++) {
      fin[i] = rgb[i]*ratio;
    }
  } else {
    ratio = (intensity-sum)/(1-sum);
    for(i=0;i<3;i++){
      fin[i] = (1-ratio)*rgb[i] + ratio;
    }
  }
}
void compute_best_intersection(double mat[MMN][4][4], double imat[MMN][4][4],
                               int *shape,
                               double p1[3], double p2[3],
                               double xyzi[4]) {
  double a, b, c; // variables of quadratic equation
  double pos[3], neg[3], ltemp;
  double pa[3]; // point a, in simplified space
  double pb[3]; // point b, in simplified space
  double dir[3];// direction vector of equation
  double tn, tp;

  xyzi[3] = -1;
  double shortest_length = 10000;
  int i,j;

  for(i=0;i<MN; i++) {
    D3d_mat_mult_pt(pa,imat[i],p1);
    D3d_mat_mult_pt(pb,imat[i],p2);

    // find direction vector
    dir[X] = pb[X] - pa[X];
    dir[Y] = pb[Y] - pa[Y];
    dir[Z] = pb[Z] - pa[Z];
    if(dir[X] == 0 && dir[Y] == 0 && dir[Z] == 0) {continue;}

    // finding parameters of quadratic equation:
    if (shape[i] == SPHERE){ // at^2 + bt + c = 0
      a = dir[X]*dir[X] + dir[Y]*dir[Y] + dir[Z]*dir[Z];
      b = 2 * (pa[X]*dir[X] + pa[Y]*dir[Y] + pa[Z]*dir[Z]);
      c =  pa[X]*pa[X] + pa[Y]*pa[Y] + pa[Z]*pa[Z] - 1;
    } else if (shape[i] == HYPERB) {
      a = dir[X]*dir[X] - dir[Y]*dir[Y] + dir[Z]*dir[Z];
      b = 2 * (pa[X]*dir[X] - pa[Y]*dir[Y] + pa[Z]*dir[Z]);
      c = pa[X]*pa[X] - pa[Y]*pa[Y] + pa[Z]*pa[Z] - 1;
    }
    // solving quadratic equation
    tp = (-b + sqrt(b*b - 4*a*c)) / (2 * a);
    tn = (-b - sqrt(b*b - 4*a*c)) / (2 * a);

    if(isnan(tp) && isnan(tn)){continue;} //they don't intersect generally

    for(j=0;j<3;j++){
      pos[j] = pa[j] + tp*dir[j];
      neg[j] = pa[j] + tn*dir[j];
    }

    if(pos[Y] < 1.5 && pos[Y] > -1.5) { // this is only necessary for hyperboloids,
                                    // but won't mess up spheres
      //test for pos being the shortest
      ltemp = tp;
      if(ltemp < shortest_length && tp >= 0){ //if shortest and t is positive
        D3d_mat_mult_pt(pos,mat[i],pos);
        xyzi[0] = pos[X];
        xyzi[1] = pos[Y];
        xyzi[2] = pos[Z];
        xyzi[3] = i;
        shortest_length = ltemp;
      }
    }

    //test for neg being the shortest
    if(neg[Y] < 1.5 && neg[Y] > -1.5) { // still only necessary for hyperboloids
      ltemp = tn;
      if (ltemp < shortest_length && tn >= 0){ // if shortest and t is positive
        D3d_mat_mult_pt(neg,mat[i],neg);
        xyzi[0] = neg[X];
        xyzi[1] = neg[Y];
        xyzi[2] = neg[Z];
        xyzi[3] = i;
        shortest_length = ltemp;
      }
    }
  }
}

int plot (int map, int *shape,
          double mat[MMN][4][4], double imat[MMN][4][4],
          double base_rgb[5][3]
	 )
{
  int p[2];
  double screen[3] = {0,0,1};
  double xyz[3] = {0,0,0};
  double th = tan(half_angle);
  double xyzi[4] = {0,0,0,-1}, shade[3]; // returns
  double origin[3] = {0,0,0};
  double norm[3];

  int i,j, obj;
  for (i=0; i < WW; i++) {
    for(j=0; j < WH; j++) {
      screen[X] = -th + i * (2*th) / WW;
      screen[Y] = -th + j * (2*th) / WH;
      screen[Z] = 1;

      compute_best_intersection(mat, imat, shape, origin, screen, xyzi);
      if(xyzi[3] == -1) {
        set_xwd_map_color(map,i,j,BG,BG,BG); // set pixel to background color
        continue;
      }

      obj = xyzi[3];

      //compute partial derivative:
      double xyz_obj[3]; //xyz in object space
      D3d_mat_mult_pt(xyz_obj, imat[obj], xyzi);

      double par[3] = {2*xyz_obj[X], 2*xyz_obj[Y], 2*xyz_obj[Z]};
      if(shape[obj] == HYPERB) {par[1] = -par[1];}
      norm[0] = imat[obj][0][0]*par[X] + imat[obj][1][0]*par[Y] + imat[obj][2][0]*par[Z];
      norm[1] = imat[obj][0][1]*par[X] + imat[obj][1][1]*par[Y] + imat[obj][2][1]*par[Z];
      norm[2] = imat[obj][0][2]*par[X] + imat[obj][1][2]*par[Y] + imat[obj][2][2]*par[Z];
      Light_Model(base_rgb[obj],origin,xyzi,norm,shade);

      set_xwd_map_color(map, i, j, shade[0], shade[1], shade[2]);
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
  int shape[MMN];
  double reflectivity[MMN];
  //---------------------------------------------------------
  // create VIEW matrix
  eye[0] = 4;
  eye[1] = 3;
  eye[2] = -10;
  up[0] = eye[0];
  up[1] = eye[1] + 1;
  up[2] = eye[2];
  coi[0] = 0; coi[1] = 0; coi[2] = 0;
  D3d_make_identity(VIEW); D3d_make_identity(iVIEW);
  D3d_view(VIEW,iVIEW,eye,coi,up);

  ls[X] = 10;
  ls[Y] = 10;
  ls[Z] = -5;
  D3d_mat_mult_pt(light_in_eye_space,VIEW,ls);
  //---------------------------------------------------------

  // keep track of matrix to define single shape
  int Tn, Ttypelist[100] ;
  double Tvlist[100] ;

  //---------------------------------------------------------
  // hyperboloid
  Tn = 0 ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] = -1 ; Tn++;
  D3d_make_movement_sequence_matrix (mat[MN],imat[MN],
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  reflectivity[MN] = 0;
  shape[MN] = HYPERB; MN++;
  //---------------------------------------------------------
  // sphere1
  Tn = 0 ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  1.1 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] = 2    ; Tn++ ;

  D3d_make_movement_sequence_matrix (mat[MN],imat[MN],
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  reflectivity[MN] = 1;
  shape[MN] = SPHERE; MN++;

  int i;
  for(i = 0; i < MN; i++){
    D3d_mat_mult(mat[i], mat[i],VIEW);
    D3d_mat_mult(imat[i], iVIEW,imat[i]);
  }

  double rgb[5][3] = {{1,1,0.5},
                      {0.5,0.1,0.6},
                      {1,0,0},
                      {0,1,0},
                      {0,0,1}};

  plot(map,shape,mat,imat,rgb);
  xwd_map_to_named_xwd_file(map, filename);

  return 1;
}
