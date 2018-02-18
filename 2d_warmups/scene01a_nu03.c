#include <FPT.h>
#include <D3d_matrix.h>

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






int plot_with_norms (double ulo, double uhi,
          int (*func)(double u, double xy[2]),
	  double mat[4][4], int norms
	 ) 
{
  double u,p[3] ;
  double twovals[3] = {0,0,0};
  double a[3]={0,0,0}, b[3] = {0,0,0}, c[2] = {0,0};

  for (u = ulo ; u <= uhi ; u += 0.001) {
    func(u, twovals) ;
    D3d_mat_mult_pt(twovals,mat,twovals);

    //normals
    if(norms && ((int)(u*100)%10 == 0)){
      func(u-0.01, a);
      D3d_mat_mult_pt(a,mat,a);

      b[0] = (a[1]-twovals[1]);
      b[1] = (a[0]-twovals[0]);
      normalize(b,2);
      scale(b,10,2);
      
      G_line(twovals[0],twovals[1],
	     twovals[0] - b[0],
	     twovals[1] + b[1]
	     );
    }else {
      p[0] = twovals[0] ; p[1] = twovals[1] ; p[2] = 0 ;
      G_point(p[0],p[1]) ;
    }
  }
  return 1;
}


int plot (double ulo, double uhi, int (*func)(double u, double xy[3]),
	  double mat[4][4]) {
  return plot_with_norms(ulo, uhi, func, mat,0);
}


int main()
{
  int i, Tn, Ttypelist[100] ;
  double Tvlist[100] ;
  double mat[4][4],imat[4][4] ;

  G_init_graphics(800,800) ;
  G_rgb(0,0,0) ;
  G_clear() ;


  
  //---------------------------------------------------------
  // circle
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   50.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  100.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  300.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  500.0 ; Tn++ ;
  D3d_make_movement_sequence_matrix (mat,imat,
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;

  G_rgb(1,0,0) ;
  plot(0.25*M_PI, 1.50*M_PI,  f1, mat) ;
  G_wait_key() ;

  //---------------------------------------------------------
  // sum4
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   30.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   60.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  300.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  300.0 ; Tn++ ;
  D3d_make_movement_sequence_matrix (mat,imat,
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;

  G_rgb(0,1,0) ;
  plot(0.50*M_PI, 1.75*M_PI,  f2, mat) ;
  G_wait_key() ;

  //---------------------------------------------------------
  // square
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  150.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  100.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  500.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  500.0 ; Tn++ ;
  D3d_make_movement_sequence_matrix (mat,imat,
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  G_rgb(0,0,1) ;
  plot(0.0, 2*M_PI,  f3, mat) ;
  G_wait_key() ;

  //---------------------------------------------------------
  // astroid
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   80.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   40.0 ; Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] =   45.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  500.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  300.0 ; Tn++ ;
  D3d_make_movement_sequence_matrix (mat,imat,
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  G_rgb(1,1,0) ;
  plot(0.0, 2*M_PI,  f4, mat) ;
  G_wait_key() ;

  //---------------------------------------------------------
  // hyperbola : right branch
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = NY ; Tvlist[Tn] =    0.0 ; Tn++ ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  100.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  100.0 ; Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] =    0.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  250.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  250.0 ; Tn++ ;
  D3d_make_movement_sequence_matrix (mat,imat,
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  G_rgb(0,1,1) ;
  plot(-1, 2,  f5, mat) ;
  G_wait_key() ;

  //---------------------------------------------------------
  // parabola
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  150.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   50.0 ; Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] =   60.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  250.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  250.0 ; Tn++ ;
  D3d_make_movement_sequence_matrix (mat,imat,
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  G_rgb(1,0,1) ;
  plot(-1, 2,  f6, mat) ;
  G_wait_key() ;

  //---------------------------------------------------------
  // lemon
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  150.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  150.0 ; Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] =   60.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  600.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  150.0 ; Tn++ ;
  D3d_make_movement_sequence_matrix (mat,imat,
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  G_rgb(1,1,1) ;
  plot(0, 2*M_PI,  f7, mat) ;
  G_wait_key() ;

  //---------------------------------------------------------
  // sphere

  

}



