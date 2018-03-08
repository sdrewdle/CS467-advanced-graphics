#include <FPT.h>
#include <D3d_matrix.h>

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

int plot (double ulo, double uhi,
          int (*func)(double u, double xy[2]),
          double mat[4][4]
	 ) 
{
  double u,p[3] ;
  double twovals[2] ;

  for (u = ulo ; u <= uhi ; u += 0.01) {
    func(u, twovals) ;
    p[0] = twovals[0] ; p[1] = twovals[1] ; p[2] = 0 ;
    D3d_mat_mult_pt(p,mat,p) ;
    G_point(p[0],p[1]) ;
  }
}

// n = (b-a) x (c-a)
int plot_normals(double ulo, double uhi,
          int (*func)(double u, double xy[2]),
          double mat[4][4]
	 ) 
{
  double u,p[3] ;
  double twovals[2] ;
  double twovalsb[2];

  double a,b;
  double xx,yy;

  
  for (u = ulo ; u <= uhi ; u += 0.1) {
    func(u, twovals) ;
    func((u-0.1), twovalsb);
    D3d_mat_mult_pt(twovals,mat,twovals);
    D3d_mat_mult_pt(twovalsb,mat,twovalsb);

    
    a = twovalsb[0] - twovals[0];
    b = twovalsb[1] - twovals[1];

    xx = twovals[0] - b;
    yy = twovals[1] + a;
    
    G_line(twovals[0],twovals[1],xx,yy) ;
  }
}





int main()
{
  int i, Tn, Ttypelist[100] ;
  double Tvlist[100] ;
  double mat[4][4],imat[4][4] ;
  double mat1[4][4], imat1[4][4];

  G_init_graphics(800,800) ;
  G_rgb(0.1, 0.1, 0.1) ;
  G_clear() ;

  //---------------------------------------------------------
  // circle 1
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   50.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  100.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  300.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  500.0 ; Tn++ ;
  D3d_make_movement_sequence_matrix (mat,imat,
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  //---------------------------------------------------------
  // circle 2
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  200.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   20.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  30.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  400.0 ; Tn++ ;
  D3d_make_movement_sequence_matrix (mat1,imat1,
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  G_rgb(.1,.7,.9) ;
  plot(0, 2*M_PI,  f1, mat) ;

  double p1[3], p2[3];
  G_rgb(1,1,1);
  
  G_wait_click(p1);
  G_fill_circle(p1[0],p1[1],3);

  double x_old = p1[0];
  double y_old = p1[1];

  G_wait_click(p2);
  G_fill_circle(p2[0],p2[1],3);

  D3d_mat_mult_pt(p1,imat,p1);  //translate back into simpler space
  D3d_mat_mult_pt(p2,imat,p2);

  double M = (p2[1]-p1[1]) / (p2[0]-p1[0]);
  double B = p1[1]  - M*p1[0];
  double a = M*M + 1;
  double b = 2 * M * B;
  double c = B*B - 1;

  double pos = (-b + sqrt(b*b - 4*a*c)) / (2 * a);
  double neg = (-b - sqrt(b*b - 4*a*c)) / (2 * a);
  double x,y, xx,yy;

  double xpos = pos; double ypos = M*xpos + B;
  double xneg = neg; double yneg = M*xneg + B;

 
  if(pow(xpos-p1[0],2) + pow(ypos-p1[1],2) < pow(xneg-p1[0],2) + pow(yneg-p1[1],2)){
    p1[0] = xpos;
    p1[1] = ypos;
  } else {
    p1[0] = xneg;
    p1[1] = yneg;
  }
  double u = atan2(p1[1],p1[0]);

  D3d_mat_mult_pt(p1,mat,p1);
  G_rgb(1,1,1);
  G_fill_circle(p1[0],p1[1],3);

  f1(u+0.001,p2);
  D3d_mat_mult_pt(p2,mat,p2);
  G_rgb(1,0,0);
  G_fill_circle(p2[0],p2[1],3);

  a = p2[0] - p1[0];
  b = p2[1] - p1[1];

  xx = p1[0] + 100*b;
  yy = p1[1] - 100*a;

  G_rgb(1,1,1);
  G_line(p1[0],p1[1],xx,yy) ;
  G_rgb(0,1,1);
  G_line(x_old,y_old, p1[0],p1[1]);

  //-------------------------
  p1[0] = old_x; p1[1] = old_y;


  M = (p2[1]-p1[1]) / (p2[0]-p1[0]);
  B = p1[1]  - M*p1[0];
  a = M*M + 1;
  b = 2 * M * B;
 c = B*B - 1;

  pos = (-b + sqrt(b*b - 4*a*c)) / (2 * a);
  neg = (-b - sqrt(b*b - 4*a*c)) / (2 * a);
  x,y, xx,yy;

  xpos = pos;  ypos = M*xpos + B;
  xneg = neg;  yneg = M*xneg + B;

 
  if(pow(xpos-p1[0],2) + pow(ypos-p1[1],2) < pow(xneg-p1[0],2) + pow(yneg-p1[1],2)){
    p1[0] = xpos;
    p1[1] = ypos;
  } else {
    p1[0] = xneg;
    p1[1] = yneg;
  }
   u = atan2(p1[1],p1[0]);

  D3d_mat_mult_pt(p1,mat1,p1);
  G_rgb(1,1,1);
  G_fill_circle(p1[0],p1[1],3);

  f1(u+0.001,p2);
  D3d_mat_mult_pt(p2,mat1,p2);
  G_rgb(1,0,0);
  G_fill_circle(p2[0],p2[1],3);

  a = p2[0] - p1[0];
  b = p2[1] - p1[1];

  xx = p1[0] + 100*b;
  yy = p1[1] - 100*a;

  G_rgb(1,1,1);
  G_line(p1[0],p1[1],xx,yy) ;
  G_rgb(0,1,1);
  G_line(x_old,y_old, p1[0],p1[1]);

  G_wait_key();
}



