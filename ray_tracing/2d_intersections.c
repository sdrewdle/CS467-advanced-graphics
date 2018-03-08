#include <FPT.h>
#include <D3d_matrix.h>

#define BG 0.1
#define X 0
#define Y 1

int MN;

double sgn(double v)
{
  if (v > 0) return  1.0  ;
  if (v < 0) return -1.0 ;
  return 0.0 ;
}

void gather_points(double p1[2], double p2[2]) {
  G_wait_click(p1);
  if(p1[0]<=5 && p1[1]<=5){exit(1);}
  G_rgb(.7,.8,.9);
  G_fill_circle(p1[0],p1[1],2);

  G_wait_click(p2);
  if(p2[0]<=5 && p2[1]<=5){exit(1);}
  G_fill_circle(p2[0],p2[1],2);
}

// circle :   x^2 + y^2 = 1
int f1 (double u, double xy[2])
{
  xy[0] = cos(u) ;  
  xy[1] = sin(u) ;
  return 1;
}

int plot (double ulo, double uhi,
          int (*func)(double u, double xy[2]),
          double mat[4][4]) 
{
  double u,p[3] ;
  double twovals[2] ;

  for (u = ulo ; u <= uhi ; u += 0.01) {
    func(u, twovals) ;
    p[0] = twovals[0] ; p[1] = twovals[1] ; p[2] = 0 ;
    D3d_mat_mult_pt(p,mat,p) ;
    G_point(p[0],p[1]) ;
  }
  return 1;
}

void draw_circles(double m[9][4][4], int mn) {
  double c[3] = {0,0,0};
  char num[4];
  int i;
  G_rgb(0,1,1);
  for(i=0; i<mn; i++) {
    c[0] = 0; c[1] = 0; c[2] = 0;
    plot(0,2*M_PI, f1, m[i]);
    D3d_mat_mult_pt(c,m[i],c);

    sprintf(num, "%d", i);
    G_draw_string(num,c[0]-6,c[1]-6);
  }
}

void compute_best_intersection(double mat[9][4][4], double imat[9][4][4],
                               double p1[3], double p2[3],
                               double ui[2]) {
  double slope, inter, a, b, c, pos[3], neg[3], ltemp;
  double pa[3]; // = {p1[0], p1[1], 0};
  double pb[3]; // = {p2[0], p2[1], 0};

  ui[1] = -1;
  double shortest_length = 10000;
  int i;

  for(i=0;i<MN; i++) {
    //p1[0] = P1[0]; p1[1] = P1[1]; p1[2] = P1[2];
    //p2[0] = P2[0]; p2[1] = P2[1]; p2[2] = P2[2];


    D3d_mat_mult_pt(pa,imat[i],p1);
    D3d_mat_mult_pt(pb,imat[i],p2);

    // finding equation of the line y=mx+b
    slope = (pb[1]-pa[1]) / (pb[0]-pa[0]);
    inter = pa[1]  - slope*pa[0];

    // finding parameters of quadratic equation
    a = slope*slope + 1;
    b = 2 * slope * inter;
    c = inter*inter - 1;

    // solving quadratic equation
    pos[0] = (-b + sqrt(b*b - 4*a*c)) / (2 * a);
    neg[0] = (-b - sqrt(b*b - 4*a*c)) / (2 * a);
    //check for nan
    if(isnan(pos[0]) && isnan(neg[0])){continue;}

    pos[1] = slope * pos[0] + inter;
    neg[1] = slope * neg[0] + inter;


    // back in normal view 
    D3d_mat_mult_pt(pos,mat[i],pos);
    D3d_mat_mult_pt(neg,mat[i],neg);

    ltemp = sqrt(pow(pos[0]-p2[0],2) + pow(pos[1]-p2[1],2));
    if(ltemp < shortest_length){
      D3d_mat_mult_pt(pos,imat[i],pos);
      ui[0] = atan2(pos[1],pos[0]);
      ui[1] = i;
      shortest_length = ltemp;

      D3d_mat_mult_pt(pos,mat[i],pos);
    }

    ltemp = sqrt(pow(neg[0]-p2[0],2) + pow(neg[1]-p2[1],2));
    if (ltemp < shortest_length){
      D3d_mat_mult_pt(neg,imat[i],neg);
      ui[0] = atan2(neg[1],neg[0]);
      ui[1] = i;
      shortest_length = ltemp;

      D3d_mat_mult_pt(neg,mat[i],neg);
    }

  }
}

int main()
{
  int Tn, Ttypelist[100] ;
  double Tvlist[100] ;
  int mn = 0;
  double mat[9][4][4],imat[9][4][4] ;

  G_init_graphics(800,800) ;
  G_rgb(0.1, 0.1, 0.1) ;
  G_clear() ;

  //---------------------------------------------------------
  // circle 1
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   50.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   70.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  300.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  500.0 ; Tn++ ;
  D3d_make_movement_sequence_matrix (mat[mn],imat[mn],
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  mn++;
  //---------------------------------------------------------
  // circle 2
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   70.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   20.0 ; Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] =   80.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  100.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  650.0 ; Tn++ ;
  D3d_make_movement_sequence_matrix (mat[mn],imat[mn],
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  mn++;
  //---------------------------------------------------------
  // circle 3
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   50.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  100.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  600.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  650.0 ; Tn++ ;
  D3d_make_movement_sequence_matrix (mat[mn],imat[mn],
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  mn++;
  //---------------------------------------------------------
  // circle 4
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  150.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   40.0 ; Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] =  -20.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  600.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  150.0 ; Tn++ ;
  D3d_make_movement_sequence_matrix (mat[mn],imat[mn],
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  mn++;
  //---------------------------------------------------------
  // circle 5
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   20.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   40.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  100.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  100.0 ; Tn++ ;
  D3d_make_movement_sequence_matrix (mat[mn],imat[mn],
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;
  mn++;
  //---------------------------------------------------------
  MN = mn;
  draw_circles(mat,mn);
  //---------------------------------------------------------

  G_rgb(1,.4,.4);
  G_line(0,5,5,5); G_line(5,5,5,0);
  double p1[3], p2[3],
    ui_pair[2] = {0,0},
    ci[3], //current intersection
    adj[3], // point adjacent to uth
    u;

  int i;
  while(1) {
    gather_points(p1,p2);

    compute_best_intersection(mat, imat, p1,p2, ui_pair);

    u = ui_pair[0];
    i = ui_pair[1];
    if(i<0){
      G_rgb(BG,BG,BG);
      G_fill_circle(p1[X], p1[Y], 3);
      G_fill_circle(p2[X], p2[Y], 3);
      continue;
    }

    ci[X] = cos(u);
    ci[Y] = sin(u);
    D3d_mat_mult_pt(ci,mat[i],ci);

    G_rgb(1,1,1);
    G_line(p1[X],p1[Y], ci[X], ci[Y]);

    f1(u+0.001,adj);
    D3d_mat_mult_pt(adj,mat[i],adj);

    G_rgb(1,.4,1);
    G_line(ci[X],ci[Y],
           ci[X] + 2000*(adj[Y] - ci[Y]),
           ci[Y] - 2000*(adj[X] - ci[X])
           );
  }
  return 1;
}



