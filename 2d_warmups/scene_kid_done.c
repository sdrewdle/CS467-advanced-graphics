#include <FPT.h>
#include <D3d_matrix.h>

double sgn(double x) {
  return (x>0) ? 1: ((x<0) ? -1 : 0);
}

double min(double a, double b) {
  if (a<b) return a;
  return b;
}

double max(double a, double b) {
  if (a>b) return a;
  return b;
}

void plot_polar(double r, double t, double m[4][4]) {
  double pt[3] = {r * cos(t), r * sin(t), 0};
  D3d_mat_mult_pt(pt,m,pt);
  G_fill_circle(pt[0],pt[1],0.75);
}

void plot_par(double x, double y, double m[4][4]) {
  double pt[3] = {x,y,0};
  D3d_mat_mult_pt (pt,m,pt);
  G_fill_circle(pt[0], pt[1],0.75);
}

void draw_circ(double t0, double t1, double m[4][4]) {
  double i;
  for (i=min(t0,t1); i<=max(t0,t1); i+=0.01) {
    plot_par(cos(i), sin(i), m);
  }
}

void draw_sum4(double t0, double t1, double m[4][4]) {
  double i;
  for (i=min(t0,t1); i<=max(t0,t1); i+=0.01) {
    plot_polar(1/sqrt(sqrt( pow(cos(i),4) + pow(sin(i),4)  )),i,m);
  }
}

void draw_square(double t0, double t1, double m[4][4]) {
  double i,c,s,mm; double a[3],b[3];
  for(i=min(t0,t1); i<=max(t0,t1); i+=0.01) {
    c = cos(i); s = sin(i);
    plot_par(sgn(c)*pow(c,2),sgn(s)*pow(s,2), m);
    if((int)(i*10)%30 == 0) {
      printf("%lf\n", i);
    }
  }
}

void drawstroid(double t0, double t1, double m[4][4]) {
  double i,c,s;
  for(i=min(t0,t1); i<max(t0,t1); i+=0.01) {
    c = cos(i); s = sin(i);
    plot_par(sgn(c)*pow(c,4),sgn(s)*pow(s,4), m);
  }
}

void tbh_i_think_youre_being_a_bit_hyperbolic(
					      double t0, double t1,
					      double m[4][4]) {
  double i,c,s;
  for(i=min(t0,t1); i<max(t0,t1); i+=0.001) {
    plot_par(cosh(i), sinh(i), m);
  }
}

void draw_parabola(double t0, double t1, double m[4][4]) {
  double i;
  for(i=min(t0,t1); i<=max(t0,t1);i+=0.001) {
    plot_par(i,i*i,m);
  }
}

void LEMON(double t0, double t1, double m[4][4]) {
  double i;
  for(i=min(t0,t1);i<=max(t0,t1);i+=0.01) {
    plot_par(pow(cos(i),3),sin(i),m);
  }
}

int main() {
  G_init_graphics(800,800);

  double m[4][4], _[4][4];
  int tl[4] = {SX, SY, TX, TY};
  double pl[4] = {50,100,300,500};
  D3d_make_movement_sequence_matrix(m,_,4,tl,pl);
  
  G_rgb(1,1,1); G_clear(); G_rgb(0,0,0);
  draw_circ(.25*M_PI, 1.5*M_PI, m);

  tl[0] = SX; tl[2] = TX;
  tl[1] = SY, tl[3] = TY;
  pl[0] = 30; pl[2] = 300;
  pl[1] = 60, pl[3] = 300;
  
  D3d_make_movement_sequence_matrix(m,_,4,tl,pl);
  draw_sum4(.5*M_PI, 1.75*M_PI, m);

  pl[0] = 150; pl[2] = 500;
  pl[1] = 100, pl[3] = 500;

  D3d_make_movement_sequence_matrix(m,_,4,tl,pl);
  draw_square(0, 2*M_PI,m);

  int tl5[5] = {SX,SY,RZ,TX,TY};
  double pl5[5] = {80,40,50,500,300};

  D3d_make_movement_sequence_matrix(m,_,5,tl5,pl5);
  drawstroid(0, 2*M_PI,m);

  int tl6[6] = {NY,SX,SY,RZ,TX,TY};
  double pl6[6] = {0,100,100,0,250,250};

  D3d_make_movement_sequence_matrix(m,_,6,tl6,pl6);
  tbh_i_think_youre_being_a_bit_hyperbolic(-1,2,m);

  tl5[0] = SX; pl5[0] = 150;
  tl5[1] = SY; pl5[1] = 50;
  tl5[2] = RZ; pl5[2] = 60;
  tl5[3] = TX; pl5[3] = 250;
  tl5[4] = TY; pl5[4] = 250;

  D3d_make_movement_sequence_matrix(m,_,5,tl5,pl5);
  draw_parabola(-1,2,m);

  tl5[0] = SX; pl5[0] = 150;
  tl5[1] = SY; pl5[1] = 150;
  tl5[2] = RZ; pl5[2] = 60;
  tl5[3] = TX; pl5[3] = 600;
  tl5[4] = TY; pl5[4] = 150;

  D3d_make_movement_sequence_matrix(m,_,5,tl5,pl5);
  LEMON(0,2*M_PI,m);
  
  G_wait_key();
}
