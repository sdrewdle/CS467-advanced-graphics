#include <FPT.h>
#include <D3d_matrix.h>

double min(double a, double b) {
  if (a<b) return a;
  return b;
}

double max(double a, double b) {
  if (a>b) return a;
  return b;
}

void draw_circ(double t0, double t1, double m[4][4]) {
  double pt[3];

  double i;
  for (i=min(t0,t1); i<=max(t0,t1); i+=0.01) {
    pt[0] = cos(i); pt[1] = sin(i);
    D3d_mat_mult_pt (pt,m,pt);
    G_fill_circle(pt[0], pt[1],0.75);
  }
}

void draw_sum4(double t0, double t1, double m[4][4]) {
  double pt[3];

  double i;
  for (i=min(t0,t1); i<=max(t0,t1); i+=0.01) {
    pt[0] = fabs(cos(i)); pt[1] = fabs(sin(i));
    D3d_mat_mult_pt (pt,m,pt);
    G_fill_circle(pt[0], pt[1],0.75);
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

  int tl1[4] = {SX, SY, TX, TY};
  double pl1[4] = {30,60,300,300};
  D3d_make_movement_sequence_matrix(m,_,4,tl1,pl1);
  draw_sum4(0*M_PI, 1.75*M_PI, m);

  G_wait_key();
}
