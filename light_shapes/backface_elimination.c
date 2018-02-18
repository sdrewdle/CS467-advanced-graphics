#include <FPT.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <D3d_matrix.h>

int window_height = 1800, window_width = 1800;
double half_angle = 30*M_PI/180;

double mean(double *a, int l) {
  int i;
  double sum = 0;
  for(i=0; i<l; i++) {
    sum += a[i];
  }
  return sum/l;
}

double largest_radius(double *x, double *y, int p) {
  double center_x = x[p], center_y = y[p],
    radius = 0, temp;
  int i;
  for(i=0;i<p;i++) {
    temp = sqrt( pow(y[i] - center_y, 2) + pow(x[i] - center_x, 2) );
    if(temp > radius) {radius = temp;}
  }
  return radius;
}

void translate(char direction, double *x, double *y, double *z, int p) {
  double dx = 0,
         dy = 0,
         dz = 0;
  if       (direction=='Q') {
    dx = -1;
  } else if(direction=='S') {
    dx = 1;
  } else if(direction=='R') {
    dz=1;
  } else if(direction=='T') {
    dz=-1;
  } else if(direction=='-') {
    dy=-1;
  } else if(direction=='=') {
    dy=1;
  }

  double m[4][4], _[4][4];
  D3d_make_identity(m);
  D3d_translate(m,_,dx,dy,dz);
  D3d_mat_mult_points(x,y,z,m,x,y,z,p+1);
}

void rotate(char axis, double *x, double *y ,double *z, int p) {
  double m[4][4], _[4][4];
  double cx = x[p], cy = y[p], cz = z[p];
  D3d_make_identity(m);
  D3d_translate(m,_,0-cx,0-cy,0-cz);
  
  if(axis == 'x'){
    D3d_rotate_x(m,_, .01);
  }else if (axis == 'y') {
    D3d_rotate_y(m,_, .01);
  }else if (axis == 'z') {
    D3d_rotate_z(m,_, .01);
  }

  D3d_translate(m,_,cx,cy,cz);
  D3d_mat_mult_points(x,y,z,m,x,y,z, p+1);
  
}

double dot(double *u, double *v) {
  return (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]);
}
double vector_length(double *v) {
  return sqrt(dot(v,v));
}
void vector_to(double x1, double y1, double z1,
	       double x2, double y2, double z2,
	       double v[3]){
  //create a vector from p2 to p1
  v[2] = z2-z1;
  v[1] = y2-y1;
  v[0] = x2-x1;

  double l = vector_length(v);
  v[2]/=l; v[1]/=l; v[0]/=l;
}

int acute(double *u, double *v) {
  return dot(u,v)>0;
}

void orthogonal(double u[3], double v[3], double c[3], double orthog[3]){
  //computes the outward facing orthogonal
  double t[3];
  D3d_x_product(t,u,v);

  if(acute(t,c)||1){
    orthog[0] = 0 - t[0];
    orthog[1] = 0 - t[1];
    orthog[2] = 0 - t[2];
  } else {
    orthog[0] = t[0];
    orthog[1] = t[1];
    orthog[2] = t[2];
  }
}

void draw (double *x, double *y,double *z,int p, int polys[][90],int q) {
  G_rgb(.1,.1,.1);
  G_clear();
  double cx = x[p], cy = y[p], cz = z[p];
  int i,j, ni, pj, nj;
  double x2d[p], y2d[p], v[3], u[3], c[3];
  for(i=0;i<p;i++){
    x2d[i] = x[i]*(window_width/2) / (z[i] * tan(half_angle)) + window_width/2;
    y2d[i] = y[i]*(window_width/2) / (z[i] * tan(half_angle)) + window_height/2;
  }
  
  double temp_x[20], temp_y[20];
  for(i=0;i<q;i++) { //each polygon
    j=polys[i][2];
    pj = polys[i][1];
    nj = polys[i][3];

    vector_to(x[pj],y[pj],z[pj], x[j],y[j],z[j], u);
    vector_to(x[nj],y[nj],z[nj], x[j],y[j],z[j], v);
    vector_to(cx,cy,cz,          x[j],y[j],z[j], c);
    
    orthogonal(u,v,c,v);
    vector_to(0,0,0, x[j],y[j],z[j], u);
    
    if(acute(u,v)){
      // make shape:
      for(j=1;j<=polys[i][0];j++) {
	temp_x[j-1] = x2d[polys[i][j]];
	temp_y[j-1] = y2d[polys[i][j]];
      }
      
      G_rgb(.9,.7,.7);
      G_polygon(temp_x,temp_y,polys[i][0]);
    }
  }
}

int main(int argc, char **argv) {
  G_init_graphics(window_width,window_height);
  int i,j,k;

  //an array of all the files
  FILE* fp_array[10];
  int p[10],
      q[10];
  int largest_p = 0;

  for(i=1;i<argc;i++) {
    fp_array[i-1] = fopen(argv[i], "r");
    if(fp_array[i-1] == NULL) {
      printf("could not open file, %s\n", argv[i]);
      exit(0);
    }

    //find num of points from each file
    fscanf(fp_array[i-1], "%d", &p[i-1]);
    if(p[i-1] > largest_p) {largest_p = p[i-1];}
  }
  double x[10][largest_p+1],
         y[10][largest_p+1],
         z[10][largest_p+1];
  int    polys[10][1000][90];
  double smallest_z = 1000;
  //loop through every file, filling arrays
  for(i=0;i<argc-1;i++) {
    //fill the x and y parallel arrays
    for(j=0; j<p[i]; j++) {
      fscanf(fp_array[i], "%lf", &x[i][j]);
      fscanf(fp_array[i], "%lf", &y[i][j]);
      fscanf(fp_array[i], "%lf", &z[i][j]);
      x[i][j]+= .01; y[i][j] += .01; z[i][j] += .01;
      if(smallest_z > z[i][j]){smallest_z = z[i][j];}
    }
    //center points
    x[i][p[i]] = mean(x[i], p[i]);
    y[i][p[i]] = mean(y[i], p[i]);
    z[i][p[i]] = mean(z[i], p[i]);
    
    //read the polygon information
    fscanf(fp_array[i], "%d", &q[i]);
    for(j=0;j<q[i];j++) {
      fscanf(fp_array[i],"%d", &polys[i][j][0]);
      for(k=1;k<=polys[i][j][0];k++) {
	fscanf(fp_array[i],"%d", &polys[i][j][k]);
      }
    }
    draw(x[i],y[i],z[i],p[i],polys[i],q[i]);
  }

  //select object
  char c;
  while(1){
    printf("choose an object\n");
    c = G_wait_key();
    
    if(c == 'q') {break;}
    if(c >= '1' && c <= '9') {
      i = c - '0' - 1;
     
      printf("choose an action\n");
      c = G_wait_key();
      if(c == 't') {         //translate
	printf("translate chosen\n");
	while(1){
	  c = G_wait_key();
	  if(c=='q') {break;}
	  translate(c,x[i],y[i],z[i],p[i]);
	  draw(x[i],y[i],z[i],p[i],polys[i],q[i]);
	}
      } else if (c == 'r') { //rotate
	printf("rotate chosen, choose direction x,y, or z\n");
	while(1) {
	  c = G_wait_key();
	  if(c=='q') {break;}
	  rotate(c,x[i],y[i],z[i],p[i]);
	  draw(x[i],y[i],z[i],p[i],polys[i],q[i]);
	}
      } else if (c == 'z') { //zoom - don't do this.....
	//do something with the half angle
      }
    } else {
      printf("%c is not a valid object\n", c);
    }

    
  }
}
