#include <FPT.h>
int numpoints[10];   
double swidth, sheight ;
double x[10][2000],y[10][2000],z[10][2000];
int i,j,k,npoly[10];
int psize[10][2000];
int con[10][2000][20];
double half=M_PI/4;
int obnum;
double light[3]={0,0,-3};

void findMidpoint(double p[], int obj){
  int i;
  p[0]=0; p[1]=0; p[2]=0;
  for(i=0;i<numpoints[obj];i++){
    p[0]+=x[obj][i]/(1.0*numpoints[obj]);
    p[1]+=y[obj][i]/(1.0*numpoints[obj]);
    p[2]+=z[obj][i]/(1.0*numpoints[obj]);
  }
}

void translate(int obj,double dx, double dy, double dz){
  double m[4][4],minv[4][4];
  D3d_make_identity(m) ;
  D3d_make_identity(minv) ;
  D3d_translate(m,minv,dx,dy,dz);
  D3d_mat_mult_points(x[obj],y[obj],z[obj],m,x[obj],y[obj],z[obj],numpoints[obj]);
}


void scale(int obj, double sx, double sy,double sz){
  double m[4][4],minv[4][4];
  D3d_make_identity(m) ;
  D3d_make_identity(minv) ;
  D3d_scale(m,minv,sx,sy,sz);
  D3d_mat_mult_points(x[obj],y[obj],z[obj],m,x[obj],y[obj],z[obj],numpoints[obj]);
}

void rotate_x(int obj, double dt){
  double m[4][4],minv[4][4],p[3];
  D3d_make_identity(m) ;
  D3d_make_identity(minv) ;
  findMidpoint(p,obj);
  D3d_translate(m,minv,-p[0],-p[1],-p[2]);
  D3d_rotate_x(m,minv,dt);
  D3d_translate(m,minv,p[0],p[1],p[2]);
  D3d_mat_mult_points(x[obj],y[obj],z[obj],m,x[obj],y[obj],z[obj],numpoints[obj]);


}

void rotate_y(int obj, double dt){
  double m[4][4],minv[4][4],p[3];
  D3d_make_identity(m) ;
  D3d_make_identity(minv) ;
  findMidpoint(p,obj);
  D3d_translate(m,minv,-p[0],-p[1],-p[2]);
  D3d_rotate_y(m,minv,dt);
  D3d_translate(m,minv,p[0],p[1],p[2]);
  D3d_mat_mult_points(x[obj],y[obj],z[obj],m,x[obj],y[obj],z[obj],numpoints[obj]);
}


void rotate_z(int obj, double dt){
  double m[4][4],minv[4][4],p[3];
  D3d_make_identity(m) ;
  D3d_make_identity(minv) ;
  findMidpoint(p,obj);
  D3d_translate(m,minv,-p[0],-p[1],-p[2]);
  D3d_rotate_z(m,minv,dt);
  D3d_translate(m,minv,p[0],p[1],p[2]);
  D3d_mat_mult_points(x[obj],y[obj],z[obj],m,x[obj],y[obj],z[obj],numpoints[obj]);
}

double dot(double x[],double y[]){
  double t=0;  
  int i;
  for(i=0;i<3;i++){
    t+=x[i]*y[i];
  }
  return(t);
}

int compare (const void *p, const void *q) {
  double *a =(double*) p;
  double *b = (double*) q;
  if ((a[0]) < (b[0])) return 1;
  if ((a[0]) > (b[0])) return -1;
  return 0;
}

void getNormal(int obj,int polynum,double n[]){
  double temp[3],p[3],temp2[3];
  findMidpoint(p,obj);
  temp[0]=x[obj][con[obj][polynum][1]]-x[obj][con[obj][polynum][0]];
  temp[1]=y[obj][con[obj][polynum][1]]-y[obj][con[obj][polynum][0]];
  temp[2]=z[obj][con[obj][polynum][1]]-z[obj][con[obj][polynum][0]];

  temp2[0]=x[obj][con[obj][polynum][2]]-x[obj][con[obj][polynum][1]];
  temp2[1]=y[obj][con[obj][polynum][2]]-y[obj][con[obj][polynum][1]];
  temp2[2]=z[obj][con[obj][polynum][2]]-z[obj][con[obj][polynum][1]];

  D3d_x_product(temp,temp,temp2);

  temp[0]=-temp[0]; temp[1]=-temp[1]; temp[2]=-temp[2];

  n[0]=temp[0]; n[1]=temp[1]; n[2]=temp[2];
}

void getReflect(double l[],double n[],double r[]){
  int i;
  for(i=0;i<3;i++){
    r[i]=2*dot(l,n)*n[i]/dot(n,n)-l[i];
  }
  
}



void normalize(double v[]){
  double x=sqrt(dot(v,v));
  int i;
  for(i=0;i<3;i++){
    v[i]=v[i]/(1.0*x);
  }
}



void shade(int obj, double sum,double intensity){
  double ratio;
  double r,g,b;
  if(obj==0){
    r=.7; g=.2; b=.7;
  }
  if(obj==1){
    r=.6; g=.6; b=.2;
  }
  if(intensity==sum){
    G_rgb(r,g,b);
  }
  else if(intensity>sum){
    ratio=intensity/sum;
    G_rgb(r*ratio,g*ratio,b*ratio);
  }
  else{
    ratio=(intensity-sum)/(1-sum);
    G_rgb(r+(1-r)*ratio,g+(1-g)*ratio,b+(1-b)*ratio);
  }

}

void itsLIT(int obj, int i){
  int j;
  double ambient=.4, specular, diff;
  double diffmax=.59;
  double spec = 1-ambient-diffmax;
  double l[3]={light[0]-x[obj][con[obj][i][0]],
	       light[1]-y[obj][con[obj][i][0]],
	       light[2]-z[obj][con[obj][i][0]]};
  double n[3];
  getNormal(obj,i,n);
  normalize(n); 
  double e[3]={-x[obj][con[obj][i][0]],
	       -y[obj][con[obj][i][0]],
	       -z[obj][con[obj][i][0]]};
  double r[3];
  normalize(l);
  if(1) {//flip normal
    for(j=0;j<3;j++) {n[j] = -n[j];}
  }
  getReflect(l,n,r);
  normalize(e);

  double cosA = dot(e,r);
  if(cosA < 0) {
    specular = 0;
  } else {
    specular = spec*pow(cosA,20);
  }
  double cosB = dot(n,l);
  if(cosB<0) {
    diff = 0;
  } else {
    diff = diffmax*dot(n,l);
  }
  double intensity=ambient+diff+specular;
  shade(obj,ambient+diffmax,intensity);
}


void drawPoly(int obj,int i){
  double xact[200],yact[200];
  for(j=0;j<psize[obj][i];j++){
    xact[j]=x[obj][con[obj][i][j]]*500/(1.0*z[obj][con[obj][i][j]]*tan(half))+500;
    yact[j]=y[obj][con[obj][i][j]]*500/(1.0*z[obj][con[obj][i][j]]*tan(half))+500;
  }
  itsLIT(obj,i);
  G_fill_polygon(xact,yact,psize[obj][i]);
  //G_rgb(0,0,0);
  //G_polygon(xact,yact,psize[obj][i]);
}




void drawObject(){
  double xactual,yactual,xnext,ynext;
  int i,j,j1,k,c;
  double pt[10][8000],d[8000];
  c=0;
  for(i=0;i<obnum;i++){
    for(j=0;j<npoly[i];j++){
      pt[i][j]+=x[i][con[i][j][0]]*x[i][con[i][j][0]];
      pt[i][j]+=y[i][con[i][j][0]]*y[i][con[i][j][0]];
      pt[i][j]+=z[i][con[i][j][0]]*z[i][con[i][j][0]];
      pt[i][j]=sqrt(pt[i][j]);
      d[c]=pt[i][j];
      c++;
    }
  }
  qsort(d,c,sizeof(double),compare);
  
  for(i=0;i<c;i++){
    for(j=0;j<obnum;j++){
      for(k=0;k<npoly[j];k++){
	if(pt[j][k]==d[i]){
	  drawPoly(j,k);
	  pt[j][k]=0;
	  break;
	}
      }
    }
  }
}


void parseZ(int obj){
  k=G_wait_key();
  while((k==65362)||(k==65364)||(k==65361)||(k==65363)){
    if(k==65362){
      translate(obj,0,0,1);
      G_rgb(0,0,0);
      G_clear();
      drawObject();
    }
    if(k==65364){
      translate(obj,0,0,-1);
      G_rgb(0,0,0);
      G_clear();
      drawObject();
    }
    if(k==65361){
      rotate_z(obj,.04);
      G_rgb(0,0,0);
      G_clear();
      drawObject();
    }
    if(k==65363){
      rotate_z(obj,-.04);
      G_rgb(0,0,0);
      G_clear();
      drawObject();
    }

    k=G_wait_key();
  }
}



void parseX(int obj){
  k=G_wait_key();
  while((k==65362)||(k==65364)||(k==65361)||(k==65363)){
    if(k==65362){
      translate(obj,1,0,0);
      G_rgb(0,0,0);
      G_clear();
      drawObject();
    }
    if(k==65364){
      translate(obj,-1,0,0);
      G_rgb(0,0,0);
      G_clear();
      drawObject();
    }
    if(k==65361){
      rotate_x(obj,.04);
      G_rgb(0,0,0);
      G_clear();
      drawObject();
    }
    if(k==65363){
      rotate_x(obj,-.04);
      G_rgb(0,0,0);
      G_clear();
      drawObject();
    }

    k=G_wait_key();
  }
}





void parseY(int obj){
  k=G_wait_key();
  while((k==65362)||(k==65364)||(k==65361)||(k==65363)){
    if(k==65362){
      translate(obj,0,1,0);
      G_rgb(0,0,0);
      G_clear();
      drawObject();
    }
    if(k==65364){
      translate(obj,0,-1,0);
      G_rgb(0,0,0);
      G_clear();
      drawObject();
    }
    if(k==65361){
      rotate_y(obj,.04);
      G_rgb(0,0,0);
      G_clear();
      drawObject();
    }
    if(k==65363){
      rotate_y(obj,-.04);
      G_rgb(0,0,0);
      G_clear();
      drawObject();
    }

    k=G_wait_key();
  }
}






int main(int argc, char **argv)
{
   FILE *q;
   for(k=0;k<argc-1;k++){
     q=fopen(argv[k+1],"r"); //r for read 
     if(q==NULL){
       printf("don't be a dummy idiot\n");
       exit(0);
     }
     fscanf(q,"%d",&numpoints[k]);
     for(i=0;i<numpoints[k];i++){
       fscanf(q,"%lf %lf %lf",&x[k][i],&y[k][i],&z[k][i]);
     }
     fscanf(q,"%d",&npoly[k]);
     for(i=0;i<npoly[k];i++){
       fscanf(q,"%d",&psize[k][i]);
       for(j=0;j<psize[k][i];j++){
	 fscanf(q,"%d",&con[k][i][j]);
       }
     }
   }
   obnum=argc-1;

   // must do this before you do 'almost' any other
   // graphical tasks 
   swidth = 1000 ;  sheight = 1000 ;
   G_init_graphics (swidth, sheight) ;
   G_rgb(0,0,0);
   G_clear();
   k=0;
   int obj=0;


   drawObject();


   k=G_wait_key();
   while((k-49>=0)&&(k-49<obnum)){
     obj=k-49;
     k=G_wait_key();
     while((k=='z')||(k=='x')||(k=='y')){
       if(k=='z'){
	 parseZ(obj);
       }
       if(k=='x'){
	 parseX(obj);
       }
       if(k=='y'){
	 parseY(obj);
       }
     }
   }
   G_close() ; // terminate graphics...probably not fatal if forgotten

}
