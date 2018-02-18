#include <FPT.h>
#include <D3d_matrix.h>

int main()
{
  int i, Tn, Ttypelist[100] ;
  double Tvlist[100] ;
  double u,x,y,p[3] ;
  double mat[4][4],imat[4][4] ;

  G_init_graphics(800,800) ;
  G_rgb(0,0,0) ;
  G_clear() ;
  G_rgb(1,1,0) ;
  
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   50.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  100.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  300.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  500.0 ; Tn++ ;
  D3d_make_movement_sequence_matrix (mat,imat,
                                     Tn,
                                     Ttypelist,
                                     Tvlist) ;

  for (u = 0.25*M_PI ; u <= 1.50*M_PI ; u += 0.01) {
    x = cos(u) ;
    y = sin(u) ;
    p[0] = x ; p[1] = y ; p[2] = 0 ;
    D3d_mat_mult_pt(p,mat,p) ;
    G_point(p[0],p[1]) ;
  }

  G_wait_key() ;
}



