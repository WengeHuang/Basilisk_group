#include "navier-stokes/centered.h"
//#include "contact.h"
#include "two-phase.h"
//#include "vof.h"
#include "tension.h"
#include "view.h"


double R = 3.;
double U = 10.;
int maxlevel = 8;
double temp = 10.;

//f[left] = ((y-8)*(y-2) <= 0 ? 1.0: 0.0);
//u.t[left] = dirichlet(0.);
//u.n[left] = dirichlet(0.);
vector h[];
double theta0 = 60;
//h.t[bottom] = contact_angle (theta0*pi/180.);

//scalar f[], f2[], * interfaces = {f, f2};   // command out two-phase.h interface

p[top] = dirichlet(0.);
p[right] = dirichlet(0.);
p[left] = dirichlet(0.);
f[bottom] = 0.;
f[top] = 0.;
//u.t[right] = dirichlet(0.);

int main(){
  mu2 = 1./30.; //gas phase
  mu1 = 100./30.; //liquid phase
  rho2 = 1.; //gas phase
  rho1 = 1000;//liquid phase
  f.sigma = 1.;
  init_grid(1 << 8);
  L0 = 10;
  X0 = -L0/2;
  Y0 = 0.;
  //f.height = h;
  run();
}

event init(i = 0){
  //refine(sq(x+5) + sq(y) < sq(R + 0.25) && sq(x+5) + sq(y) > sq(R - 0.25) && level < maxlevel);
  
  //fraction(f, (sq(R) - sq(x) - sq(y-1.5))*(sq(x) + sq(y+0.6) - sq(1.2)));
  fraction(f2, (sq(x) + sq(y) - sq(1.2))*(- sq(x) - sq(y) + sq(1.6)));
  fraction(f, (sq(x) + sq(y) - sq(1.2))*(- sq(x) - sq(y) + sq(2.)));
  //foreach()
    //u.y[] = f2[]*U*10;
}

/*
event add_velocity(i++){
   foreach(){
    if (x >= 0.){
    u.x[] = f2[]*U*cos(atan(y/(x)))/sqrt(sq(x) + sq(y))*fm.x[];
    u.y[] = f2[]*U*sin(atan(y/(x)))/sqrt(sq(x) + sq(y))*fm.y[];
    }
    else{
    u.x[] = -f2[]*U*cos(atan(y/(x)))/sqrt(sq(x) + sq(y))*fm.x[];
    u.y[] = -f2[]*U*sin(atan(y/(x)))/sqrt(sq(x) + sq(y))*fm.y[];
    }
    }
    
}*/

event acceleration (i++) {
  face vector av = a;
  foreach_face(x){
    if (x > 0.){
    av.x[] += f2[]*U*cos(atan(y/x));
   }
    else if (x < 0.){
    av.x[] += -f2[]*U*cos(atan(y/x));
    }
   }
   foreach_face(y){
    if (x > 0.){
    av.y[] += f2[]*U*sin(atan(y/x));
   }
    else if (x < 0.){
    av.y[] += -f2[]*U*sin(atan(y/x));
    }
   }
}

event adapt(i++)
  adapt_wavelet((scalar *){u,f}, (double []){0.02, 0.02, 0.001}, maxlevel-1);
  event viewer(t=0.05;t+=0.1;t<=temp){
  clear();
  //view(psi = -pi/2, tx = 0.5);
  view(ty = -0.5);
  box();
  //cells();
  draw_vof("f", lw = 3);
  draw_vof("f", filled = 1, fc = {0.9, 0.4, 0.2});
  draw_vof("f2", lw = 3);
  squares("u.y",linear = true);
  save("radial_10.mp4");
}














