//#include "axi.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "reduced.h"





double R = 3.;
double U = 1.;
int maxlevel = 8;
double temp = 10.;

//f[left] = ((y-8)*(y-2) <= 0 ? 1.0: 0.0);
//u.t[left] = dirichlet(0.);
//u.n[left] = dirichlet(0.);


f[right] = 0.;
u.t[right] = dirichlet(0.);



vector h[];
double theta0 = 150;




int main(){
  mu2 = 1./30.; //gas phase
  mu1 = 100./30.; //liquid phase
  rho2 = 1.; //gas phase
  rho1 = 100;//liquid phase
  f.sigma = 100.;
  init_grid(1 << 8);
  L0 = 10;
  X0 = -L0/2;
  Y0 = 0.;

  //double Bond = 0.308;
  //G.x = -Bond;
 
  //f.height = h;
  run();
}

event init(i = 0){
  refine(sq(x+5) + sq(y) < sq(R + 0.25) && sq(x+5) + sq(y) > sq(R - 0.25) && level < maxlevel);
  fraction(f, sq(R) - sq(x+5) - sq(y-5));
  //foreach()
    //u.x[] = -f[]*U*0.5;
}


event CA_waves(i++){

    if (i % 50 >= 1 && i % 50 <= 25) {
        // Perform action 1 for every range of 5 steps
           h.t[left] = contact_angle (theta0*pi/180.);
           f.height = h;
    } else if (i % 50 >= 26 && i % 50 <= 50) {
        // Perform action 2 for every range of 5 steps
           h.t[left] = contact_angle ((theta0-60)*pi/180.);
           f.height = h;
    }
}


/*
event contact_angl1(i = 0; i +=10)
{   

   h.t[left] = contact_angle (theta0*pi/180.);
   f.height = h;
}

event contact_angl2(i = 21; i +=10)
{   

   h.t[left] = contact_angle ((theta0-60)*pi/180.);
   f.height = h;
}*/

event adapt(i++)
  adapt_wavelet((scalar *){u,f}, (double []){0.02, 0.02, 0.001}, maxlevel-1);
 
  event viewer(t=0.05;t+=0.1;t<=temp*0.1){
  clear();
  //view(psi = -pi/2, theta = 0.1, phi = 0.2);
  view(psi = -pi/2, tx = 0.5);
  box();
  //cells();
  draw_vof("f", lw = 3);
  draw_vof("f", filled = 1, fc = {0.9, 0.4, 0.2});
  //mirror({0,1}) {
  //  draw_vof("f");
  //  draw_vof("f", filled = 1, fc = {0.9, 0.4, 0.2});
    //squares("u.x",linear = true);
  //}
  save("jump_capillary.mp4");
}
