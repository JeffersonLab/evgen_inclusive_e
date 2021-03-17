#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "F1F221Wrapper.hh"

static const double kDEG = 3.14159265358979323846 / 180.0;
static const double degrad = 180.0/M_PI;

int main()
{
    const double M = 0.9383;
    const double proton_mass=0.9383;
    double F1, F2,F1n,F2n, F1p, F2p, F1d, F2d, r;
    double F1be, F2be, rbe;
    double F1dqe, F2dqe;
    double  W2, Nu,q2,w2,nu;
    double x,theta/*,Ep*/,W,xj;
    //const double E=10.618;
    const double E=4.74;
    F1F221Wrapper *model = new F1F221Wrapper();
   //const int n=17;
   // double x[n]  = {0.185,0.215,0.245,0.275,0.305,0.335,0.365,0.395,0.425,0.455,0.49,0.53,0.57,0.61,0.655,0.705,0.76};
  // double Q2[n]  = {2.709,3.059,3.447,3.88,4.313,4.708,5.122,5.549,6.024,6.397,6.924,7.473,8.044,8.616,9.242,9.998,10.689};
   const int n=9;
  // double  x[n]  = {0.214,0.299,0.456,0.494,0.533,0.579,0.629,0.686,0.745};
   double  Ep[n]  = {0.599,0.798,1.118,1.188,1.257,1.336,1.416,1.504,1.593};
   double  Q2[n]  = {1.659,2.209,3.094,3.285,3.472,3.694,3.909,4.149,4.387};
 //  double Q2; 
  double xs1=0.0;
	printf("{ ");
    for(int i=0;i<n;i++){
		Nu=E-Ep[i];
		//Nu=Q2[i]/(2*proton_mass*x[i]);
		//Ep=E-Nu;
		//Nu=Q2[i]/(2*proton_mass*x[i]);
                //theta=acos(1-Q2[i]/(2*E*Ep));
		x=Q2[i]/(2*proton_mass*Nu);
                //theta=acos(1-Q2[i]/(2*E*Ep[i]));
                theta=45*kDEG;
                W2 = M * M + 2. * M * Nu - Q2[i];
            model->GetF1F2IN21(2, 3, Q2[i], W2, F1d, F2d);
            model->GetF1F2QE21(2, 3, Q2[i], W2, F1dqe, F2dqe);
            F1= F1d;
            F2= F2d;
       // printf(" %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf \n", x, Q2[i], Ep[i],theta*degrad,W2,Nu,E,F1d, F2d,F1dqe,F2dqe);
     xs1=(2./137.*Ep[i]/Q2[i]*cos(theta/2))*(2./137.*Ep[i]/Q2[i]*cos(theta/2));
     //        printf("%5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.8lf\n", x[i], Q2[i], Ep,theta*degrad,W2,Nu,E,F1, F2,1000*xs1);
    xs1 = xs1 * (2. / proton_mass * F1 * tan(fabs(theta) / 2.) * tan(fabs(theta) / 2.) + F2 / Nu);
     //      printf("%5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.13lf\n", x[i], Q2[i], Ep,theta*degrad,W2,Nu,E,F1, F2,1000*xs1);
    xs1 = xs1 * 389.379;

   double   dxs=1000*xs1;
            //printf("%5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf \n", x[i], Q2[i], Ep,theta*degrad,W2,Nu,F1, F2,1000*xs1);
  printf("%4.3f,", dxs);
}
	printf("},");
  /*  for (int i = 1; i <= 2; i++) {
        q2 = 1.0 * i;
        for (int j = 1; j <= 5; j++) {
            w2 = 1. + 0.5 * j;
            nu = (w2 - M * M + q2) / 2. / M;
            xj= q2 / 2. / M / nu;
            model->GetF1F2IN21(0., 1., q2, w2, F1n, F2n);
            model->GetF1F2IN21(1., 1., q2, w2, F1p, F2p);
            model->GetF1F2IN21(1., 2., q2, w2, F1d, F2d);
            model->GetF1F2IN21(4., 9., q2, w2, F1be, F2be);
            printf("%5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf \n", q2, w2, xj, F2n, F2p, F2d, F1p, F1be);
        }
    }*/

 /*  for (int i = 1; i <= 3; i++) {
        Q2 = 1.0 * i;
    for(int j=1;j<7;j++){
        W2 =0.5*j;
        model->GetF1F2QE21(1., 2., Q2, W2, F1dqe, F2dqe);
        printf("%5.3lf %5.3lf %5.3lf %5.3lf\n",W2, Q2, F1dqe, F2dqe);
    }
  }*/
}
