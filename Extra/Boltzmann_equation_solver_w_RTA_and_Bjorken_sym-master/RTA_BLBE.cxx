// RTA_BLBE - A solver for the Relaxation Time Approximation (RTA) in the Bjorken Limit (BL) (i.e. using the Bjorken symmetry) of the Boltzmann Equation (BE)
// Author: Gojko Vujanovic 
// Licence: GNU General Public Licence version 3, see details at http://www.gnu.org/licenses/gpl-3.0.html
// Copywritht (c) 2018 
// This code was used to obtain the solution in Physcal Review C 97, no.6, 064909 (2018)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string>
#include <string.h>
#include <iostream>
using namespace std;

int const four_pi_eta_bar=1.0; // Note that eta_bar is specific shear viscosity (i.e. shear viscosity over entropy density or eta/s).
                              // The variable four_pi_eta_bar is really just 4*pi*eta/s. 

const double EPSILON=2.220446049e-16;

// Numerical solution to the RTA BLBE uses Gauss-Legendre quadrature.  

const int Np_pts=32;      // Number of points in Gauss-Legendre quadrature.
double roots_p[Np_pts];   // Roots of the Gauss-Legendre quadrature.
double weights_p[Np_pts]; // Weights of the Gauss-Legendre quadrate.

// The temperaure (or energy density) proper time (tau) profile of the Boltzmann equation in the RTA_BL limit  

const int N_tau=1025;                                     // Number of tau steps calculated.
double const tau_min=1.0;                                // Minimum proper time tau.      
double const tau_max=5.0;                               // Maximum proper time tau.
double const dtau=(tau_max-tau_min)/((double) (N_tau-1)); // Proper time step.  
 
double tau[N_tau]; // Contains the proper time grid used in calculating the temperaure (or energy density) profile.

// First piece of the 3-momentum integral of the Boltzmann distributuin in Eq. (10) of PRC 97 064909 (2018), giving the energy density.  
double I_1(double k, double a, double beta, double Delta)
{
 double result;
 if(a>=1.0){result=24.0*M_PI*k/(beta*beta*beta*beta);}
 else{
       result=4.0*M_PI*k*a*a/(1.0-a*a)*(2.0*beta*Delta*(1.0-2.0*a*a)+
                                        exp(-beta*Delta)*(
                                                          6.0*(1.0-a*a)+2.0*beta*Delta*(2.0-a*a)+beta*beta*Delta*Delta
                                                         )
                                       )*1.0/(beta*beta*beta*beta);
 }
 return(result);
}

// Second piece of the 3-momentum integral of the Boltzmann distributuin in Eq. (10) of PRC 97 064909 (2018), giving the energy density. 
// Like the first piece, the second piece can also be expressed in closed using simple functions, when there are no Electric fields (E-fields).
// The function below is written in anticipation of electric fields, where the variable Delta is in Eq. (42) of PRD 99, 016009 (2019), 
// and is zero if there are no E-fields. In the case with E-fields, the closed form expression exists using Exponential Integral functions.
// Using pre-existing librarires to evaluate the Exponential Integral evaluation, the user is bound to the precision and speed of that library. 
// So, the Exponentail Integral is done directly using a numerical integral allowing the full control on the performance of the code.
// k=degeneracy_factor/pow(2.0*M_PI,3.0). Herein degeneracy_factor=2 was chosen.
// beta is the inverse temperature, while a is a ratio of proper times defined in line 247 and 275 taking values between 0 and 1.
double I_2(double k, double a, double beta, double Delta)
{
 double result;

 if(a>=1.0){result=24.0*M_PI*k/(beta*beta*beta*beta);return(result);}
 else{
       double a2=a*a;
       double rev_a2=1.0-a2;
       double norm=2.0*M_PI*k*a/(rev_a2*sqrt(rev_a2));

       double Delta2=Delta*Delta;
       
       double result=0.0;
       int i;
       for(i=0;i<Np_pts;i++){
          double p=(1.0+roots_p[i])/(1.0-roots_p[i]);
          double x_plus =(p*rev_a2+a2*Delta)/sqrt(p*p*rev_a2+a2*Delta2);
          double x_minus=(p*rev_a2-a2*Delta)/sqrt(p*p*rev_a2+a2*Delta2);
          double integrand=2.0/((1.0-roots_p[i])*(1.0-roots_p[i]))*exp(-beta*p)*(p*p*p*rev_a2+p*a2*Delta2)*(asin(x_plus)+asin(x_minus));
                 result+=integrand*weights_p[i];
       }
       result*=norm;
       return(result);
 }
}

// Sum of the two pieces of the 3-momentum integral of the Boltzmann distributuin in Eq. (10) of PRC 97 064909 (2018), giving the energy density.  
double I(double k, double a, double beta, double Delta)
{
 double result=I_1(k,a,beta,Delta)+I_2(k,a,beta,Delta);
 return(result);
}

// Linear interpolation function.
double interp_linear(double tau, double *y)
{
 int itau;
 double tau_frac;

 itau=((int) floor((tau-tau_min)/dtau));
 tau_frac = (tau-tau_min)/dtau - ((double) itau);

 double result=(1.0-tau_frac)*y[itau]+tau_frac*y[itau+1];
 if(isnan(tau_frac)==1 || isnan(y[itau])==1 || isnan(y[itau+1])==1){
    printf("tau=%.10e itau=%d tau_frac=%10e y[itau]=%.10e y[itau+1]=%.10e\n", tau, itau, tau_frac, y[itau], y[itau+1]);
    exit(0);
 }

 return(result);
}


// Relaxation time function, which in general can be proper time (tau) and tempertaure (or beta, i.e. inverse-temperature) dependent. 
double tau_R(double tau, double *beta)
{
  return (5.0*interp_linear(tau, beta)*((double) four_pi_eta_bar)/(4.0*M_PI));
}

// Dapming function in Eq. (10) of PRC 97 064909 (2018).
double D(double tau_2, double tau_1, double *beta)
{
 int itau;
 double result=0.0;
 if((tau_2-tau_1)<=EPSILON){result=1.0;}
 else{
      for(itau=0;itau<Np_pts;itau++){
         double tau=(tau_2-tau_1)/2.0*roots_p[itau]+(tau_2+tau_1)/2.0;
         result+=1.0/tau_R(tau,beta)*weights_p[itau];
      }
      result*=(tau_2-tau_1)/2.0;
      result=exp(-result);
 }
 return(result);
}

// Loading of Gauss-Legendre roots and weights and calculating the proper time grid over which the energy density (or temperature) will be evaluated.
int load_data(void)
{

  roots_p[0]=-0.99726386184948156354;    weights_p[0]=0.0070186100094700966128;    
  roots_p[1]=-0.98561151154526833540;    weights_p[1]=0.016274394730905670603;    
  roots_p[2]=-0.96476225558750643078;    weights_p[2]=0.025392065309262059452;    
  roots_p[3]=-0.93490607593773968917;    weights_p[3]=0.034273862913021433104;    
  roots_p[4]=-0.89632115576605212396;    weights_p[4]=0.042835898022226680660;    
  roots_p[5]=-0.84936761373256997013;    weights_p[5]=0.050998059262376176194;    
  roots_p[6]=-0.79448379596794240696;    weights_p[6]=0.058684093478535547148;    
  roots_p[7]=-0.73218211874028968039;    weights_p[7]=0.065822222776361846844;    
  roots_p[8]=-0.66304426693021520097;    weights_p[8]=0.072345794108848506232;    
  roots_p[9]=-0.58771575724076232904;    weights_p[9]=0.078193895787070306478;    
  roots_p[10]=-0.50689990893222939002;   weights_p[10]=0.083311924226946755224;   
  roots_p[11]=-0.42135127613063534537;   weights_p[11]=0.087652093004403811140;   
  roots_p[12]=-0.33186860228212764979;   weights_p[12]=0.091173878695763884708;   
  roots_p[13]=-0.239287362252137074540;  weights_p[13]=0.093844399080804565636;   
  roots_p[14]=-0.144471961582796493480;  weights_p[14]=0.095638720079274859422;   
  roots_p[15]=-0.048307665687738316235;  weights_p[15]=0.096540088514727800566;   
  roots_p[16]=0.048307665687738316235;   weights_p[16]=0.096540088514727800566;   
  roots_p[17]=0.14447196158279649348;    weights_p[17]=0.095638720079274859422;   
  roots_p[18]=0.23928736225213707454;    weights_p[18]=0.093844399080804565636;   
  roots_p[19]=0.33186860228212764979;    weights_p[19]=0.091173878695763884708;   
  roots_p[20]=0.42135127613063534537;    weights_p[20]=0.087652093004403811140;   
  roots_p[21]=0.50689990893222939002;    weights_p[21]=0.083311924226946755224;   
  roots_p[22]=0.58771575724076232904;    weights_p[22]=0.078193895787070306478;   
  roots_p[23]=0.66304426693021520097;    weights_p[23]=0.072345794108848506232;   
  roots_p[24]=0.73218211874028968039;    weights_p[24]=0.065822222776361846844;   
  roots_p[25]=0.79448379596794240696;    weights_p[25]=0.058684093478535547148;   
  roots_p[26]=0.84936761373256997013;    weights_p[26]=0.050998059262376176194;   
  roots_p[27]=0.89632115576605212396;    weights_p[27]=0.042835898022226680660;   
  roots_p[28]=0.93490607593773968917;    weights_p[28]=0.034273862913021433104;   
  roots_p[29]=0.96476225558750643078;    weights_p[29]=0.025392065309262059452;   
  roots_p[30]=0.98561151154526833540;    weights_p[30]=0.016274394730905670603;   
  roots_p[31]=0.99726386184948156354;    weights_p[31]=0.0070186100094700966128;  

  // Generates the proper time (tau) grid. In general, it doesn't have to be even-spaced. 
  int itau;
  for(itau=0;itau<N_tau;itau++){
      tau[itau]=tau_min+dtau*((double) itau);
  }

 return 0;
}








    
int main(void)
{
  FILE *fp_beta; 

  // Loads the proper time Gauss-Legendre quadrature points and calculates the proper time (tau) grid. 
  int num_error;
  if( (num_error = load_data()) != 0) {
    fprintf(stderr, "Error loading data (%d)!\n", num_error);
    return 1;
  }




  const double k=2.0/pow(2.0*M_PI,3.0);

  int itau;

  // The solution to the  Eq. (10) of PRC 97 064909 (2018) will be done iteratively. 
  // Below is the starting guess for the temperature profile.
  fp_beta=fopen("e_vs_tau.dat", "w");
  double const beta_0=1.0/0.6; // Inverse temperature at tau_min. 
  double beta[N_tau];             
  for(itau=0;itau<N_tau;itau++){

      // First guess for the temeprature profile to be iterated over.
      beta[itau]=beta_0*pow(tau[itau]/tau_min,4.0/3.0); 
      fprintf(fp_beta, "%d %.10e %.10e\n", 0, tau[itau], 1.0/beta[itau]); // Saving the first guess to disk.

  }
  fclose(fp_beta);
  fp_beta = NULL;


  double energy_density[N_tau];
  energy_density[0]=48.0*M_PI*k/pow(beta[0],4.0);

  int n;
  const int N_iter=201;
  char file_name[20];
  for(n=0;n<N_iter;n++){ // Loops over 100 interations until the solution to Eq. (10) of PRC 97 064909 (2018).

	 if (n % 25 == 0){
		 sprintf(file_name,"e_vs_tau_%i.dat",n);
		 fp_beta=fopen(file_name,"w");
	 }
     // This loop calculates the updated proper time profile.
     for(itau=1;itau<N_tau;itau++){

        /***************************************************************************/
        //                                                                         //
        // The second term in  Eq. (10) of PRC 97 064909 (2018) is computed first. // 
        //                                                                         //
        /***************************************************************************/

        // Note that s is the integration variable in the second term in Eq. (10) of PRC 97 064909 (2018), which was denoted by tau' there. 
        int is;// is gives the index labeling the locations along which the proper time integral is performed.

        double s_min=tau_min;   // Lower bound of the integral in the second term of Eq. (10) of PRC 97 064909 (2018).
        double s_max=tau[itau]; // Upper bound of the integral in the second term of Eq. (10) of PRC 97 064909 (2018).

        double ans=0.0;         // Initializing the variable storing the result of the proper time integral in Eq. (10) of PRC 97 064909 (2018).
        for(is=0;is<Np_pts;is++){

           double s=(s_max-s_min)/2.0*roots_p[is]+(s_max+s_min)/2.0; // Mapping the Gauss-Legendre roots to the integration region.
           double a_s=s/tau[itau];                                   
           double beta_s=interp_linear(s,beta);                      // Linearly interpolates the inverse temeprature to the location s.
           double Delta_s=0.0;                                       // This variable was created in anticipation of including the Electric field, and is zero otherwise.
           double D_eval_s=D(tau[itau],s,beta)/tau_R(s,beta);        // Evaluates the damping function in the second term of Eq. (10) of PRC 97 064909 (2018).
           if(isnan(D_eval_s)==1){                                   // Check that the damping function was correctly evalated. Exits if not.
              printf("D_eval_s is a nan. s=%.10e a_s=%.10e beta_s=%.10e Delta_s=%.10e\n", s, a_s, beta_s, Delta_s);
              exit(0);
           }
           double I_eval_s=I(k,a_s,beta_s,Delta_s);                  // Evaluates the 3-momentum integral of the Boltzmann distribution present in the second term of Eq. (10) of PRC 97 064909 (2018).
           if(isnan(I_eval_s)==1){                                   // Checks that the 3-momentum integral of the Boltzmann distribtion was correctly evalated. Code exits if not.
              printf("I_eval_s is a nan. s=%.10e a_s=%.10e beta_s=%.10e Delta_s=%.10e D_eval_s=%.10e I_1_eval_s=%.10e I_2_eval_s=%.10e\n",
                                         s, a_s, beta_s, Delta_s, D_eval_s, I_1(k,a_s,beta_s,Delta_s), I_2(k,a_s,beta_s,Delta_s)
                    );
              exit(0);
           }
           ans+=D_eval_s*I_eval_s*weights_p[is];

        }
        ans*=(s_max-s_min)/2.0;



        /***********************************************************************/ 
        //                                                                     //
        // The first term in Eq. (10) of PRC 97 064909 (2018) is now evaluated.// 
        //                                                                     //
        /***********************************************************************/  

        double a_tau=tau_min/tau[itau];
        double beta_tau=interp_linear(tau_min,beta);    // Determines inverse temperature at tau_min. 
        double Delta_tau=0.0;
        double D_eval=D(tau[itau],tau_min,beta);        // Evaluates the damping function in the first term of Eq. (10) of PRC 97 064909 (2018).
        if(isnan(D_eval)==1){                           // Check that the damping function was correctly evalated. Exits if not.
          printf("D_eval_tau is a nan. tau=%.10e a_tau=%.10e beta_tau=%.10e Delta_tau=%.10e\n", tau[itau], a_tau, beta_tau, Delta_tau);
          exit(0);
        }
        double I_eval=I(k, a_tau, beta_tau, Delta_tau); // Evaluates the 3-momentum integral of the Boltzmann distribution present in the first term of Eq. (10) of PRC 97 064909 (2018).
        if(isnan(I_eval)==1){                           // Checks that the 3-momentum integral of the Boltzmann distribtion was correctly evalated. Code exits if not.
           printf("I_eval_tau is a nan. tau=%.10e a_tau=%.10e beta_tau=%.10e Delta_tau=%.10e D_eval=%.10e\n", tau[itau], a_tau, beta_tau, Delta_tau, D_eval);
           exit(0);
        }
        energy_density[itau]=D_eval*I_eval+ans;         // Stores the energy density profile after an iteration of Eq. (10) of PRC 97 064909 (2018).
     }
  
     // Updating the temperature profile after an iteration of Eq. (10) of PRC 97 064909 (2018), and printing the updated temperature profile to file. 
     if (n % 25 == 0) {
		 for(itau=0;itau<N_tau;itau++){
         beta[itau]=pow(48.0*M_PI*k/energy_density[itau],0.25);     
		 double e = 6.0/(pow(beta[itau],4.0)*M_PI*M_PI*pow(0.19731,3.0));
         fprintf(fp_beta, "%d %.10e %.10e %.10e\n", n+1, tau[itau], e, 1./beta[itau]);
			}	
		fclose(fp_beta); 
		fp_beta = NULL;
	 }
    
     // Done an iteration of the temperature profile.
     printf("done iteration=%d\n",n+1); 
  }
  

 return 0;
}    
