//-------------------------------------CAPTION + LOG file--------------------------------------------------

// This is a program written for CMB lensing - Cosmic Shear Project
// This program calculate the correlated angular power spectrum of CMB lensing and shear.
// Farbod Hassani & Shant Baghram

// *Version:1  *Date:25 April  *Filename:crosslens-25Apr2015-SB.cpp *Edited by:Shant Baghram
// *Version:2  *Date:19 June 2015  *Filename:crosslens-19AJune2015-FH.cpp *Edited by:Farbod Hassani
// *Version:3  *Date:20 June 2015  *Filename:crosslens-20AJune2015-FH.cpp *Edited by:Farbod Hassani Convergence power spectrum is added.


//------------------------------------LIBRARIES-------------------------------------------------

// In this section the required libraries are added by their headers name.


#include "Cosmology-Subroutine.h"   // This is a sub-routine written by me for cosological functions
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "iostream"
#include "stdlib.h"
#include "time.h"
#include <complex>
#include "cmath"
#include "cstdlib"
#include <string>
using namespace std;




//-------------------------------Defines ---------------------------------------------------------

#define PI 3.141592


//-------------------------------Main Program-------------------------------------------------



int main ()

{






//****************************** Parameter Definition *******************

double Omega_m,Omega_r,Omega_L,Omega_cdm,Omega_b, Omega_k, Omega_nu, A_inf, ns,  kpivot, sigma8 , h , Apresent, A_m, Neutrinonum, ksigma, neff,CC,w_DE;  // Cosmological parameters


// **** PLANCK Cosmological parameters *****
//Ref: Planck 2015 results. XIII. Cosmological parameters - arXiv:1502.01589 - TT,TE,EE+lowP+lensing best parameters

// h=0.6727; // Hubble parameter  h= H_0 / 100
h=0.6727; // Simulation

Omega_cdm=0.226;   //Cold Dark matter density parameter is fixed
Omega_b= 0.045;      // The baryon density parameter is fixed
// Omega_m=Omega_b+Omega_cdm;      // Matter  desity parameter value set : Best value
Omega_m=0.271;      // Matter  desity parameter value set : Best value

Omega_L=1. - Omega_m;    // Cosmological constant value set

//Omega_b=0.4/square(h);
// Omega_nu=0.0062/(square(h));   // Neutrino density
// Neutrinonum=3.046; //Number of degenerate massive neutrino species
Omega_r=0.;
Omega_k=0.000;

A_inf=1*2.26990e-9;    // COBE normalization best parameter
kpivot=0.05;       // Pivot point wavenumber in Mpc^{-1}
ns=0.966;         //primordial power-spectrum's spectral index
 A_m=(8*square(PI)*A_inf*pow(kpivot,1.0-ns)*pow(3000,4)*Inverse(h))/(25*square(Omega_m));// Matter Fluctuation
 // A_m=0;// Matter Fluctuation

    // w_DE=-1;   //LCDM model
sigma8=0.8150;      // Sigma - 8 Value from data
// Apresent= 2.1*1e6;  // present time
// Theta_CMB=2.7255;
// zCMB=1090;
// //normkCMB=Inverse(0.345259);
// ksigma= 0.306;               // The wavenumber which has a unity variance "sigma" in LCDM model
// neff=-1.550;                 //Effective spectral index in LCDM
// CC=0.384;                    //Effective second derivative of power-spectrum
// double bias;
// double lnorm;
// bias=1.9;

    double Sigma1;
    Sigma1= Sigmaa_BBKS(8/h, Omega_r,  Omega_cdm,   Omega_b ,   Omega_L , A_m,  h,  ns, 0);
    cout<<" Sigma= "<< sqrt(Sigma1)<<endl;




    return(0);

    }
