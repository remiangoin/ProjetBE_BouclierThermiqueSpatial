#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void) {
    // Données (SI)
    double Text  = 1273.15;       // K
    double T0    = 150.0;         // K
    double Tstar = 700; // K  cible : 150°C
    double rho   = 191;        // kg/m^3
    double Cp    = 538.7;        // J/(kg.K) 
    double k     = 0.06457;       // W/(m.K)
    double t     = 15.0 * 60.0;   // s
    double h     = Resolution1DParaConstants(k,rho,Cp,Tstar,Text,T0,t); // m 
    return EXIT_SUCCESS;
}


int Resolution1DAnalytique(k,rho,Cp,Tstar,Text,T0,t){

}

int Resolution1DParaConstants(k,rho,Cp,Tstar,Text,T0,t) {
   
}

int Resolution1DParaVariaMethodeDecouplee(k,rho,Cp,Tstar,Text,T0,t) {

}

