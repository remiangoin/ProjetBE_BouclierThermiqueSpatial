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

    // Diffusivité
    double alpha = k / (rho * Cp);

    // Facteur A (doit être 0 < A < 1)
    double A = ((Tstar - Text) * M_PI) / (4.0 * (T0 - Text));

    if (A <= 0.0 || A >= 1.0) {
        fprintf(stderr, "Erreur: A doit être dans (0,1). Vérifie Text, T0, Tstar.\n");
        return EXIT_FAILURE;
    }

    // Solution fermée pour h
    double h = (M_PI / 2.0) / sqrt(-(1.0 / (alpha * t)) * log(A));

    // Affichage
    printf("alpha = %.3e m^2/s\n", alpha);
    printf("A      = %.6f\n", A);
    printf("h      = %.6f m  (%.2f cm)\n", h, 100.0 * h);

    return EXIT_SUCCESS;
}

