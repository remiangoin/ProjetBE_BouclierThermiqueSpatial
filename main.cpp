#include <Eigen/Dense>
#include <iostream>
#include <cstdio>
#include <thread>
#include <chrono>
#include <algorithm>
#include <cmath>


// Pour éviter d'écrire Eigen:: partout
using Eigen::VectorXd;
using Eigen::MatrixXd;

// --- série de Fourier (Dirichlet) ---
MatrixXd solve_heat_series(
    double k, double rho, double Cp,
    double Text, double T0,
    double L,
    const VectorXd& x,
    const VectorXd& t,
    int Nmax
) {
    const double PI = 3.14159265358979323846;
    const double alpha = k / (rho * Cp);

    MatrixXd Theta(t.size(), x.size());
    Theta.setConstant(T0); // remplissage initial

    for (int j = 0; j < t.size(); ++j) {
        for (int i = 0; i < x.size(); ++i) {
            double sum = 0.0;
            for (int n = 0; n <= Nmax; ++n) {
                double m = 2.0*n + 1.0;
                double kn = (m*PI) / (2.0*L);
                sum += (1.0/m) * std::sin(kn*x[i]) * std::exp(-alpha*kn*kn*t[j]);
            }
            Theta(j,i) = Text + (4.0*(T0-Text)/PI)*sum;
        }
    }
    return Theta;
}


MatrixXd construction_matrices_1D(
    double k, double rho, double Cp)

    return K;
}
// --- conduction 1D par volumes finis (Euler implicite) ---
// Hypothèses: k, rho, Cp constants. CL: Dirichlet T=Text aux deux bords.
// x peut être non uniforme; on interprète x[i] comme centres de cellules (y compris bords).
MatrixXd solve_heat_Constant_volumeMethod(
    double k, double rho, double Cp,
    double Text, double T0,
    double L,
    const VectorXd& x,
    const VectorXd& t,
    int /*Nmax non utilisé en VF*/
) {
    const int Nx = static_cast<int>(x.size());
    const int Nt = static_cast<int>(t.size());
    MatrixXd Theta(Nt, Nx);

    // --- paramètres physiques
    const double alpha = k / (rho * Cp); // (pas utilisé directement, mais utile pour checks)

    // --- géométrie locale : tailles de cellules et distances entre centres
    // Cette méthode permet de gérer un maillage non uniforme. (où le pas dx n est pas variable et doit être calculé)
    VectorXd dx_cell(Nx);        // Δx_i (volume de contrôle)
    VectorXd d_face(Nx - 1);     // d_{i+1/2} = x[i+1]-x[i]
    for (int i = 0; i < Nx - 1; ++i) d_face[i] = x[i + 1] - x[i];
    dx_cell[0]       = 0.5 * (x[1] - x[0]);
    dx_cell[Nx - 1]  = 0.5 * (x[Nx - 1] - x[Nx - 2]);
    for (int i = 1; i < Nx - 1; ++i) dx_cell[i] = 0.5 * (x[i + 1] - x[i - 1]);

    // --- transmissibilités diffusives aux faces : T_{i+1/2} = k / d_{i+1/2}
    VectorXd Tface(Nx - 1);
    for (int i = 0; i < Nx - 1; ++i) Tface[i] = k / d_face[i];

    // --- "transmissibilités" de bord (demi-cellule) pour Dirichlet
    const double TL = Text, TR = Text;
    const double Tbord_L = k / (dx_cell[0]);           // = 2k/(x[1]-x[0]) si x régulier aux bords
    const double Tbord_R = k / (dx_cell[Nx - 1]);      // = 2k/(x[N-1]-x[N-2])

    // --- matrices d'assemblage (denses pour simplicité ; tri-diagonales en pratique) (on peut sinon utiliser des sparse)
    MatrixXd K = MatrixXd::Zero(Nx, Nx);     // conductance (diffusion)
    VectorXd b_bord = VectorXd::Zero(Nx);    // source due aux CL

    // Intérieur : K tri-diagonale
    for (int i = 1; i < Nx - 1; ++i) {
        K(i, i - 1) += -Tface[i - 1];
        K(i, i)     +=  (Tface[i - 1] + Tface[i]);
        K(i, i + 1) += -Tface[i];
    }

    // Bords Dirichlet via demi-cellule : ajout sur diagonale + RHS
    // Bord gauche (i = 0)
    if (Nx >= 2) {
        // flux entre 0 et 1 (face 1/2)
        K(0, 0)     +=  Tface[0];
        K(0, 1)     += -Tface[0];
        // lien avec le "réservoir" TL via demi-cellule
        K(0, 0)     +=  Tbord_L;
        b_bord[0]   +=  Tbord_L * TL;
    } else {
        // Cas dégénéré Nx==1
        K(0,0)      +=  Tbord_L + Tbord_R;
        b_bord[0]   +=  Tbord_L*TL + Tbord_R*TR;
    }

    // Bord droit (i = Nx-1)
    if (Nx >= 2) {
        // flux entre N-2 et N-1 (face N-1/2)
        K(Nx - 1, Nx - 2) += -Tface[Nx - 2];
        K(Nx - 1, Nx - 1) +=  Tface[Nx - 2];
        // lien avec "réservoir" TR via demi-cellule
        K(Nx - 1, Nx - 1) +=  Tbord_R;
        b_bord[Nx - 1]    +=  Tbord_R * TR;
    }

    // --- masse diagonale
    VectorXd Mdiag(Nx);
    for (int i = 0; i < Nx; ++i) Mdiag[i] = rho * Cp * dx_cell[i];

    // --- état initial : T0 partout mais on force les bords à Text (cohérent Dirichlet)
    VectorXd Tn = VectorXd::Constant(Nx, T0);
    Tn[0] = TL; Tn[Nx - 1] = TR;
    Theta.row(0) = Tn.transpose();

    // --- boucle en temps (Euler implicite : (M/Δt + K) T^{n+1} = M/Δt T^n + b_bord)
    for (int n = 1; n < Nt; ++n) {
        const double dt = std::max(1e-12, t[n] - t[n - 1]);

        // Assemble A = M/dt + K
        MatrixXd A = K;
        for (int i = 0; i < Nx; ++i) A(i, i) += Mdiag[i] / dt; //Therme de gauche que multiplie T^{n+1}

        // RHS
        VectorXd rhs = b_bord;
        for (int i = 0; i < Nx; ++i) rhs[i] += (Mdiag[i] / dt) * Tn[i]; //Therme complet de droite

        // Résolution (A est SPD pour diffusion pure + Dirichlet) : LDLT/LLT conviennent
        Eigen::LLT<MatrixXd> llt(A);
        VectorXd Tnp1 = llt.solve(rhs); // Résolution du système linéaire via une méthode spéciale pour matrice SPD

        // Par sécurité (numérique), on impose exactement Dirichlet aux bords
        Tnp1[0] = TL; Tnp1[Nx - 1] = TR;

        Theta.row(n) = Tnp1.transpose(); // sauvegarde du vecteur Tnp1(x) dans une matrice T(x,t)
        Tn.swap(Tnp1);  // Tn = Tnp1 pour la prochaine itération (swap plus efficace que affectation)
    }

    return Theta;
}


// --- affichage avec gnuplot ---
int affichage(const MatrixXd& Theta, const VectorXd& x, const VectorXd& t,
              double T0, double Text, double L) {
    const int Nt = Theta.rows();
    const int Nx = Theta.cols();

    FILE* gp = _popen("gnuplot -persist", "w");
    if (!gp) { std::perror("gnuplot"); return EXIT_FAILURE; }

    double ymin = std::min(T0, Text);
    double ymax = std::max(T0, Text);
    double pad = 0.05*(ymax - ymin);
    ymin -= pad; ymax += pad;

    std::fprintf(gp, "set term wxt size 900,600\n");
    std::fprintf(gp, "set title 'Theta(x,t) — evolution temporelle'\n");
    std::fprintf(gp, "set xlabel 'x (m)'\nset ylabel 'Temperature (K)'\n");
    std::fprintf(gp, "set xrange [0:%g]\n", L);
    std::fprintf(gp, "set yrange [%g:%g]\n", ymin, ymax);
    std::fprintf(gp, "set grid\n");

    for (int j = 0; j < Nt; ++j) {
        std::fprintf(gp, "unset label 1\n");
        std::fprintf(gp, "set label 1 sprintf('t = %.2f s', %f) at graph 0.02, 0.95 front\n", t[j], t[j]);

        std::fprintf(gp, "plot '-' with lines lw 2 title 'profil a t'\n");
        for (int i = 0; i < Nx; ++i)
            std::fprintf(gp, "%g %g\n", x[i], Theta(j,i));
        std::fprintf(gp, "e\n");
        std::fflush(gp);

        std::this_thread::sleep_for(std::chrono::milliseconds(60));
    }

    std::fprintf(gp, "pause -1 'Fin. Appuyez sur Entrée...'\n");
    std::fflush(gp);
    _pclose(gp);
    return EXIT_SUCCESS;
}

int main() {
    // Données
    const double Text  = 1273.15;
    const double T0    = 150.0;
    const double rho   = 191.0;
    const double Cp    = 538.7;
    const double k     = 0.06457;
    const double t_end = 15.0*60.0;
    const double L     = 0.03;

    const int Nx = 101;
    const int Nt = 200;
    const int Nmax = 100;

    VectorXd x = VectorXd::LinSpaced(Nx, 0.0, L);
    VectorXd t = VectorXd::LinSpaced(Nt, 0.0, t_end);

    MatrixXd Theta = solve_heat_series_Constant_volumeMethod(k, rho, Cp, Text, T0, L, x, t, Nmax);

    affichage(Theta, x, t, T0, Text, L);

    return EXIT_SUCCESS;
}
