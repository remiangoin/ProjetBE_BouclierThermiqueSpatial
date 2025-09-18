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
    const double L     = 0.01;

    const int Nx = 101;
    const int Nt = 200;
    const int Nmax = 100;

    VectorXd x = VectorXd::LinSpaced(Nx, 0.0, L);
    VectorXd t = VectorXd::LinSpaced(Nt, 0.0, t_end);

    MatrixXd Theta = solve_heat_series(k, rho, Cp, Text, T0, L, x, t, Nmax);

    affichage(Theta, x, t, T0, Text, L);

    return EXIT_SUCCESS;
}

