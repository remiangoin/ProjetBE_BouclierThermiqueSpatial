#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <thread>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include <iostream>


using Matrix = std::vector<std::vector<double>>;

// --- util ---
std::vector<double> linspace(double a, double b, std::size_t n) {
    std::vector<double> v(n);
    if (n == 1) { v[0] = a; return v; }
    const double h = (b - a) / double(n - 1);
    for (std::size_t i = 0; i < n; ++i) v[i] = a + h * double(i);
    return v;
}

// --- série de Fourier (Dirichlet) ---
Matrix solve_heat_series(
    double k, double rho, double Cp,
    double Text, double T0,
    double L,
    const std::vector<double>& x_grid,
    const std::vector<double>& t_grid,
    int Nmax
) {
    const double PI = 3.14159265358979323846;
    const double alpha = k / (rho * Cp);
    Matrix Theta(t_grid.size(), std::vector<double>(x_grid.size(), Text));
    for (std::size_t j = 0; j < t_grid.size(); ++j) {
        const double t = t_grid[j];
        for (std::size_t i = 0; i < x_grid.size(); ++i) {
            const double x = x_grid[i];
            double sum = 0.0;
            for (int n = 0; n <= Nmax; ++n) {
                const double m = 2.0*n + 1.0;
                const double kn = (m*PI) / (2.0*L);
                sum += (1.0/m) * std::sin(kn*x) * std::exp(-alpha*kn*kn*t);
            }
            Theta[j][i] = Text + (4.0*(T0-Text)/PI)*sum;
        }
    }
    return Theta;
}
int affichage(Matrix Theta, const double T0, const double Text, const double L) {
    const std::size_t Nt = Theta.size();
    const std::size_t Nx = Theta[0].size();
    auto x = linspace(0.0, L, Nx);
    auto t = linspace(0.0, 15.0*60.0, Nt);
 // --- ANIMATION avec gnuplot ---
    // Ouvre un pipe vers gnuplot
    FILE* gp = _popen("gnuplot -persist", "w"); // sous Windows
    if (!gp) { std::perror("gnuplot"); return EXIT_FAILURE; }

    // Axes fixes
    double ymin = std::min(T0, Text);
    double ymax = std::max(T0, Text);
    // (optionnel) marge visuelle
    double pad = 0.05*(ymax - ymin);
    ymin -= pad; ymax += pad;

    std::fprintf(gp, "set term wxt size 900,600\n");
    std::fprintf(gp, "set title 'Theta(x,t) — evolution temporelle'\n");
    std::fprintf(gp, "set xlabel 'x (m)'\nset ylabel 'Temperature (K)'\n");
    std::fprintf(gp, "set xrange [0:%g]\n", L);
    std::fprintf(gp, "set yrange [%g:%g]\n", ymin, ymax);
    std::fprintf(gp, "set grid\n");

    for (std::size_t j = 0; j < Nt; ++j) {
        // Affiche l’instant
        std::fprintf(gp, "unset label 1\n");
        std::fprintf(gp, "set label 1 sprintf('t = %.2f s', %f) at graph 0.02, 0.95 front\n", t[j], t[j]);

        // Envoie les données de la courbe x -> Theta[j][x]
        std::fprintf(gp, "plot '-' with lines lw 2 title 'profil a t'\n");
        for (std::size_t i = 0; i < Nx; ++i)
            std::fprintf(gp, "%g %g\n", x[i], Theta[j][i]);
        std::fprintf(gp, "e\n");
        std::fflush(gp);

        // cadence (ms)
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

    const std::size_t Nx = 101;
    const std::size_t Nt = 200;   // plus d’images => plus fluide
    const int Nmax = 100;

    auto x = linspace(0.0, L, Nx);
    auto t = linspace(0.0, t_end, Nt);
    auto Theta = solve_heat_series(k, rho, Cp, Text, T0, L, x, t, Nmax);

    affichage(Theta, T0, Text, L);

    return EXIT_SUCCESS;
}
