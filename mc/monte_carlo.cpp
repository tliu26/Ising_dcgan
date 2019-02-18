#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

using namespace std;

int main(int argc, char *argv[])
{
    /* Check input file */
    if (argc != 2)
    {
        cerr << "Please specify an input file" << endl;
        exit(1);
    }
    string inFileName = argv[1];
    cout << "Reading from " << inFileName << " ..." << endl << endl;

    /* Initialize parameters */
    int Nx, Ny, N, Niter;
    double J, T, beta, W4, W2;
    string outdir, J_str, T_str;

    /* Read input file and calculate parameters */
    ifstream inFile(inFileName);
    string line;
    if (inFile.is_open())
    {
        while (getline(inFile, line))
        {
            istringstream iss(line);
            string var_name, var_v;
            iss >> var_name >> var_v;
            if (var_name == "Nx") Nx = stoi(var_v);
            if (var_name == "Ny") Ny = stoi(var_v);
            if (var_name == "J")
            {
                J = stod(var_v);
                J_str = var_v;
            }
            if (var_name == "T")
            {
                T = stod(var_v);
                T_str = var_v;
            }
            if (var_name == "Niter") Niter = stoi(var_v);
            if (var_name == "Outdir") outdir = var_v;
        }
    }
    else
    {
        cerr << "Unable to open input file" << endl << endl;
        exit(1);
    }
    struct stat sb;
    if (stat(outdir.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode))
    {
        cerr << "Output directory not valied" << endl << endl;
        exit(1);
    }
    
    /* calculate some parameters */
    N = Ny * Nx;
    beta = 1/T;
    W4 = exp(-beta * 2 * J * 4);
    W2 = exp(-beta * 2 * J * 2);

    /* set simulation variables */
    int S[N], total_spin;
    int xi, yi;
    int xi_p, xi_m, yi_p, yi_m;
    double rnd;
    double dE, prob, avg_spin;
    int s_sums;
    int naccept = 0;

    cout.precision(10);

    /* random number generator */
    random_device rd;
    unsigned int rn = rd();
    mt19937 mt(rn);
    uniform_real_distribution<double> dis1(0.0, 1.0);
    uniform_int_distribution<> dis2x(0, Nx-1);
    uniform_int_distribution<> dis2y(0, Ny-1);

    /* open files for writing output */
    string dirName = outdir + "/" + "N" + to_string(Ny) + "x" + to_string(Nx) + "_J"
    + J_str + "_T" + T_str + "_Niter" + to_string(Niter) + "_rn" + to_string(rn);
    mkdir(dirName.c_str(), 0744);
    ofstream S_configs_f(dirName + "/S_configs.dat", ios::binary); // lattice configuration at each ite
    ofstream avg_spin_f(dirName + "/avg_spin.dat", ios::binary); // average spin at each ite
    cout << "Writing to " << dirName << " ..." << endl << endl;

    /* initialize lattice and count total spins */
    total_spin = 0;
    for (int i=0; i<N; i++)
    {
        rnd = dis1(mt);
        if (rnd < 0.5) S[i] = 1;
        else S[i] = -1;
        total_spin += S[i];
    }

    /* monte carlo */
    for (int ite=0; ite<Niter; ite++)
    {
        for (int i=0; i<N; i++)
        {
            /* select a random lattice site */
            xi = dis2x(rd);
            yi = dis2y(rd);

            /* indices of nearest neighbors with pbc */
            xi_p = (xi      + 1) % Nx;
            xi_m = (xi + Nx - 1) % Nx;
            yi_p = (yi      + 1) % Ny;
            yi_m = (yi + Ny - 1) % Ny;

            /* calculate change in energy given by 2J * S_ij * (sum of spins over n.n.'s) */
            s_sums = S[yi * Nx + xi] *
            (S[yi * Nx + xi_p] + S[yi * Nx + xi_m] + S[yi_p * Nx + xi] + S[yi_m * Nx + xi]);
            if (s_sums == 2) prob = W2;
            else if (s_sums == 4) prob = W4;
            if (s_sums <= 0 || dis1(mt) < prob)
            {
                S[yi * Nx + xi] = -1 * S[yi * Nx + xi];
                naccept ++;

                /* update quantities */
                dE = 2 * J * s_sums;
                total_spin += 2 * S[yi * Nx + xi];
            }
        }
        /* measurements */
        avg_spin = (double) total_spin / N;

        /* write to file */
        avg_spin_f.write((char *) &avg_spin, sizeof(avg_spin));
        S_configs_f.write((char *) &S, sizeof(S));

        if (ite%500 == 0)
        {
            cout << "ite = " << ite << endl;
            cout << avg_spin << endl;
        }
    }

    cout << endl << "naccept = " << (double) naccept/Niter << endl;
}