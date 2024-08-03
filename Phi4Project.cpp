#include <iostream>
#include <fstream>   //used for writing data
#include <sstream>   //used for reading csv data
#include <cmath>     //used for useful map functions
#include <vector>    //used for storing field vals and data in a vector
#include <random>    //used for random generator
#include <chrono>    //used for seeding random genrator
#include<algorithm>  //used for min function
#include <numeric>   //accumulate function 
#include <iterator>
#include <iomanip>
#include <chrono>    // For high-resolution clock

using namespace std;

/* Data structure to store all the parameters of the algorithms */
struct hmc_params_t {
    double tlength{ 1 }; /* trajectory length */
    int nstep{ 10 };        /* leapfrog steps per trajectory */
    int ntherm{ 100 };      /* number of thermalization steps about 10/20% of total steps */
    int ntraj{ 20000 };        /* number of trajectories after thermalization */
};

struct act_params_t {
    /*
    double kappa{ 0.185825 };
    double lambda{ 1.1689 };
    */
    double kappa{ 0.400173 / 2.0 };
    double lambda{ 0.5 };
    int d{ 3 };
};

static hmc_params_t hmc_params;
static act_params_t act_params;

class Hopf {
private:
    int D_;
    int L_;
    int V_;

public:
    vector<vector<int>> hop;
    Hopf(int D, int L) : D_(D), L_(L), V_(static_cast<int>(pow(L, D))) {
        hop.resize(V_, vector<int>(2 * D_));

        int x, y, Lk;
        int xk, k, dxk;

        /* Go through all the points */
        for (x = 0; x < V_; x++) {
            Lk = V_;
            y = x;

            /* Go through the components k */
            for (k = D_ - 1; k >= 0; k--) {
                /* Begin iteration from top dimension, successively scanning through layers
                   of the lattice like reading a book - hence lexicographically */

                Lk /= L_;                       /* pow(L, k) */
                xk = y / Lk;                    /* kth component */
                y = y - xk * Lk;                /* y <- y % Lk */

                /* Forward - the first D components */
                if (xk < L_ - 1) { dxk = Lk; }
                else { dxk = Lk * (1 - L_); }
                hop[x][k] = x + dxk;

                /* Backward - the last D components */
                if (xk > 0) { dxk = -Lk; }
                else { dxk = Lk * (L_ - 1); }
                hop[x][k + D_] = x + dxk;
            }
        }
    }
};

void fillWithGaussian(vector<double>& vec, mt19937& gen) {
    // Function to fill a vector with Gaussian random numbers
    normal_distribution<double> dist(0, 1.0);

    generate(vec.begin(), vec.end(), [&]() {
        return dist(gen);
        });
}

void fillWithUniform(vector<double>& vec) {
    // Function to fill a vector with Gaussian random numbers
    static random_device rd;
    static mt19937 gen(rd());
    static uniform_real_distribution<double> dist(-1, 1);

    generate(vec.begin(), vec.end(), [&]() {
        return dist(gen);
        });
}

double CalculateAverage(const vector<double>& vec) {
    if (vec.empty()) {
        throw invalid_argument("Vector is empty. Cannot calculate average.");
    }

    // Calculate the sum of the elements
    double sum = std::accumulate(vec.begin(), vec.end(), 0.0);

    // Calculate and return the average
    return sum / vec.size();
}

class Phi4 {
private:
    int D_;
    int L_;
    int V_;

public:
    vector<double> phi;
    vector<double> mom;
    vector<double> dummy_mom;

    vector<vector<int>> hop;
    mt19937 gen;

    vector<vector<int>> indices_t;

    Phi4(int D, int L) : D_(D), L_(L), V_(static_cast<int>(pow(L, D))) {
        phi.assign(V_, 0);

        Hopf hopf(D_, L_);
        hop = hopf.hop;

        mom.assign(V_, 0);
        dummy_mom.assign(V_, 0);

        /* 
        random_device rd; //seed
        mt19937 gen(10001);
        */
    }

    double Magnetisation() {
        /* Loop over all sites */
        double m = 0;
        for (int i = 0; i < V_; i++) {
            m += phi[i];
        }
        return m;
    }

    double TwoPtSusc() {
        return phi[0] * Magnetisation();
    }

    double G_0() {
        return pow(Magnetisation(), 2);
    }
   
    /*
    std::vector<int> index_to_coordinates(int i) {
        std::vector<int> coordinates(D_);
        std::vector<int> dimensions(D_, L_);


        for (int d = 0; d < D_; ++d) {
            int product = 1;
            for (int k = 0; k < d; ++k) {
                product *= dimensions[k];
            }
            coordinates[d] = (i / product) % dimensions[d];
        }

        return coordinates;
    }
    */

    vector<int> index_to_coordinates(int i) { 
        vector<int> coordinates(D_);
        int product;

        for (int d = 0; d < D_; ++d) {
            product = pow(L_, d);
            coordinates[d] = (i / product) % L_ ;
        }

        return coordinates;
    }

    int coordinates_to_index(const std::vector<int>& coordinates) {
        int index = 0;
        int product=1;

        for (int d = 0; d < D_; ++d) {
            index += coordinates[d] * product;
            product = pow(L_, d);
        }

        return index;
    }
    
    void Spatial_coords() {

        vector<vector<int>> coords;

        coords.resize(L_, vector<int>(V_/L_));
        vector<int> index_to_sum_at_t(L_);


        vector<int> coordinates(D_);

        for (int t = 0; t < L_; t++) {

            index_to_sum_at_t.resize(V_ / L_, 0);
            /* Go through all the spacial points */
            for (int x = 0; x < L_; x++) {
                for (int y = 0; y < L_; y++) {
                    coordinates = { t, x, y };
                    index_to_sum_at_t.push_back(coordinates_to_index(coordinates));
                }
            }
            coords[t] = index_to_sum_at_t;
        }

        indices_t = coords;
    }

    vector<double> G_t() { // Specified to D = 2+1 spacetime dimensions
        vector<double> G(L_, 0);
        vector<int> coordinates(D_);
        vector<int> indices;
        
        for (int t = 0; t < L_; t++) {
            indices = indices_t[t];
            /* Go through all the spacial points */
            for (int i = 0; i < indices.size(); i++) {
                G[t] += phi[0]*phi[indices[i]];
            }
        }

        return G;
    }

    double Correlations() {
        /* Loop over all sites */
        double W = 0;
        double phin;

        for (int i = 0; i < V_; i++) {
            /* Sum over neighbors in positive direction */
            phin = 0;
            for (int mu = 0; mu < D_; mu++) {
                phin += phi[hop[i][mu]];
            }
            W += phi[i] * phin;
        }
        return 2 * W;
    }

    double Action() {
        double phin, S, phi2;
        double kappa = act_params.kappa;
        double lambda = act_params.lambda;

        S = 0;

        /* Loop over all sites */
        for (int i = 0; i < V_; i++) {
            /* Sum over neighbors in positive direction */
            phin = 0;
            for (int mu = 0; mu < D_; mu++) {
                phin += phi[hop[i][mu]];
            }
            phi2 = phi[i] * phi[i];
            S += -2 * kappa * phin * phi[i] + phi2 + lambda * (phi2 - 1.0) * (phi2 - 1.0);
        }

        return S;
    }

    double Hamiltonian() {
        double H = 0;

        /* Loop over all sites */
        for (int i = 0; i < V_; i++) {
            H += pow(mom[i], 2);
        }
        H = 0.5 * H + Action();
        return H;  // Add missing return statement
    }

    void Move_phi(double eps) {
        for (int i = 0; i < V_; i++) {
            phi[i] += eps * mom[i];
        }
    }

    double Force(int i) {
        int mu;

        double phin, f;
        double kappa = act_params.kappa;
        double lambda = act_params.lambda;

        /* Sum over neighbors in positive direction and negative directions */
        double Jx = 0;
        for (mu = 0; mu < D_; mu++) {
            Jx += phi[hop[i][mu]] + phi[hop[i][D_ + mu]];
        }

        f = 2 * kappa * Jx - 2 * phi[i] - 4 * lambda * (phi[i] * phi[i] - 1.0) * phi[i];

        //cout << "Force: " << f << endl;

        return f;
    }

    void Move_mom(double eps) {
        for (int i = 0; i < V_; i++) {
            mom[i] += eps * Force(i);
        }
    }

    void MomentumRotation(double eps, double gamma) {
        fillWithGaussian(dummy_mom, gen);
        double c1 = exp(-gamma * eps);
        double c2 = pow(1-c1*c1, 0.5);

        for (int i = 0; i < V_; i++) {
            mom[i] = c1*mom[i] + c2* dummy_mom[i];
        }
    }

    void Leapfrog(double tau, int nstep) {
        double eps{ tau / nstep };
        for (int j = 1; j < nstep; j++) {
            Move_phi(eps / 2);
            Move_mom(eps);
            Move_phi(eps / 2);
        }
    }

    void Omelyan(double tau, int nstep) {
        double eps{ tau / nstep };
        double zeta = 0.1931833;

        for (int j = 1; j < nstep; j++) {
            Move_phi(eps * zeta);
            Move_mom(eps / 2.0);
            Move_phi(eps * (1 - 2 * zeta));
            Move_mom(eps / 2.0);
            Move_phi(eps * zeta);
        }
    }

    void Omelyan4(double tau, int nstep) {
        double eps{ tau / nstep };
        double r1{ 0.08398315262876693 }, r2{ 0.2539785108410595 }, r3{ 0.6822365335719091 }, r4{ -0.03230286765269967 };

        for (int j = 1; j < nstep; j++) {
        Move_mom(eps * r1);
        Move_phi(eps * r2);
        Move_mom(eps * r3);
        Move_phi(eps * r4);
        Move_mom(eps * (1 / 2 - r1 - r3));
        Move_phi(eps * (1 - 2 * (r2 + r4)));
        Move_mom(eps * (1 / 2 - r1 - r3));
        Move_phi(eps * r4);
        Move_mom(eps * r3);
        Move_phi(eps * r2);
        Move_mom(eps * r1);
        }
    }

    void HybridDynamics() {
        int ntraj = hmc_params.ntraj;
        int nstep = hmc_params.nstep;
        int emm = 100;

        double eps = 0.1;
        double tau = emm * eps;

        for (int i = 0; i < ntraj; i++) {
            fillWithGaussian(mom, gen);
            Leapfrog(tau, nstep);
        }
    }

    double GenerateUniformRandom() {
        // Function to generate a uniform random number in the range [0, 1)
        static random_device rd; // Seed for the random number engine
        static mt19937 gen(rd()); // Mersenne Twister random number engine
        static uniform_real_distribution<> dis(0.0, 1.0); // Uniform distribution in the range [0, 1)
        return dis(gen); // Generate the random number
    }

    void SaveFinalConfiguration() {
        ofstream file;
        file.open("FinalConfig.csv");
        //file << "phi" << "," << "mom" << endl;
        for (int i = 0; i < phi.size(); i++) {
            file << phi[i] << "," << mom[i] << endl;
        }
        file.close();
    }

    void LoadFinalConfiguration() {
        vector <double> vec1{};
        vector <double> vec2{};


        // Open the file in read mode
        ifstream inFile("FinalConfig.csv");

        // Check if the file is opened successfully
        if (!inFile.is_open()) {
            cerr << "Failed to open file: " << "FinalConfig.csv" << endl;
            return;
        }

        string line;
        while (getline(inFile, line)) {
            stringstream ss(line);
            string item;
            double val1, val2;

            // Read first value
            if (getline(ss, item, ',')) {
                val1 = stod(item);
            }
            else {
                cerr << "Error reading first value in line: " << line << endl;
                continue;
            }

            // Read second value
            if (getline(ss, item, ',')) {
                val2 = stod(item);
            }
            else {
                std::cerr << "Error reading second value in line: " << line << std::endl;
                continue;
            }

            // Add values to vectors
            vec1.push_back(val1);
            vec2.push_back(val2);
        }

        // Close the file
        inFile.close();

        phi = vec1;
        mom = vec2;
    }

    void SaveRNGState() {

        std::ofstream ofs("RNGstate.dat", std::ios::out | std::ios::binary);
        if (ofs) {
            ofs << gen;
            ofs.close();
        }
        else {
            std::cerr << "Error opening file for writing: " << "RNGstate.dat" << std::endl;
        }
    }

    void LoadRNGState() {
        std::ifstream ifs("RNGstate.dat", std::ios::in | std::ios::binary);
        if (ifs) {
            ifs >> gen;
            ifs.close();
        }
        else {
            std::cerr << "Error opening file for reading: " << "RNGstate.dat" << std::endl;
        }
    }

    vector<double> generateEvenlySpaced(int L) {
        std::vector<double> vec(L);
        double step = 1 / (L - 1);

        for (int i = 0; i < L; ++i) {
            vec[i] = i * step;
        }

        return vec;
    }

    void TestLeapfrog() {
        double tau = 1;
        int L = 1000;
        vector<double> eps_vec = generateEvenlySpaced(L);

        ofstream file;
        file.open("dHdata.csv");
        file << "dH" << "," << "eps" << endl;

        fillWithGaussian(mom, gen);
        vector<double> momorig{ mom };

        fillWithGaussian(phi, gen);
        vector<double> phiorig{ phi };

        double H0{ Hamiltonian() };

        for (int i = 2; i < L; i += 1) {
            phi = phiorig;
            mom = momorig;

            double eps = eps_vec[i];
            Leapfrog(0.1, i);
            double dH{ Hamiltonian() - H0 };

            file << dH << "," << (0.1 / i) << endl;

        }
        file.close();

    }

    void BinderCumulant(const vector<double> mag) {
        // Square all elements using std::transform
        vector<double> mag4 = mag;
        transform(mag4.begin(), mag4.end(), mag4.begin(), [](double x) {
            return std::pow(x, 4);
            });

        double bincum = CalculateAverage(mag4) / pow(CalculateAverage(mag4), 2);

        cout << "Binder Cumulant: " << bincum << endl;
    }

    bool Acceptance(double &delta) {
        if (delta <= 0) {
            return true;
        }
        else {
            if (exp(-delta) > GenerateUniformRandom()) {
                //cout << "accept" << endl;
                return true;
            }
            else {
                //cout << "reject" << endl;
                return false;
            }
        }
    }

    void HybridMonteCarlo() {
        int ntraj = hmc_params.ntraj;
        int ntherm = hmc_params.ntherm;
        int nstep = hmc_params.nstep;
        double tau = hmc_params.tlength;

        /* Observables */
        double M, M2, chi2, H0;
        double expdH{ 0 }, accept_rate{ 0.0 }, dH{ 0 }, Mag_average{ 0 };
        int accpept_count{ 0 };


        /* File management */
        ofstream file;
        file.open("LatticeDATA.csv");

        //file.open("MagDATA.csv");
        //file << "M, " << "M2, " << "Chi2, " << "E, " << "dH " << endl;
        vector <double> mag_untherm;

        //fillWithUniform(phi);

        for (int i = 1; i < ntherm; i++) {
            fillWithGaussian(mom, gen);

            vector<double> phiold{ phi };
            H0 =  Hamiltonian() ;

            Omelyan(tau, nstep);

            dH = Hamiltonian() - H0 ;

            if (!Acceptance(dH)) {
                phi = phiold;
            }
            else {
                accpept_count += 1;
                /* Observables */
                /*
                M = Magnetisation();
                chi2 = phi[0] * M;
                M2 = M * M;
                file << M << "," << M2 << "," << chi2 << "," << H0 << "," << dH << endl;
                */
            }
        }

        vector <double> mag;
        int reject_count = 0;
        for (int i = 1; i < ntraj; i++) {
            fillWithGaussian(mom, gen);

            vector<double> phiold{ phi };
            H0 = Hamiltonian();

            Omelyan(tau, nstep);

            dH = Hamiltonian() - H0 ;

            if (!Acceptance(dH)) {
                phi = phiold;
            }
            else {
                accpept_count += 1;
                /* Observables */
                /*
                M = Magnetisation();
                chi2 = phi[0] * M;
                M2 = M * M;
                expdH += exp(-dH);
                Mag_average = Mag_average + M;
                file << M << "," << M2 << "," << chi2 << "," << H0 << "," << dH << endl;
                */
                for (int j = 0; j < phi.size(); j++) { file << phi[j] << ","; }
                file << endl;
                
            }
        }
        /*
        expdH = expdH / ntraj;
        Mag_average = Mag_average / ntraj;
        accept_rate = double(accpept_count) / double(ntraj+ntherm);
        cout << "<M>= " << Mag_average << endl;
        cout << "exp(dH) " << expdH << endl;
        cout << "Acceptance rate: " << accept_rate << endl;
        */

        file.close();
    }

    void SMD() {
        int ntraj = hmc_params.ntraj;
        int ntherm = hmc_params.ntherm;
        int nstep = hmc_params.nstep;
        double tau = hmc_params.tlength;
        double eps{ tau / nstep };
        double gamma = 0.16;

        /* Observables */
        double M, M2, chi2, H0;
        double expdH{ 0 }, accept_rate{ 0.0 }, dH{ 0 }, Mag_average{ 0 };
        int accpept_count{ 0 };


        /* File management */
        ofstream file;
        file.open("MagDATA_smd.csv");
        file << "M, " << "M2, " << "Chi2, " << "E, " << "dH " << endl;
        vector <double> mag_untherm;

        //fillWithUniform(phi);

        vector <double> mag;
        int reject_count = 0;
        for (int i = 1; i < ntraj; i++) {
            MomentumRotation(eps, gamma);

            vector<double> phiold{ phi };

            H0 = Hamiltonian();

            Omelyan(tau, nstep);

            dH = Hamiltonian() - H0;

            if (!Acceptance(dH)) {
                phi = phiold;
            }
            else {
                if (i % 4 == 0) {
                    accpept_count += 1;
                    /* Observables */
                    M = Magnetisation();
                    chi2 = phi[0] * M;
                    M2 = M * M;
                    expdH += exp(-dH);
                    Mag_average = Mag_average + M;

                    file << M << "," << M2 << "," << chi2 << "," << H0 << "," << dH << endl;
                }
            }
        }

        expdH = expdH / ntraj;
        Mag_average = Mag_average / ntraj;
        accept_rate = double(accpept_count) / double(ntraj + ntherm);


        cout << "SMD Algorithm " << endl;
        cout << "<M>= " << Mag_average << endl;
        cout << "exp(dH) " << expdH << endl;
        cout << "Acceptance rate: " << accept_rate << endl;


        file.close();
    }

    void HybridMonteCarlo_forCorrellator() {
        int ntraj = hmc_params.ntraj;
        int ntherm = hmc_params.ntherm;
        int nstep = hmc_params.nstep;
        double tau = hmc_params.tlength;

        /* Observables */
        double H0;
        double expdH{ 0 }, accept_rate{ 0.0 }, dH{ 0 }, Mag_average{ 0 };
        int accpept_count{ 0 };
        vector<double> G;

        /* Thermailsation Steps */
        for (int i = 1; i < ntherm; i++) {
            fillWithGaussian(mom, gen);

            vector<double> phiold{ phi };
            H0 = Hamiltonian();

            Leapfrog(tau, nstep);

            dH = Hamiltonian() - H0;

            if (!Acceptance(dH)) {
                phi = phiold;
            }
        }

        cout<<"Finshed thermalising"<<endl;

        /* File management */
        ofstream file;
        file.open("CorrelatorDATA.csv");
        
        /* Find spacial coords to sum over */
        Spatial_coords();

        int reject_count = 0;
        for (int i = 1; i < ntraj; i++) {
            fillWithGaussian(mom, gen);

            vector<double> phiold{ phi };
            H0 = Hamiltonian();

            Leapfrog(tau, nstep);

            dH = Hamiltonian() - H0;

            if (!Acceptance(dH)) {
                phi = phiold;
            }
            else {     
                accpept_count += 1;
                expdH += exp(-dH);

                G = G_t();
                for (int t = 0; t < L_; t++) { file << G[t] << ", "; } file << endl; 
            }
        }

        expdH = expdH / ntraj;
        
        accept_rate = double(accpept_count) / double(ntraj + ntherm);

        cout << "exp(dH) " << expdH << endl;
        cout << "Acceptance rate: " << accept_rate << endl;

        file.close();
    }

    void SMD_forCorrellator() {
        int ntraj = hmc_params.ntraj;
        int ntherm = hmc_params.ntherm;
        int nstep = hmc_params.nstep;
        double tau = hmc_params.tlength;
        double eps{ tau / nstep };
        double gamma = 0.16;

        /* Observables */
        double H0;
        double expdH{ 0 }, accept_rate{ 0.0 }, dH{ 0 }, Mag_average{ 0 };
        int accpept_count{ 0 };
        vector<double> G;

        /* File management */
        ofstream file;
        file.open("CorrelatorDATA.csv");

        //fillWithUniform(phi);

        vector <double> mag;
        int reject_count = 0;
        for (int i = 1; i < ntraj; i++) {
            MomentumRotation(eps, gamma);

            vector<double> phiold{ phi };

            H0 = Hamiltonian();

            Omelyan(tau, nstep);

            dH = Hamiltonian() - H0;

            if (!Acceptance(dH)) {
                phi = phiold;
            }
            else {
                accpept_count += 1;
                expdH += exp(-dH);

                G = G_t();
                for (int t = 0; t < L_; t++) { file << G[t] << ", "; } file << endl;
            }
        }

        
        accept_rate = double(accpept_count) / double(ntraj);


        cout << "SMD Algorithm " << endl;
        cout << "exp(dH) " << expdH << endl;
        cout << "Acceptance rate: " << accept_rate << endl;


        file.close();
    }
};

int main() {    
    auto start = std::chrono::high_resolution_clock::now(); // start the timer

    Phi4 phi4(3, 16);

    //phi4.LoadFinalConfiguration();
    //phi4.LoadRNGState();

    phi4.HybridMonteCarlo();
    //phi4.HybridMonteCarlo_forCorrellator();

    //phi4.SaveRNGState();
    //phi4.SaveFinalConfiguration();

    auto end = chrono::high_resolution_clock::now(); // Stop the timer and measure how long
    chrono::duration<double, ratio<60>> duration = end - start;
    cout << "Time taken to run: " << duration.count() << " mins" << std::endl;
    
    return 0;
}