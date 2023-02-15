#include <vector>
#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <random>
#include <chrono>
#include <fstream>
#include <omp.h>

// Decomposition stuff
const int segmentation = 2;
const int N_threads = segmentation*segmentation;
const int N=20; //World Size = NxN = N^2

// Each local grid has sidelength
const int side_length = N / segmentation;
const long int iter = 64*1e5/N_threads;

// Make signal vector indicating whether a core is working or not
std::vector<int> signal(N_threads, 0);

// TEST: collisions counter
int N_collisions = 0;

const double J = 1; // J>0!!!
const double h = 0; // No external magnetic field 
std::mt19937 gen(42);     //Setting the seed for the random number generator



class Lattice {
public:
    // The size of the Lattice
    double Energy;
    double Magnetization;

    // NB: it is up to the calculation to interpret this vector in two dimension.
    
    //std::vector<int> data;
    // NB: this vector has the same length as `spin` and must be interpreted in two dimension as well.
    std::vector<int> data;
    std::vector<double> History_Energy;
    std::vector<double> History_Magnetization;


    Lattice(double Energy, double Magnetization, std::vector<int> data, std::vector<double> History_Energy, std::vector<double> History_Magnetization)         : 
    Energy(Energy), Magnetization(Magnetization), data(std::move(data)), History_Energy(std::move(History_Energy)), History_Magnetization(std::move(History_Magnetization)) {}

};

double checksum(std::vector<double> const& v){
    double sumTotal = 0;
    for(long unsigned int k= (v.size() - 500); k < v.size(); ++k){
        sumTotal += v[k];            
    }
    return sumTotal / 500;
}

void Energy_and_Magnetization(Lattice &lattice) {
    // Calc of the init energy and magnetization
    double E=0;
    double M=0;
    
    //#pragma omp for collapse(2) reduction(+:E,M)
    //#pragma omp taskloop collapse(2) reduction(+:E,M)
    for (int i = 0; i<N; ++i)
    for (int j = 0; j<N; ++j){
        
        int center = lattice.data[i*N+j];
        int upper = lattice.data[((N+i-1)%N)*N+j];
        int lower = lattice.data[((N+i+1)%N)*N+j];
        int right = lattice.data[i*N+((N+j+1)%N)];
        int left = lattice.data[i*N+((N+j-1)%N)];
        
        E += ( -J * center * (upper + lower + right + left) - h * center);
        M += center;
    }
    
    lattice.Energy = E/2;
    lattice.Magnetization = M;
}

void MC_step(Lattice &lattice,double T) {
    std::uniform_int_distribution<> Matrix_index(0, side_length-1);
    std::uniform_real_distribution<> Unif_float(0.0, 1.0);
    
    // Set up private variables
    int center;
    int upper;
    int lower;
    int right;
    int left;
    double init_energy;
    double flipped_energy;
    double E_diff;
    float r;
    
    // Get information about my thread
    int my_thread = omp_get_thread_num();
    div_t divresult;
    divresult = div(my_thread, segmentation);
    int x_thread = divresult.rem;
    int y_thread = divresult.quot;
    
    // Get a x,y on matrix:
    int x = Matrix_index(gen) + x_thread*side_length;
    int y = Matrix_index(gen) + y_thread*side_length;
    
    // Check if index is in boundary
    if (x % side_length == 0 || x % side_length == side_length || y % side_length == 0 || y % side_length == side_length ) {
        
        // Set its own signal to 1
        signal[my_thread] = 1;
        
        // Get neighbouring core indices
        int neighbour_up = (segmentation+x_thread-1)%segmentation * segmentation + y_thread;
        int neighbour_down = (segmentation+x_thread+1)%segmentation * segmentation + y_thread;
        int neighbour_left = x_thread*segmentation + ((segmentation+y_thread-1)%segmentation);
        int neighbour_right = x_thread*segmentation + ((segmentation+y_thread+1)%segmentation);
        
        // Check if neighbouring cores are in boundary
        if (signal[neighbour_up] == 1 || signal[neighbour_down] == 1 || signal[neighbour_left] == 1 || signal[neighbour_right] == 1) {
            
            // TEST COLLISIONS
            N_collisions += 1;
            
            // While in bound
            while (x % side_length == 0 || x % side_length == side_length || y % side_length == 0 || y % side_length == side_length ) {
            
                // Resample
                x = Matrix_index(gen) + x_thread*side_length;
                y = Matrix_index(gen) + y_thread*side_length; 
            }
            
            // Change its own signal back to 0
            signal[my_thread] = 0;
        }
    }
    
    // Calculate current energy:
    center = lattice.data[x*N+y];
    upper = lattice.data[((N+x-1)%N)*N+y];
    lower = lattice.data[((N+x+1)%N)*N+y];
    right = lattice.data[x*N+((N+y+1)%N)];
    left = lattice.data[x*N+((N+y-1)%N)];
    
    init_energy = ( -J * center * (upper + lower + right + left) - h * center);
    flipped_energy = (-1* init_energy);
    
    // Calc difference
    E_diff = (flipped_energy - init_energy);

    r = Unif_float(gen);
    if (exp(-E_diff/T) >= r) {
        lattice.data[x*N+y] = (-1* center);   
        lattice.Energy += E_diff;
        lattice.Magnetization += (-2*center);
    }
}
        

int main(int argc, char **argv) {
    long int iter = 1e5;
    srand(42);     //Setting the seed for the random number generator

    
    std::vector <double> T_vec = {3.27};
    
    // ***********************************************************
    // Initialize Lattice: //
    double Init_Energy = 0;
    double Init_Magnetization = 0  ;
    std::vector <double> history_vec_en;
    std::vector <double> history_vec_mag;

    std::vector<int> Spins;
    for (int i =0; i < N*N; i++){
            Spins.push_back(2*(rand() % 2 )-1);
    }
    // ***********************************************************
    
    
    for (double T: T_vec) {
    
        Lattice lattice = Lattice(Init_Energy, Init_Magnetization,Spins, history_vec_en,history_vec_mag);
        Energy_and_Magnetization(lattice); // Calculating the initial 
        lattice.History_Energy.push_back(lattice.Energy);
        lattice.History_Magnetization.push_back(lattice.Magnetization);

        
        auto begin = std::chrono::steady_clock::now();

        omp_set_num_threads(N_threads);

        // Set up private variables for energy calculation

        #pragma omp parallel reduction(+:N_collisions)
        {
            
            for (long int i = 1; i<iter+1; ++i){
                MC_step(lattice,T);

                #pragma omp single
                {
                    Energy_and_Magnetization(lattice); // Calculating the initial
                    lattice.History_Energy.push_back(lattice.Energy);
                    lattice.History_Magnetization.push_back(lattice.Magnetization);
                }
            } 
        }
        
        // TEST COLLISIONS
        std::cout<< "Total no. collisions"<<N_collisions<<std::endl;

        auto end = std::chrono::steady_clock::now();
        std::cout << "elapsed time  : " << (end - begin).count() / 1000000000.0 << " sec" << std::endl;
        std::cout << "Checksum Energy, (mean of last 500)  : " <<  checksum(lattice.History_Energy) << std::endl;
        std::cout << "Checksum Magnetization, (mean of last 500)  : " <<  checksum(lattice.History_Magnetization) << std::endl;

        std::string S = std::to_string(T).substr(0, std::to_string(T).find(".")+3);
        std::string A = "Results-";
        std::string txt = ".txt";    


        std::ofstream myfile;
        myfile.open(A+S+txt);
        myfile << "History" << ";" << "Magnetiztion" ;
        myfile << "\n" ;
        for (int i = 0; i < iter; ++i) {
            myfile << lattice.History_Energy[i] << ";" << lattice.History_Magnetization[i];
            myfile << "\n" ;
        }
        myfile.close();
    }

    
    return 0;
}