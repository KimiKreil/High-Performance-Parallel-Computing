#include <vector>
#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <random>
#include <chrono>
#include <fstream>



const int N=20; //World Size = NxN = N^2
const double J = 1; // J>0!!!
const double h = 0; // No external magnetic field 
std::mt19937 gen(42);     //Setting the seed for the random number generator
const long int iter = 1e6;

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
    for (int i = 0; i<N; ++i)
    for (int j = 0; j<N; ++j){
        
        int center = lattice.data[i*N+j];
        int upper = lattice.data[((N+i-1)%N)*N+j];
        int lower = lattice.data[((N+i+1)%(N))*N+j];
        int right = lattice.data[i*N+((N+j+1)%N)];
        int left = lattice.data[i*N+((N+j-1)%N)];
        
        
        E += ( -J * center * (upper + lower + right + left) - h * center);
        M += center;
    }

    lattice.Energy = E/2;
    lattice.Magnetization = M;
}

void MC_step(Lattice &lattice,double T) {
    std::uniform_int_distribution<> Matrix_index(0, N-1);
    std::uniform_real_distribution<> Unif_float(0.0, 1.0);

    
        
    // Get a x,y on matrix:
    int x = Matrix_index(gen);
    int y = Matrix_index(gen);
    

    // Calculate current energy:
    int center = lattice.data[x*N+y];
    int upper = lattice.data[((N+x-1)%N)*N+y];
    int lower = lattice.data[((N+x+1)%N)*N+y];
    int right = lattice.data[x*N+((N+y+1)%N)];
    int left = lattice.data[x*N+((N+y-1)%N)];
    
    
    double Init_Energy = ( -J * center * (upper + lower + right + left) - h * center);
    double Flipped_Energy = (-1* Init_Energy);
    
    // Calc difference
    double E_Diff = (Flipped_Energy - Init_Energy);

    float r = Unif_float(gen);

    if (exp(-E_Diff/T) > r) {
        lattice.data[x*N+y] = (-1* center);   
        lattice.Energy += E_Diff;
        lattice.Magnetization += (-2*center);
        
    }
}
        


int main(int argc, char **argv) {
    auto begin = std::chrono::steady_clock::now();

    srand (42);     //Setting the seed for the random number generator

    
    std::vector <double> T_vec = {2.27};

    
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

    
   
    for (long int i = 1; i<iter+1; ++i){
        MC_step(lattice,T);
        lattice.History_Energy.push_back(lattice.Energy);
        lattice.History_Magnetization.push_back(lattice.Magnetization);
    }
          
  
    
    
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
    for (int i = 0; i < iter; ++i,i += 1) {
        myfile << lattice.History_Energy[i] << ";" << lattice.History_Magnetization[i];
        myfile << "\n" ;}
    
    myfile.close();
    }

    
    return 0;
}