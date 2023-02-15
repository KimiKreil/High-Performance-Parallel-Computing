
/* This is our assignment for the SIR model */


/* Including various libraries */
#include <iostream>
#include <vector>
#include <fstream>

/* Defining constants */
const double N=1000;
const double beta = 0.2;
const double gamma = 1.0/10.0;

/* Define number of days to simulate */
const double N_days=200;

/* Define time steps to test */
std::vector<double> Times = {1.0,0.5,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001}; //{1.0,0.1,0.01,0.001,0.0001};

/* Define a function calculating next step of S_(n+1) */
void update_vectors( std::vector<double> &S, std::vector<double> &I, std::vector<double> &R, double t_step  ) {
    
    /* Take the previous S value, aka the last element in the vector. 
    Same for I and R */
    double S_now = S.back();
    double I_now = I.back();
    double R_now = R.back();
    
    /* Discretized equation for dS ------------------------ */
    double S_next = S_now - ( beta * I_now * (S_now/N) )*t_step;
        
    /* Append S_next to the S vector */
    S.push_back(S_next);
    
    
    /* Discretized equation for dI ------------------------ */
    double I_next = I_now + ( beta * I_now * (S_now/N) - gamma * I_now )*t_step;
    
    /* Append I_next to the I vector */
    I.push_back(I_next);
    
    
    /* Discretized equation for dR ------------------------ */
    double R_next = R_now + ( gamma * I_now )*t_step;
    
    /* Append I_next to the I vector */
    R.push_back(R_next);  
}

/* Update the vector N_ite times */
int main() {
    
    /* Define name of file with results */
    std::ofstream myfile("./sir_t_analysis.txt");
    
    /* Define the header */
    myfile << "t_step N_ite S I R \n";
    
    /* Make a loop over the t_steps we will use */
    for (double t_step : Times) {
        
        /* Initial parameters */
        std::vector<double> I = {1.0};
        std::vector<double> R = {0.0};
        std::vector<double> S = {N-I.back()};
        
        /* Get number of iterations */
        int N_ite = N_days/t_step;
        
        /* Make a for loop to update the vectors and export values to file */
        for (int i = 0; i < N_ite; ++i) {  
            update_vectors(S, I, R, t_step);
        }   
        
        /* Insert to file the last value of S, I, R */
            myfile << t_step << ' '
                   << N_ite << ' '
                   << S.back() << ' '
                   << I.back() << ' '
                   << R.back() << '\n'; 
    }
    
    return 0;
}

