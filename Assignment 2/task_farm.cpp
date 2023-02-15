/*
  Assignment: Make an MPI task farm. A "task" is a randomly generated integer.
  To "execute" a task, the worker sleeps for the given number of milliseconds.
  The result of a task should be send back from the worker to the master. It
  contains the rank of the worker
*/

#include <iostream>
#include <random>
#include <chrono>
#include <thread>
#include <array>

// To run an MPI program we always need to include the MPI headers
#include <mpi.h>

/* Define constants and random seed to make it reproducible */
const int NTASKS=5000;  // number of tasks 
const int RANDOM_SEED=1234;

/* Define the master function, the one who manages the tasks */
void master (int nworker) {
    
    /* two arrays that we define in one line, static size so not exactly like vectors */
    std::array<int, NTASKS> task, result;

    // set up a random number generator
    std::random_device rd;
    std::default_random_engine engine(RANDOM_SEED);
    
    // make a distribution of random integers in the interval [0:30]
    std::uniform_int_distribution<int> distribution(0, 30);
    
    /* Select one random number that represents the task.*/
    for (int& t : task) {
        t = distribution(engine);   // set up some "tasks"
    }
    
    /* Loop over workers to send out first tasks */
    for (int i = 0; i < nworker; i++ ) {
        
        /* Define task to send */
        int xsend = task[i];
            
        /* Send out a task to each worker */ 
        //data_pointer, count, datatype, tag =index, communicator
        MPI_Send(&xsend, 1, MPI_INT, i+1, i, MPI_COMM_WORLD); 
    }
    
    /* Define how long we are in the task list */
    int task_to_send = nworker;
    
    /* Looping as long as there are tasks to complete */
    while (task_to_send < NTASKS) { //IT WAS <= BEFORE
        
        /* Receive messages */
        //data_pointer, count, datatype, source, tag, comm, status
        int xrecv;
        MPI_Status status;
        MPI_Recv(&xrecv, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        
        /* Extract tag and source */    
        int tag = status.MPI_TAG;
        int source = status.MPI_SOURCE;
        
        /* Store results */
        result[tag] = xrecv;
        
        /* Define message to send */
        int xsend = task[task_to_send];
        
        /* Send message */
        MPI_Send(&xsend, 1, MPI_INT, source, task_to_send, MPI_COMM_WORLD);
        
        /* +1 to task index thing to complete */
        task_to_send++;
    }
    
    /* Receive the last tasks from the workers (since we increment task_to_send before we receive)*/
    for (int i = 0; i < nworker; i++ ) { 
        
        /* Receive */
        int xrecv;
        MPI_Status status;
        MPI_Recv(&xrecv, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        
        /* Extract tag and source */    
        int tag = status.MPI_TAG;
        int source = status.MPI_SOURCE;
        
        /* Store results */
        result[tag] = xrecv;
        
        /* Send signal to stop. 1 = work is done, 0 = more work to do */
        int xsend = 0; 
        MPI_Send(&xsend, 1, MPI_INT, source, NTASKS, MPI_COMM_WORLD); //SIGNAL CHANGED TO NTASKS
        
    }

    // Print out a status on how many tasks were completed by each worker
    for (int worker=1; worker<=nworker; worker++) {
        int tasksdone = 0; int workdone = 0;
        for (int itask=0; itask<NTASKS; itask++)
        if (result[itask]==worker) {
            tasksdone++;
            workdone += task[itask];
        }
        std::cout << "Master: Worker " << worker << " solved " << tasksdone << 
                    " tasks\n";    
    }
}

// call this function to complete the task. It sleeps for task milliseconds
void task_function(int task) {
    std::this_thread::sleep_for(std::chrono::milliseconds(task));
}

void worker (int rank) {
    
    /* Define signal to 0, there is work to do when beginning */
    int SIGNAL = 0;
    
    /* Loop as long as signal is zero */
    while (SIGNAL == 0) {
        
        /* Receive message from master */
        int xrecv;
        MPI_Status status;
        MPI_Recv(&xrecv, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); //SET TAG TO NTASK FOR SIGNAL LATER
        
        /* Extract tag */
        int tag = status.MPI_TAG;
        
        if (tag == NTASKS){
            return; 
        }
    
        /* Perform task */
        task_function(xrecv);
        
        /* Send result (rank) to master */
        MPI_Send(&rank, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
    
    // CHANGE TO IRECV AND ISEND, WHILE SIGNAL MAYBE WITH TAG FROM MASTER TO WORKER
        
    }   
}

int main(int argc, char *argv[]) {
    /* allocate memory for them */
    int nrank, rank;
    
    /* argc stuff is for synchronizing arguments on all ranks */
    MPI_Init(&argc, &argv);                // set up MPI
    /* save rank and nrank value at allocated memory */
    MPI_Comm_size(MPI_COMM_WORLD, &nrank); // get the total number of ranks
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // get the rank of this process
    
    /* The same code runs on multiple computers, so first we need to check if the computer it runs on is a master or worker */
    if (rank == 0)       // rank 0 is the master
        master(nrank-1); // there is nrank-1 worker processes
    else                 // ranks in [1:nrank] are workers
        worker(rank);

    MPI_Finalize();      // shutdown MPI
}
