#include <mpi.h>
#include <stdio.h>

int sum_even(int i, int max);

int main(int argc, char** argv) {

  // Initiallise MPI environment
  MPI_Init(NULL,NULL);

  // Get processor rank
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int count = rank + 1;

  if(rank == 0) {
    MPI_Send(&count, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  } else if(rank == 1) {
    MPI_Recv(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
      MPI_STATUS_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);
  } 


  if(rank == 0) { 
    int sum = sum_even(0, 20);
           
    printf("The sum of even numbers between 0 and 20 is %d.\n",sum);
  }

  // Finalise MPI environment
  MPI_Finalize();

  return 0;

}


int sum_even(int i, int max) {

  int result = 0;

  while(i <= max) {
    if(i % 2 == 0){
      result += i;
    }
    i++;
  }
  
  return result;
}

