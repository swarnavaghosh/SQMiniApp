// Swarnava Ghosh
// NCCS, ORNL (ghoshs@ornl.gov)
static char help[] = "Simulation Package for Linear Scaling Ab-initio Real-space Calculations (SPARC-LS) \n\
options:\n\
-name name of file\n";

#include "lsdft.h"
#include "ilsdft.h"
#include <petsctime.h>
#include <mpi.h>
#include <time.h>
//#undef __FUNCT__
//#define __FUNCT__ "main"
///////////////////////////////////////////////////////////////////////////////////////////////
//                              main: the main function                                      //
///////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{

   int ierr;
    LSDFT_OBJ lsdft;
    int rank;
    
    PetscInitialize(&argc, &argv, (char *) 0, help);
    PetscLogDouble v1,v2,elapsed_time;
    PetscLogDouble w1,w2,total_time;
     clock_t begin0 = clock();
     PetscTime(&w1);

     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     
     Get_Input(&lsdft);
     Initialize(&lsdft);

      GetInfluencingAtoms_Local(&lsdft);
      printf("done influencing atoms local\n");

       GetInfluencingAtoms_Nonlocal(&lsdft);
    printf("done influencing atoms nonlocal\n");
    
      CalculateNonlocalProjectors(&lsdft);
    printf("Nonlocal projectors calculated\n");
    
    //MPI_Barrier(MPI_COMM_WORLD);

    clock_t begin = clock();
    PetscTime(&v1);

    //time_t start,end;
    //time(&start);
    NodesAndWeights_LSSGQ(&lsdft);
    //time(&end);
    //double time_spent = double(end - start);
    PetscTime(&v2);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
     printf("rank=%d, Time spent in SQ=%lf \n",rank,time_spent);
     elapsed_time = v2 - v1;   
     printf("rank=%d, Petsc Time spent in SQ=%lf \n",rank,elapsed_time);

    MPI_Barrier(MPI_COMM_WORLD);
     printf("Done LSSGQ \n");
     PetscTime(&w2);
     clock_t end0 = clock();
     double total_time_spent = (double)(end0 - begin0) / CLOCKS_PER_SEC;
     printf("rank=%d, Total time spent=%lf \n",rank,total_time_spent);
     total_time = w2 - w1;
     printf("rank=%d, Petsc Total time spent=%lf \n",rank,total_time);
     PetscFinalize();
    
  return 0;
}
