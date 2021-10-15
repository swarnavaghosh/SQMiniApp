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

     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     
     Get_Input(&lsdft);
     Initialize(&lsdft);

      GetInfluencingAtoms_Local(&lsdft);
      printf("done influencing atoms local\n");

       GetInfluencingAtoms_Nonlocal(&lsdft);
    printf("done influencing atoms nonlocal\n");
    
      CalculateNonlocalProjectors(&lsdft);
    printf("Nonlocal projectors calculated\n");

    // PetscLogDouble t1,t2;
    // PetscTime(&t1);
     clock_t begin = clock();
    NodesAndWeights_LSSGQ(&lsdft);
     clock_t end = clock();
    // PetscTime(&t2);
     double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
     printf("rank=%d, Total time spent=%lf \n",rank,time_spent);

    MPI_Barrier(MPI_COMM_WORLD);
     printf("Done LSSGQ \n");
    
   PetscFinalize();
   
  return 0;
}
