/*
  |
  | file name: scf.cc
  |
  | Description: This file contains the functions required for self consistent field
  |
  | Authors: Swarnava Ghosh
  |
  
  |-------------------------------------------------------------------------------------------*/
#include "lsdft.h"
#include "petscsys.h"
#include <cmath>
// #include "mkl_lapacke.h"
// #include "mkl.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// PetscScalar SCF_LSSGQ(LSDFT_OBJ *pLsdft)
// {
//    NodesAndWeights_LSSGQ(pLsdft);
//   return 0;
// }
  


///////////////////////////////////////////////////////////////////////////////////////////
// Calculate LSSGQ weights
//////////////////////////////////////////////////////////////////////////////////////////
PetscScalar NodesAndWeights_LSSGQ(LSDFT_OBJ *pLsdft)
{
    // first form a local hamiltonian 
    // then go over all the nodes in the processor domain and calculate electron density using LSSGQ
  PetscInt pos, LIpos, k, j, i,Kmid,Jmid,Imid;
    PetscInt xcor, ycor, zcor, lxdim, lydim, lzdim;
    // AO aodmda1;
    // AO aodmda2;
    // DMDAGetAO(pLsdft->da, &aodmda1);
    // DMDAGetAO(pLsdft->daloc, &aodmda2);
    PetscInt pos2, Nmid;
    int rank;
    PetscScalar ***arrVkm1, ***arrVk;
    int ctr=0,tot;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   
    DMDAGetCorners(pLsdft->da, &xcor, &ycor, &zcor, &lxdim, &lydim, &lzdim);
   
    printf("rank=%d, number of FD points=%d \n",rank,lxdim*lydim*lzdim);
    
    Kmid=(pLsdft->nz_loc-1)/2;
    Jmid=(pLsdft->ny_loc-1)/2;
    Imid=(pLsdft->nx_loc-1)/2;

     Nmid = Kmid *pLsdft->nx_loc * pLsdft->ny_loc + Jmid * pLsdft->nx_loc + Imid;    // this is natural ordering for the central node

     tot=lxdim*lydim*lzdim;
    for (k = zcor; k < zcor + lzdim; k++)
        for (j = ycor; j < ycor + lydim; j++)
            for (i = xcor; i < xcor + lxdim; i++)
    	      {
                 
		//	printf("$ starting LSSGQ at rank=%d, at node (k=%d,j=%d,i=%d) ctr=%d out of %d \n",rank,k,j,i,ctr,tot);

                LIpos = (k - zcor) * lxdim * lydim + (j - ycor) * lxdim + (i - xcor);
	
               	//FormVeffNodal(pLsdft,k,j,i);
		// assign constant values (0.333) for now
		VecSet(pLsdft->VeffNodal,0.333);
                         
                MatCopy(pLsdft->LaplacianOprloc,pLsdft->LapPlusVeffOprloc,SAME_NONZERO_PATTERN);
		//2. Add Veff to laplacian                                                                                                                                            
                MatDiagonalSet(pLsdft->LapPlusVeffOprloc,pLsdft->VeffNodal,ADD_VALUES);
		
                CalculateSpectralNodesAndWeights(pLsdft, Nmid, LIpos,k,j,i);
		ctr++;
	
            }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Calculate LSSGQ weights
//////////////////////////////////////////////////////////////////////////////////////////
PetscScalar CalculateSpectralNodesAndWeights(LSDFT_OBJ *pLsdft, int p, int LIp,int Kp, int Jp, int Ip)
{
    // this is the routine that calculates spectral weights and nodes
    // this routine will do Lanczos

    PetscInt N_qp;              // number of quadrature points
    // PetscInt Nd;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    N_qp = pLsdft->N_qp;
    // Nd = pLsdft->Nd;

    //PetscInt Nx, Ny, Nz;
    //    PetscScalar Lk, Mk, Lkp1, Mkp1, deltaL, deltaM;
    int k;
    // Vec Vk, iVk;
    //Vec Vkm1, iVkm1;
    //Vec Vkp1, iVkp1;
    PetscScalar *a, *b;
      k=0;

      //

      PetscMalloc(sizeof(PetscScalar)*(N_qp+1), &a);
      PetscMalloc(sizeof(PetscScalar)*(N_qp+1), &b);

    //zero out vectors
    VecZeroEntries(pLsdft->Vk);
    VecZeroEntries(pLsdft->Vkm1);
    VecZeroEntries(pLsdft->Vkp1);

    VecSetValue(pLsdft->Vkm1, p, 1.0, INSERT_VALUES);   // index of p might be incorrect
    //     printf("value set in Vkm1, p=%d, LIp=%d \n",p,LIp);
    //Mult_HamiltonianVector(pLsdft, &pLsdft->Vkm1, &pLsdft->Vk,Kp,Jp,Ip);
    MatMult(pLsdft->LapPlusVeffOprloc,pLsdft->Vkm1,pLsdft->Vk);
     
    //     printf("Hamiltonian Vector mult done \n");
    VecDot(pLsdft->Vkm1, pLsdft->Vk, &a[0]);
    //    printf("Here a\n");
    VecAXPY(pLsdft->Vk, -a[0], pLsdft->Vkm1);
    VecNorm(pLsdft->Vk, NORM_2, &b[0]);
    //   printf("Here b\n");
    VecScale(pLsdft->Vk, 1.0 / b[0]);
    //     printf("Here c\n");
   

     //   printf("rank=%d, a[%d]=%lf,b[%d]=%lf \n",rank,k,a[k],k,b[k]);      

    for (k = 0; k < N_qp; k++) {
      //Mult_HamiltonianVector(pLsdft, &pLsdft->Vk, &pLsdft->Vkp1,Kp,Jp,Ip);  
      MatMult(pLsdft->LapPlusVeffOprloc,pLsdft->Vk,pLsdft->Vkp1);  
      VecDot(pLsdft->Vk, pLsdft->Vkp1, &a[k + 1]);
        VecAXPY(pLsdft->Vkp1, -a[k + 1], pLsdft->Vk);
        VecAXPY(pLsdft->Vkp1, -b[k], pLsdft->Vkm1);
        VecCopy(pLsdft->Vk, pLsdft->Vkm1);
        VecNorm(pLsdft->Vkp1, NORM_2, &b[k + 1]);
        VecCopy(pLsdft->Vkp1, pLsdft->Vk);
        VecScale(pLsdft->Vk, 1.0 / b[k + 1]);

	//	printf("rank=%d, a[%d]=%lf,b[%d]=%lf \n",rank,k,a[k],k,b[k]);      
    }
   
    /*
     * Call function to find eigenvalue of Tridiagonal matrix
     *
     */
    TridiagEigenVecSolve_NodesAndWeights(pLsdft, a, b, N_qp, LIp);  // todo this function

    // free a,b
    PetscFree(a);
    PetscFree(b);


    return 0;
}


