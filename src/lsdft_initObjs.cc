/*=============================================================================================
  |
  | file name: initObjs.cc
  |
  | Description: This file contains the functions required for initializing variables
  |
  | Authors: Swarnava Ghosh,
  |
  | 
  |-------------------------------------------------------------------------------------------*/
#include "lsdft.h"
//#include "sddft.h"
#include "petscsys.h"
#include <petsctime.h> 
// #include "mkl_lapacke.h"
// #include "mkl.h" 
#include <iostream>
using namespace std;
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define max(a,b,c) (a>b&&a>c?a:b>c?b:c)
///////////////////////////////////////////////////////////////////////////////////////////////
//                                Get_Input: reads the input filename                        //
///////////////////////////////////////////////////////////////////////////////////////////////
void Get_Input(LSDFT_OBJ *pLsdft)
{
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-name", pLsdft->file, sizeof(pLsdft->file), PETSC_NULL);
    //PetscOptionsGetString(PETSC_NULL,"-name",pLsdft->file,sizeof(pLsdft->file),PETSC_NULL);
    //PetscOptionsGetString(PetscOptions options,const char pre[],const char name[],char string[],size_t len,PetscBool  *set)
}

///////////////////////////////////////////////////////////////////////////////////////////////
//                              Initialize: calls functions related to initialization        //
///////////////////////////////////////////////////////////////////////////////////////////////
void Initialize(LSDFT_OBJ *pLsdft)
{
    Read_Parameters(pLsdft);
    printf("parameters read \n");
    Read_Ion(pLsdft);
    printf("ions read \n");   
    Read_Relax(pLsdft);
    printf("relaxparameters read \n");
    Read_Pseudopotential(pLsdft);
    printf("done1 \n");
    Create_Objects(pLsdft);
    printf("done2 \n");
   
    // Calculate_IntWts(pLsdft);
    // printf("Integration weights calculated \n");
   
    Create_Laplacian(pLsdft);   // laplacian
    printf("done3 \n");
    //  MatView(pLsdft->laplaceOpr,PETSC_VIEWER_STDOUT_WORLD);
    Create_LaplacianLocal(pLsdft);  // laplacian for local calculation
    // //printf("Local laplacian created \n");
    //  Create_Gradient(pLsdft);

    // // int rank;
    // // MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    //  if(pLsdft->precond==1)
    //    {
    // 	 Create_HelmholtzOperator(pLsdft);
    // 	 Initialize_HelmholtzSolver(pLsdft);
    //    }
    
    Calculate_SplineFit(pLsdft);
    printf("Spline fit done \n");
    //Create_Orbitals(pLsdft);

   
    Calculate_NonlocalCoefficients(pLsdft);
    printf("Nonlocal coeff done \n");
  
    Calculate_PseudochargeCutoff(pLsdft);
    printf("Pseudocharge cutoff calculation done\n");

     // /*
    //  * Boundary conditions
    //  */
    //  // create supercell of atoms
    //  Calculate_BoundaryInformation(pLsdft);
      Create_SupercellAtomList(pLsdft);
     printf("supercell created");
       Calculate_NonlocalIndex(pLsdft);
    printf("Nonlocal index done \n");
    // Create_VeffNodalBoundary(pLsdft);

    //  exit(1);
    

}

///////////////////////////////////////////////////////////////////////////////////////////////
//                         Create_Objects: Creates the PETSc objects                         //
///////////////////////////////////////////////////////////////////////////////////////////////
PetscErrorCode Create_Objects(LSDFT_OBJ *pLsdft)
{
    PetscInt n_x = pLsdft->numPoints_x;
    PetscInt n_y = pLsdft->numPoints_y;
    PetscInt n_z = pLsdft->numPoints_z;

    PetscInt SCn_x,SCn_y,SCn_z;

    PetscInt nx_loc=pLsdft->nx_loc;
    PetscInt ny_loc=pLsdft->ny_loc;
    PetscInt nz_loc=pLsdft->nz_loc;

    PetscInt o = pLsdft->order;
    int MAX_ITS_ANDERSON = pLsdft->MixingHistory;
    PetscInt xcor, ycor, zcor, lxdim, lydim, lzdim;
    int i, a;
    Mat A;
    PetscMPIInt comm_size;
    PetscScalar *velMag, mass_sum = 0.0;
    double s, x, y, vsum_x = 0.0, vsum_y = 0.0, vsum_z = 0.0;
    int Index = 0;

    AO aodmdaloc;

    MPI_Comm_size(PETSC_COMM_WORLD, &comm_size);

        /*
         * Periodic boundary condition
         */
        if (comm_size == 1) {
            DMDACreate3d(PETSC_COMM_SELF, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                DMDA_STENCIL_STAR, n_x, n_y, n_z, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, o, 0, 0, 0, &pLsdft->da);
        } else {

	  DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                DMDA_STENCIL_STAR, n_x, n_y, n_z, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, o, 0, 0, 0, &pLsdft->da);
	    }

	   DMSetFromOptions(pLsdft->da);
	   DMSetUp(pLsdft->da);

        DMDAGetCorners(pLsdft->da, &xcor, &ycor, &zcor, &lxdim, &lydim, &lzdim);
	SCn_x=lxdim+ 2*ceil(pLsdft->Rcut/pLsdft->delta_x);
	SCn_y=lydim+ 2*ceil(pLsdft->Rcut/pLsdft->delta_y);
	SCn_z=lzdim+ 2*ceil(pLsdft->Rcut/pLsdft->delta_z);

	//printf("SCn_x=%d\n",SCn_x);

	    // create a local dmda for veff_local
	//  DMDACreate3d(PETSC_COMM_SELF,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DMDA_STENCIL_STAR,n_x,n_y,n_z,
        //   PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,o,0,0,0,&pLsdft->daVeff);

    // local dmda for operators and vectors for LSSGQ
      DMDACreate3d(PETSC_COMM_SELF, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_STAR,
            nx_loc, ny_loc, nz_loc, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, o, 0, 0, 0, &pLsdft->daloc);
      
      // local dmda for local supercell
      DMDACreate3d(PETSC_COMM_SELF, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_STAR,
            SCn_x, SCn_y, SCn_z, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, o, 0, 0, 0, &pLsdft->daloccell);
      
       DMSetFromOptions(pLsdft->daloc);
	   DMSetUp(pLsdft->daloc);

	    DMSetFromOptions(pLsdft->daloccell);
	   DMSetUp(pLsdft->daloccell);
      


    // allocate memory for storing processor domains
       pLsdft->ProcDomains = (PetscInt *)calloc(comm_size*7,sizeof(PetscInt));
      pLsdft->OwnershipRange = (PetscInt *)calloc(comm_size*2,sizeof(PetscInt));
      //printf("before ao\n");
      DMDAGetAO(pLsdft->da,&pLsdft->aodmda1);
      //printf("after ao\n");

    DMCreateGlobalVector(pLsdft->da, &pLsdft->elecDensRho);
    //printf("edens created \n");
    VecDuplicate(pLsdft->elecDensRho, &pLsdft->SuperposAtRho);
    VecDuplicate(pLsdft->elecDensRho, &pLsdft->chrgDensB);
    VecDuplicate(pLsdft->elecDensRho, &pLsdft->chrgDensB_TM);
    VecDuplicate(pLsdft->elecDensRho, &pLsdft->potentialPhi);
    VecDuplicate(pLsdft->elecDensRho, &pLsdft->Phi_c);
    VecDuplicate(pLsdft->elecDensRho, &pLsdft->twopiRhoPB);
    VecDuplicate(pLsdft->elecDensRho, &pLsdft->twopiBTMmBPS);
    VecDuplicate(pLsdft->elecDensRho, &pLsdft->Veff);
    VecDuplicate(pLsdft->elecDensRho, &pLsdft->bjVj);
    VecDuplicate(pLsdft->elecDensRho, &pLsdft->bjVj_TM);
    VecDuplicate(pLsdft->elecDensRho, &pLsdft->Vxc);
    VecDuplicate(pLsdft->elecDensRho, &pLsdft->tempVec);
    VecDuplicate(pLsdft->elecDensRho, &pLsdft->PoissonRHSAdd);
    VecDuplicate(pLsdft->elecDensRho, &pLsdft->IntWts);

    VecDuplicate(pLsdft->elecDensRho, &pLsdft->xkprev);
    VecDuplicate(pLsdft->elecDensRho, &pLsdft->xk);
    VecDuplicate(pLsdft->elecDensRho, &pLsdft->fkprev);

    //DMCreateLocalVector(pLsdft->daVeff, &pLsdft->Veff_loc);
    //DMCreateGlobalVector(pLsdft->daVeff, &pLsdft->Veff_loc);

    // printf("SG1\n");
    DMCreateGlobalVector(pLsdft->daloc, &pLsdft->Vk);
    // DMCreateLocalVector(pLsdft->daloc, &pLsdft->Vk);
    // printf("SG2\n");
    VecDuplicate(pLsdft->Vk, &pLsdft->Vkm1);
    VecDuplicate(pLsdft->Vk, &pLsdft->Vkp1);
    VecDuplicate(pLsdft->Vk, &pLsdft->VeffVktemp);
    VecDuplicate(pLsdft->Vk, &pLsdft->VeffNodal);
    VecDuplicate(pLsdft->Vk, &pLsdft->DMcolgrad_x);
    VecDuplicate(pLsdft->Vk, &pLsdft->DMcolgrad_y);
    VecDuplicate(pLsdft->Vk, &pLsdft->DMcolgrad_z);
    VecDuplicate(pLsdft->Vk, &pLsdft->DMcol);
    VecDuplicate(pLsdft->Vk, &pLsdft->Tk);
    VecDuplicate(pLsdft->Vk, &pLsdft->Tkp1);
    VecDuplicate(pLsdft->Vk, &pLsdft->Tkm1);
    
    // VecDuplicate(pLsdft->Vk, &pLsdft->VeffNodalBoundary);//todo
    //  DMCreateGlobalVector(pLsdft->daloccell, &pLsdft->VeffNodalBoundary);
    // printf("before gboundary\n");
      DMCreateLocalVector(pLsdft->daloccell, &pLsdft->VeffNodalBoundary);
      //  printf("after gboundary\n");

    // create vec scatter context
    //  VecScatterCreate(pLsdft->Veff,NULL,pLsdft->Veff_loc,NULL,&pLsdft->VeffScatterCtx);
    DMDAGlobalToNaturalAllCreate(pLsdft->da, &pLsdft->VeffScatterCtx);
    //VecScatterCreateToAll(pLsdft->Veff,&pLsdft->VeffScatterCtx,&pLsdft->Veff_loc);
    //  printf("here1 \n");

    PetscMalloc(sizeof(Vec) * (MAX_ITS_ANDERSON), &pLsdft->Xk);
    PetscMalloc(sizeof(Vec) * (MAX_ITS_ANDERSON), &pLsdft->Fk);
    PetscMalloc(sizeof(Vec) * (MAX_ITS_ANDERSON), &pLsdft->XpbF);
    //  printf("here3 pLsdft->MixingHistory=%d \n", pLsdft->MixingHistory);
    for (i = 0; i < MAX_ITS_ANDERSON; i++) {
      //   printf("Inside loop \n");
        VecDuplicate(pLsdft->elecDensRho, &pLsdft->Xk[i]);
	//  printf("Ai=%d", i);
        VecDuplicate(pLsdft->elecDensRho, &pLsdft->Fk[i]);
	// printf("Bi=%d", i);
        VecDuplicate(pLsdft->elecDensRho, &pLsdft->XpbF[i]);
	// printf("Ci=%d", i);
    }
    // printf("here3a \n");

    if (comm_size == 1) {
        DMCreateMatrix(pLsdft->da, &pLsdft->laplaceOpr);    // for poisson
        DMSetMatType(pLsdft->da, MATSEQBAIJ);
    } else {
      //  printf("a1\n");
        DMCreateMatrix(pLsdft->da, &pLsdft->laplaceOpr);
        //printf("a2\n");
        DMSetMatType(pLsdft->da, MATMPIBAIJ);
	// printf("here3b \n");
    }
    // printf("here4 \n");
    DMCreateMatrix(pLsdft->daloc, &pLsdft->LaplacianOprloc);
    // create local gradient operators
    DMCreateMatrix(pLsdft->daloc,&pLsdft->gradient_x);
      DMCreateMatrix(pLsdft->daloc,&pLsdft->gradient_y);
      DMCreateMatrix(pLsdft->daloc,&pLsdft->gradient_z); 
    DMSetMatType(pLsdft->daloc, MATSEQBAIJ);
    
    A = pLsdft->laplaceOpr;

    //  if(pLsdft->kptParalFlag==1)
    //  KSPCreate(pLsdft->kpoint_group_comm,&pLsdft->ksp);
    //else
    KSPCreate(PETSC_COMM_WORLD, &pLsdft->ksp);

    KSPSetType(pLsdft->ksp, KSP_TYPE);
    KSPSetOperators(pLsdft->ksp, A, A);
    KSPSetTolerances(pLsdft->ksp, pLsdft->KSPTOL, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(pLsdft->ksp);
    // printf("here5 \n");

    /*
     * allocate memory for lssgq weights and nodes
     */
    //DMDAGetCorners(pLsdft->da, &xcor, &ycor, &zcor, &lxdim, &lydim, &lzdim);
    // printf("lxdim=%d,lydim=%d,lzdim=%d \n", lxdim, lydim, lzdim);
    PetscMalloc(sizeof(PetscScalar *) * (lxdim * lydim * lzdim), &pLsdft->lssgq_weights);
    PetscMalloc(sizeof(PetscScalar *) * (lxdim * lydim * lzdim), &pLsdft->lssgq_lambda);
    assert(pLsdft->lssgq_weights != NULL && pLsdft->lssgq_lambda != NULL);
    for (i = 0; i < (lxdim * lydim * lzdim); i++) {
        PetscMalloc(sizeof(PetscScalar) * pLsdft->N_qp, &pLsdft->lssgq_weights[i]);
        PetscMalloc(sizeof(PetscScalar) * pLsdft->N_qp, &pLsdft->lssgq_lambda[i]);
    }

    pLsdft->lambda_max = (PetscScalar *)calloc((lxdim * lydim * lzdim),sizeof(PetscScalar));
    pLsdft->lambda_min = (PetscScalar *)calloc((lxdim * lydim * lzdim),sizeof(PetscScalar));

    pLsdft->ChebCoeff = (PetscScalar *)calloc((pLsdft->NplCC+1),sizeof(PetscScalar));

    // pLsdft->lssgq_weights= (PetscScalar **)malloc(sizeof(PetscScalar *)*(lxdim*lydim*lzdim));
    // pLsdft->lssgq_weights= (PetscScalar **)malloc(sizeof(PetscScalar *)*(lxdim*lydim*lzdim));
    // assert(pLsdft->lssgq_weights != NULL && pLsdft->lssgq_lambda != NULL);

    //  printf("here2 \n");

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//                     Create_Laplacian: Initializes the -(1/2)*Laplacian operator           //
///////////////////////////////////////////////////////////////////////////////////////////////
PetscErrorCode Create_Laplacian(LSDFT_OBJ *pLsdft)
{
    PetscInt i, j, k, l, colidx, gxdim, gydim, gzdim, xcor, ycor, zcor, lxdim, lydim, lzdim;
    MatStencil row;
    MatStencil *col;
    PetscScalar *val;
    Mat A = pLsdft->laplaceOpr;
    PetscInt o = pLsdft->order;
#if _DEBUG
    PetscTruth flg;
#endif

    DMDAGetInfo(pLsdft->da, 0, &gxdim, &gydim, &gzdim, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    DMDAGetCorners(pLsdft->da, &xcor, &ycor, &zcor, &lxdim, &lydim, &lzdim);

    MatScale(pLsdft->laplaceOpr, 0.0);
    //MatScale(pLsdft->HamiltonianOpr,0.0);

    PetscMalloc(sizeof(MatStencil) * (o * 6 + 1), &col);
    PetscMalloc(sizeof(PetscScalar) * (o * 6 + 1), &val);

    for (k = zcor; k < zcor + lzdim; k++)
        for (j = ycor; j < ycor + lydim; j++)
            for (i = xcor; i < xcor + lxdim; i++) {
                row.k = k;
                row.j = j, row.i = i;
                colidx = 0;
                col[colidx].i = i;
                col[colidx].j = j;
                col[colidx].k = k;
                val[colidx++] = pLsdft->coeffs_x[0] + pLsdft->coeffs_y[0] + pLsdft->coeffs_z[0];
                for (l = 1; l <= o; l++) {
                    col[colidx].i = i;
                    col[colidx].j = j;
                    col[colidx].k = k - l;
                    val[colidx++] = pLsdft->coeffs_z[l];
                    col[colidx].i = i;
                    col[colidx].j = j;
                    col[colidx].k = k + l;
                    val[colidx++] = pLsdft->coeffs_z[l];
                    col[colidx].i = i;
                    col[colidx].j = j - l;
                    col[colidx].k = k;
                    val[colidx++] = pLsdft->coeffs_y[l];
                    col[colidx].i = i;
                    col[colidx].j = j + l;
                    col[colidx].k = k;
                    val[colidx++] = pLsdft->coeffs_y[l];
                    col[colidx].i = i - l;
                    col[colidx].j = j;
                    col[colidx].k = k;
                    val[colidx++] = pLsdft->coeffs_x[l];
                    col[colidx].i = i + l;
                    col[colidx].j = j;
                    col[colidx].k = k;
                    val[colidx++] = pLsdft->coeffs_x[l];
                }
                MatSetValuesStencil(A, 1, &row, 6 * o + 1, col, val, ADD_VALUES);
                // MatSetValuesStencil(pLsdft->HamiltonianOpr,1,&row,6*o+1,col,val,ADD_VALUES);
            }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    //MatAssemblyBegin(pLsdft->HamiltonianOpr, MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(pLsdft->HamiltonianOpr, MAT_FINAL_ASSEMBLY);

#if _DEBUG
    MatIsSymmetric(A, 1.e-15, &flg);
    assert(flg);
#endif

#ifdef _DEBUG
    MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Is symmetric: %d\n", flg);
    PetscPrintf(PETSC_COMM_WORLD, "Istart: %d\n", Istart);
    PetscPrintf(PETSC_COMM_WORLD, "Iend: %d\n", Iend);
#endif


    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//                     Create_Laplacian: Initializes the -(1/2)*Laplacian operator           //
///////////////////////////////////////////////////////////////////////////////////////////////
PetscErrorCode Create_LaplacianLocal(LSDFT_OBJ *pLsdft)
{
    PetscInt i, j, k, l, colidx, gxdim, gydim, gzdim, xcor, ycor, zcor, lxdim, lydim, lzdim;
    MatStencil row;
    MatStencil *col;
    PetscScalar *val;
    //Mat A = pLsdft->LaplacianOprloc;
    PetscInt o = pLsdft->order;
#if _DEBUG
    PetscTruth flg;
#endif

    // we store the laplacian directly into the hamiltonian operator

    DMDAGetInfo(pLsdft->daloc, 0, &gxdim, &gydim, &gzdim, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    DMDAGetCorners(pLsdft->daloc, &xcor, &ycor, &zcor, &lxdim, &lydim, &lzdim);

    MatScale(pLsdft->LaplacianOprloc, 0.0);
    //MatScale(pLsdft->HamiltonianOpr,0.0);

    PetscMalloc(sizeof(MatStencil) * (o * 6 + 1), &col);
    PetscMalloc(sizeof(PetscScalar) * (o * 6 + 1), &val);

    for (k = zcor; k < zcor + lzdim; k++)
        for (j = ycor; j < ycor + lydim; j++)
            for (i = xcor; i < xcor + lxdim; i++) {
                row.k = k;
                row.j = j;
                row.i = i;
                colidx = 0;
                col[colidx].i = i;
                col[colidx].j = j;
                col[colidx].k = k;
                val[colidx++] = pLsdft->coeffs_x[0] + pLsdft->coeffs_y[0] + pLsdft->coeffs_z[0];
                for (l = 1; l <= o; l++) {
                    col[colidx].i = i;
                    col[colidx].j = j;
                    col[colidx].k = k - l;
                    val[colidx++] = pLsdft->coeffs_z[l];
                    col[colidx].i = i;
                    col[colidx].j = j;
                    col[colidx].k = k + l;
                    val[colidx++] = pLsdft->coeffs_z[l];
                    col[colidx].i = i;
                    col[colidx].j = j - l;
                    col[colidx].k = k;
                    val[colidx++] = pLsdft->coeffs_y[l];
                    col[colidx].i = i;
                    col[colidx].j = j + l;
                    col[colidx].k = k;
                    val[colidx++] = pLsdft->coeffs_y[l];
                    col[colidx].i = i - l;
                    col[colidx].j = j;
                    col[colidx].k = k;
                    val[colidx++] = pLsdft->coeffs_x[l];
                    col[colidx].i = i + l;
                    col[colidx].j = j;
                    col[colidx].k = k;
                    val[colidx++] = pLsdft->coeffs_x[l];
                }
                // MatSetValuesStencil(A,1,&row,6*o+1,col,val,ADD_VALUES);
                MatSetValuesStencil(pLsdft->LaplacianOprloc, 1, &row, 6 * o + 1, col, val, ADD_VALUES);
            }
    //MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(pLsdft->LaplacianOprloc, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(pLsdft->LaplacianOprloc, MAT_FINAL_ASSEMBLY);

    /// comment out later
    // int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    //if(rank==6)
    //  MatView(pLsdft->LaplacianOprloc,PETSC_VIEWER_STDOUT_SELF);

     // duplicate the local laplacian                                                                                                                                                  
    MatDuplicate(pLsdft->LaplacianOprloc,MAT_COPY_VALUES,&pLsdft->LapPlusVeffOprloc);
    
    return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////
//                         Create_Gradient: Initializes the Gradient operators              //
///////////////////////////////////////////////////////////////////////////////////////////////
PetscErrorCode Create_Gradient(LSDFT_OBJ *pLsdft)
{

    PetscInt i, j, k, l, colidx, gxdim, gydim, gzdim, xcor, ycor, zcor, lxdim, lydim, lzdim;
    MatStencil row;
    MatStencil *col_x;
    MatStencil *col_y;
    MatStencil *col_z;
    PetscScalar *val_x, *val_y, *val_z;
    PetscInt o = pLsdft->order;

    DMDAGetInfo(pLsdft->daloc, 0, &gxdim, &gydim, &gzdim, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    DMDAGetCorners(pLsdft->daloc, &xcor, &ycor, &zcor, &lxdim, &lydim, &lzdim);

    MatScale(pLsdft->gradient_x, 0.0);
    MatScale(pLsdft->gradient_y, 0.0);
    MatScale(pLsdft->gradient_z, 0.0);

    PetscMalloc(sizeof(MatStencil) * (o * 2 + 1), &col_x);
    PetscMalloc(sizeof(MatStencil) * (o * 2 + 1), &col_y);
    PetscMalloc(sizeof(MatStencil) * (o * 2 + 1), &col_z);
    PetscMalloc(sizeof(PetscScalar) * (o * 2 + 1), &val_x);
    PetscMalloc(sizeof(PetscScalar) * (o * 2 + 1), &val_y);
    PetscMalloc(sizeof(PetscScalar) * (o * 2 + 1), &val_z);

    for (k = zcor; k < zcor + lzdim; k++)
        for (j = ycor; j < ycor + lydim; j++)
            for (i = xcor; i < xcor + lxdim; i++) {
                row.k = k;
                row.j = j, row.i = i;
                colidx = 0;

                col_x[colidx].i = i;
                col_x[colidx].j = j;
                col_x[colidx].k = k;
                col_y[colidx].i = i;
                col_y[colidx].j = j;
                col_y[colidx].k = k;
                col_z[colidx].i = i;
                col_z[colidx].j = j;
                col_z[colidx].k = k;
                val_x[colidx] = pLsdft->coeffs_grad_x[0];
                val_y[colidx] = pLsdft->coeffs_grad_y[0];
                val_z[colidx] = pLsdft->coeffs_grad_z[0];
                colidx++;

                for (l = 1; l <= o; l++) {
                    col_x[colidx].i = i - l;
                    col_x[colidx].j = j;
                    col_x[colidx].k = k;
                    col_y[colidx].i = i;
                    col_y[colidx].j = j - l;
                    col_y[colidx].k = k;
                    col_z[colidx].i = i;
                    col_z[colidx].j = j;
                    col_z[colidx].k = k - l;
                    val_x[colidx] = -pLsdft->coeffs_grad_x[l];
                    val_y[colidx] = -pLsdft->coeffs_grad_y[l];
                    val_z[colidx] = -pLsdft->coeffs_grad_z[l];
                    colidx++;

                    col_x[colidx].i = i + l;
                    col_x[colidx].j = j;
                    col_x[colidx].k = k;
                    col_y[colidx].i = i;
                    col_y[colidx].j = j + l;
                    col_y[colidx].k = k;
                    col_z[colidx].i = i;
                    col_z[colidx].j = j;
                    col_z[colidx].k = k + l;
                    val_x[colidx] = pLsdft->coeffs_grad_x[l];
                    val_y[colidx] = pLsdft->coeffs_grad_y[l];
                    val_z[colidx] = pLsdft->coeffs_grad_z[l];
                    colidx++;
                }
                MatSetValuesStencil(pLsdft->gradient_x, 1, &row, 2 * o + 1, col_x, val_x, ADD_VALUES);
                MatSetValuesStencil(pLsdft->gradient_y, 1, &row, 2 * o + 1, col_y, val_y, ADD_VALUES);
                MatSetValuesStencil(pLsdft->gradient_z, 1, &row, 2 * o + 1, col_z, val_z, ADD_VALUES);
            }

    MatAssemblyBegin(pLsdft->gradient_x, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(pLsdft->gradient_x, MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(pLsdft->gradient_y, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(pLsdft->gradient_y, MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(pLsdft->gradient_z, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(pLsdft->gradient_z, MAT_FINAL_ASSEMBLY);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//                    Calculate_SplineFit: Calculates and stores the spline fit              //
///////////////////////////////////////////////////////////////////////////////////////////////
void Calculate_SplineFit(LSDFT_OBJ *pLsdft)
{
    int at, i, l;
    for (at = 0; at < pLsdft->Ntype; at++) {

        pLsdft->psd[at].UDeltaV = (PetscScalar **) malloc(sizeof(PetscScalar *) * ((pLsdft->psd[at].lmax) + 1));
        pLsdft->psd[at].SplineFitUDeltaV =
            (PetscScalar **) malloc(sizeof(PetscScalar *) * ((pLsdft->psd[at].lmax) + 1));
        pLsdft->psd[at].SplineFitU = (PetscScalar **) malloc(sizeof(PetscScalar *) * ((pLsdft->psd[at].lmax) + 1));
        pLsdft->psd[at].SplineFitVloc = (PetscScalar *) malloc(sizeof(PetscScalar) * pLsdft->psd[at].size);
        pLsdft->psd[at].SplineFitIsoAtomDen = (PetscScalar *) malloc(sizeof(PetscScalar) * pLsdft->psd[at].size);

        getYD_gen(pLsdft->psd[at].RadialGrid, pLsdft->psd[at].Vloc, pLsdft->psd[at].SplineFitVloc,
            pLsdft->psd[at].size);
        getYD_gen(pLsdft->psd[at].RadialGrid, pLsdft->psd[at].uu, pLsdft->psd[at].SplineFitIsoAtomDen,
            pLsdft->psd[at].size);

        for (l = 0; l <= pLsdft->psd[at].lmax; l++) {

            pLsdft->psd[at].UDeltaV[l] = (PetscScalar *) malloc(sizeof(PetscScalar) * pLsdft->psd[at].size);
            assert(pLsdft->psd[at].UDeltaV[l] != NULL);

            pLsdft->psd[at].SplineFitUDeltaV[l] = (PetscScalar *) malloc(sizeof(PetscScalar) * pLsdft->psd[at].size);
            assert(pLsdft->psd[at].SplineFitUDeltaV[l] != NULL);

            pLsdft->psd[at].SplineFitU[l] = (PetscScalar *) malloc(sizeof(PetscScalar) * pLsdft->psd[at].size);
            assert(pLsdft->psd[at].SplineFitU[l] != NULL);

            /*
             * compute nonlocal projectors from pseudopotential and pseudowavefunction
             */
            for (i = 0; i < pLsdft->psd[at].size; i++) {
                pLsdft->psd[at].UDeltaV[l][i] =
                    pLsdft->psd[at].U[l][i] * (pLsdft->psd[at].V[l][i] - pLsdft->psd[at].Vloc[i]);
            }
            /*
             * compute spline fit to the nonlocal projectors and pseudowavefunctions
             */
            getYD_gen(pLsdft->psd[at].RadialGrid, pLsdft->psd[at].UDeltaV[l], pLsdft->psd[at].SplineFitUDeltaV[l],
                pLsdft->psd[at].size);

            getYD_gen(pLsdft->psd[at].RadialGrid, pLsdft->psd[at].U[l], pLsdft->psd[at].SplineFitU[l],
                pLsdft->psd[at].size);

        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////
//              Calculate_PseudochargeCutoff: calculate pseudocharge density cutoff          //
///////////////////////////////////////////////////////////////////////////////////////////////
void Calculate_PseudochargeCutoff(LSDFT_OBJ *pLsdft)
{
    PetscLogDouble t2;

    PetscScalar ***pVpsArray, ***pBJArray, ***weights, rmax, MaxRadius = 15.0;
    PetscScalar tableR[MAX_TABLE_SIZE], tableVps[MAX_TABLE_SIZE], *YD = NULL;
    PetscScalar x, y, z, r, Dtemp, Bint, Rcut, error, val;
    PetscInt i, j, k, at, p, a, lloc;
    PetscInt Npts_x = ceil(MaxRadius / pLsdft->delta_x);
    PetscInt Npts_y = ceil(MaxRadius / pLsdft->delta_y);
    PetscInt Npts_z = ceil(MaxRadius / pLsdft->delta_z);
    PetscInt Ncube_x, Ncube_y, Ncube_z;
    PetscInt o = pLsdft->order;
    FILE *fOutFile;
    int count;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        if ((fOutFile = fopen((const char *) pLsdft->OutFilename, "a")) == NULL) {
            cout << "Couldn't open output file...\n" << pLsdft->OutFilename;
            cout << "Exiting...\n";
            exit(1);
        }
    }

    for (at = 0; at < pLsdft->Ntype; at++) {
        /*
         * initial cutoff radius
         */
        lloc = pLsdft->localPsd[at];

        Rcut = pLsdft->psd[at].rc[lloc];

        count = 0;

        PetscMalloc(sizeof(PetscScalar **) * (Npts_z + 2 * o), &pVpsArray);
        if (pVpsArray == NULL) {
            cout << "Memory allocation fail";
            exit(1);
        }
        for (i = 0; i < (Npts_z + 2 * o); i++) {
            PetscMalloc(sizeof(PetscScalar *) * (Npts_y + 2 * o), &pVpsArray[i]);
            if (pVpsArray[i] == NULL) {
                cout << "Memory allocation fail";
                exit(1);
            }
            for (j = 0; j < (Npts_y + 2 * o); j++) {
                PetscMalloc(sizeof(PetscScalar) * (Npts_x + 2 * o), &pVpsArray[i][j]);
                if (pVpsArray[i][j] == NULL) {
                    cout << "Memory allocation fail";
                    exit(1);
                }
            }
        }

        /*
         * we assume that the atom is situated on the origin and go over the grid
         */
        PetscMalloc(sizeof(PetscScalar **) * (Npts_z), &pBJArray);
        if (pBJArray == NULL) {
            cout << "Memory allocation fail";
            exit(1);
        }
        for (i = 0; i < (Npts_z); i++) {
            PetscMalloc(sizeof(PetscScalar *) * (Npts_y), &pBJArray[i]);
            if (pBJArray[i] == NULL) {
                cout << "Memory allocation fail";
                exit(1);
            }
            for (j = 0; j < (Npts_y); j++) {
                PetscMalloc(sizeof(PetscScalar) * (Npts_x), &pBJArray[i][j]);
                if (pBJArray[i][j] == NULL) {
                    cout << "Memory allocation fail";
                    exit(1);
                }
            }
        }

        /*
         * calculate pseudocharge density from the pseudopotential
         */
        for (k = 0; k < (Npts_z + 2 * o); k++)
            for (j = 0; j < (Npts_y + 2 * o); j++)
                for (i = 0; i < (Npts_x + 2 * o); i++) {
                    x = (i - o) * pLsdft->delta_x;
                    y = (j - o) * pLsdft->delta_y;
                    z = (k - o) * pLsdft->delta_z;
                    r = sqrt(x * x + y * y + z * z);
                    SplineEvaluate(pLsdft->psd[at].RadialGrid, pLsdft->psd[at].Vloc, pLsdft->psd[at].size, &r,
                        &pVpsArray[k][j][i], 1, pLsdft->psd[at].SplineFitVloc);
                }

        /*
         * calculate pseudocharge density from the pseudopotential
         */
        for (k = o; k < Npts_z + o; k++)
            for (j = o; j < Npts_y + o; j++)
                for (i = o; i < Npts_x + o; i++) {
                    pBJArray[k - o][j - o][i - o] = 0.0;
                    pBJArray[k - o][j - o][i - o] =
                        pVpsArray[k][j][i] * (pLsdft->stencil_coeffs_x[0] + pLsdft->stencil_coeffs_y[0] +
                        pLsdft->stencil_coeffs_z[0]);
                    for (a = 1; a <= o; a++) {
                        pBJArray[k - o][j - o][i - o] +=
                            (pVpsArray[k][j][i - a] + pVpsArray[k][j][i + a]) * pLsdft->stencil_coeffs_x[a] +
                            (pVpsArray[k][j - a][i] + pVpsArray[k][j + a][i]) * pLsdft->stencil_coeffs_y[a] +
                            (pVpsArray[k - a][j][i] + pVpsArray[k + a][j][i]) * pLsdft->stencil_coeffs_z[a];
                    }
                }
        /*
         * now go over different radii by incrementing radius by 1 bohr, unless relative
         * error in charge is below tolerence
         */
        do {
            Ncube_x = ceil(Rcut / pLsdft->delta_x);
            Ncube_y = ceil(Rcut / pLsdft->delta_y);
            Ncube_z = ceil(Rcut / pLsdft->delta_z);
            Bint = 0.0;
            for (k = 0; k < Ncube_z; k++)
                for (j = 0; j < Ncube_y; j++)
                    for (i = 0; i < Ncube_x; i++) {
                        val = pBJArray[k][j][i];
                        if (k == 0)
                            val = 0.5 * val;
                        if (k == Ncube_z - 1)
                            val = 0.5 * val;
                        if (j == 0)
                            val = 0.5 * val;
                        if (j == Ncube_y - 1)
                            val = 0.5 * val;
                        if (i == 0)
                            val = 0.5 * val;
                        if (i == Ncube_x - 1)
                            val = 0.5 * val;

                        Bint += val;
                    }
            Bint = Bint * 8.0 * pLsdft->delVol;
            error = fabs(Bint + pLsdft->noe[at]) / fabs(pLsdft->noe[at]);
            Rcut = Rcut + 1.0;
        } while (Rcut <= 15.0 && error > pLsdft->PseudochargeRadiusTOL);

        /*
         * store cutoff radius
         */
        pLsdft->CUTOFF_x[at] = Ncube_x * pLsdft->delta_x;
        pLsdft->CUTOFF_y[at] = Ncube_y * pLsdft->delta_y;
        pLsdft->CUTOFF_z[at] = Ncube_z * pLsdft->delta_z;

        if (rank == 0) {
            if ((fabs(pLsdft->delta_x - pLsdft->delta_y) <= 1e-15) && (fabs(pLsdft->delta_x - pLsdft->delta_z) <= 1e-15)
                && (fabs(pLsdft->delta_y - pLsdft->delta_z) <= 1e-15)) {
                fprintf(fOutFile, "Atom type %d: pseudocharge radius   : % E (Bohr)\n", at + 1, pLsdft->CUTOFF_x[at]);
            } else {
                fprintf(fOutFile, "Atom type %d: pseudocharge length x : % E (Bohr)\n", at + 1, pLsdft->CUTOFF_x[at]);
                fprintf(fOutFile, "Atom type %d: pseudocharge length y : % E (Bohr)\n", at + 1, pLsdft->CUTOFF_y[at]);
                fprintf(fOutFile, "Atom type %d: pseudocharge length z : % E (Bohr)\n", at + 1, pLsdft->CUTOFF_z[at]);
            }
        }
        for (i = 0; i < Npts_z + 2 * o; i++) {
            for (j = 0; j < Npts_y + 2 * o; j++) {
                PetscFree(pVpsArray[i][j]);
            }
            PetscFree(pVpsArray[i]);
        }
        PetscFree(pVpsArray);

        for (i = 0; i < Npts_z; i++) {
            for (j = 0; j < Npts_y; j++) {
                PetscFree(pBJArray[i][j]);
            }
            PetscFree(pBJArray[i]);
        }
        PetscFree(pBJArray);
        PetscFree(YD);

    }

    /*
     * end of initialization
     */
    PetscTime(&t2);
    pLsdft->timeInitialization = (t2 - pLsdft->timeInitialization_start);
    if (rank == 0) {
        fprintf(fOutFile, "Timing                                                 \n");
        fprintf(fOutFile, "Initialization                     : % E (sec)\n", pLsdft->timeInitialization);
        fclose(fOutFile);
    }
    return;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//         Calculate_NonlocalCoefficients: Calculate nonlocal coefficients                     //
///////////////////////////////////////////////////////////////////////////////////////////////
PetscScalar Calculate_NonlocalCoefficients(LSDFT_OBJ *pLsdft)
{

    PetscScalar rmax, rtemp, Rcut;
    PetscScalar dr = 1e-6;
    PetscErrorCode ierr;
    PetscScalar Dexact;
    int lmax, lloc, l, m, at;
    PetscScalar tableR[MAX_TABLE_SIZE], tableUlDeltaV[4][MAX_TABLE_SIZE], tableU[4][MAX_TABLE_SIZE];
    PetscScalar UlDelVl, Ul, Dtemp, delta;

    delta = max(pLsdft->delta_x, pLsdft->delta_y, pLsdft->delta_z);
    /*
     * first allocate memory for storing the coefficients
     */

    /*
     * loop over different types of atoms
     */
    for (at = 0; at < pLsdft->Ntype; at++) {
        lmax = pLsdft->psd[at].lmax;
        lloc = pLsdft->localPsd[at];

        pLsdft->psd[at].Gamma = (PetscScalar *) malloc(sizeof(PetscScalar) * (lmax + 1));

        for (l = 0; l <= lmax; l++) {
            pLsdft->psd[at].Gamma[l] = 0.0;
            rtemp = dr;
            Dexact = 0.0;
            if (l != lloc) {
                Rcut = pLsdft->psd[at].rc[l];
                while (rtemp <= (Rcut + delta)) {
                    SplineEvaluate(pLsdft->psd[at].RadialGrid, pLsdft->psd[at].UDeltaV[l], pLsdft->psd[at].size, &rtemp,
                        &UlDelVl, 1, pLsdft->psd[at].SplineFitUDeltaV[l]);
                    SplineEvaluate(pLsdft->psd[at].RadialGrid, pLsdft->psd[at].U[l], pLsdft->psd[at].size, &rtemp, &Ul,
                        1, pLsdft->psd[at].SplineFitU[l]);

                    Dexact += UlDelVl * Ul * rtemp * rtemp * dr;
                    rtemp = rtemp + dr;
                }
                Dexact = (Dexact - 0.5 * UlDelVl * Ul * (rtemp - dr) * (rtemp - dr) * dr) / (pLsdft->delVol);
                pLsdft->psd[at].Gamma[l] = 1.0 / Dexact;
            }
        }
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//          Calculate_nonlocalindex: evaluates the indices for storing inner product         //
//                                           in an array                                     //
///////////////////////////////////////////////////////////////////////////////////////////////
void Calculate_NonlocalIndex(LSDFT_OBJ *pLsdft)
{

    PetscInt at, l, m, n, poscnt, start, end, index, posIdx;
    pLsdft->T_indexlength = 0;
    pLsdft->T_length = 0;
    //TO BE CHANGED LATER:
    pLsdft->Nstates = 1;

    /*
     * find the length of the array storing the inner product with nonlocal pseudopotential operator
     */
    for (at = 0; at < pLsdft->Ntype; at++) {
        // start = (int) floor(pLsdft->startPos[at] / 3);
        // end = (int) floor(pLsdft->endPos[at] / 3);
	
	start = (int) floor(pLsdft->startPosSC[at] / 3);
        end = (int) floor(pLsdft->endPosSC[at] / 3);

        //if(pLsdft->psd[at].lmax!=0)
        {
            pLsdft->T_indexlength += pLsdft->noaSC[at];
            for (poscnt = start; poscnt <= end; poscnt++) {
                for (n = 0; n < pLsdft->Nstates; n++) {
                    for (l = 0; l <= pLsdft->psd[at].lmax; l++) {
                        if (l != pLsdft->localPsd[at]) {
                            for (m = -l; m <= l; m++)
                                pLsdft->T_length++;
                        }
                    }
                }
            }
        }
    }
    /*
     * store the indices
     */
    pLsdft->IP_length = (int) (pLsdft->T_length / pLsdft->Nstates);
    pLsdft->T_index = (PetscInt *) calloc(pLsdft->T_indexlength, sizeof(PetscInt));
    pLsdft->IP_index = (PetscInt *) calloc(pLsdft->T_indexlength, sizeof(PetscInt));


    index = 0;
    posIdx = 0;
    for (at = 0; at < pLsdft->Ntype; at++) {
        // start = (int) floor(pLsdft->startPos[at] / 3);
        // end = (int) floor(pLsdft->endPos[at] / 3);
      start = (int) floor(pLsdft->startPosSC[at] / 3);
        end = (int) floor(pLsdft->endPosSC[at] / 3);

        // if(pLsdft->psd[at].lmax!=0)
        {
            for (poscnt = start; poscnt <= end; poscnt++) {
                pLsdft->T_index[index++] = posIdx;
                for (n = 0; n < pLsdft->Nstates; n++) {
                    for (l = 0; l <= pLsdft->psd[at].lmax; l++) {
                        if (l != pLsdft->localPsd[at]) {
                            for (m = -l; m <= l; m++)
                                posIdx++;
                        }
                    }
                }
            }
        }
    }

    /*
     * store the indices for 1 inner product
     */
    index = 0;
    posIdx = 0;
    for (at = 0; at < pLsdft->Ntype; at++) {
        // start = (int) floor(pLsdft->startPos[at] / 3);
        // end = (int) floor(pLsdft->endPos[at] / 3);

	 start = (int) floor(pLsdft->startPosSC[at] / 3);
        end = (int) floor(pLsdft->endPosSC[at] / 3);

        // if(pLsdft->psd[at].lmax!=0)
        {
            for (poscnt = start; poscnt <= end; poscnt++) {
                pLsdft->IP_index[index++] = posIdx;
                for (l = 0; l <= pLsdft->psd[at].lmax; l++) {
                    if (l != pLsdft->localPsd[at]) {
                        for (m = -l; m <= l; m++)
                            posIdx++;
                    }
                }
            }
        }
    }

}

// ///////////////////////////////////////////////////////////////////////////////////////////////
// //                 CalculateNonlocalProjectors: Calculate nonlocal projectors                //
// ///////////////////////////////////////////////////////////////////////////////////////////////
// PetscScalar CalculateNonlocalProjectors(LSDFT_OBJ *pLsdft)
// {
//     PetscInt at, poscnt, l, m, k, j, i, nz, nx, ny;
//     PetscScalar x, y, z, r, zz, yy, SpHarmonic, pUlDeltaVl;
//     int rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     printf("ok1\n");
//     pLsdft->Projector_nonlocal = (PROJECTOR_OBJ *) malloc(sizeof(PROJECTOR_OBJ) * pLsdft->Ntype);
//     assert(pLsdft->Projector_nonlocal != NULL);
//     printf("ok2\n");
//     for (at = 0; at < pLsdft->Ntype; at++) {
//         pLsdft->Projector_nonlocal[at].Chi =
//             (double ******) malloc(sizeof(double *****) * (pLsdft->AtomOverlap_nonlocal[at].Natoms));
//         assert(pLsdft->Projector_nonlocal[at].Chi != NULL);
// 	printf("ok3\n");
//         for (poscnt = 0; poscnt < pLsdft->AtomOverlap_nonlocal[at].Natoms; poscnt++) {
//             nz = pLsdft->AtomOverlap_nonlocal[at].zend[poscnt] - pLsdft->AtomOverlap_nonlocal[at].zstart[poscnt] + 1;
//             ny = pLsdft->AtomOverlap_nonlocal[at].yend[poscnt] - pLsdft->AtomOverlap_nonlocal[at].ystart[poscnt] + 1;
//             nx = pLsdft->AtomOverlap_nonlocal[at].xend[poscnt] - pLsdft->AtomOverlap_nonlocal[at].xstart[poscnt] + 1;

//             pLsdft->Projector_nonlocal[at].Chi[poscnt] =
//                 (double *****) malloc(sizeof(double ****) * (pLsdft->psd[at].lmax + 1));
//             assert(pLsdft->Projector_nonlocal[at].Chi[poscnt] != NULL);
//             for (l = 0; l <= pLsdft->psd[at].lmax; l++) {
//                 if (l != pLsdft->localPsd[at]) {
//                     pLsdft->Projector_nonlocal[at].Chi[poscnt][l] =
//                         (double ****) malloc(sizeof(double ***) * (2 * l + 1));
//                     assert(pLsdft->Projector_nonlocal[at].Chi[poscnt][l] != NULL);
//                     for (m = -l; m <= l; m++) {
//                         pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m + l] =
//                             (double ***) malloc(sizeof(double **) * nz);
//                         assert(pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m + l] != NULL);
//                         for (k = 0; k < nz; k++) {
//                             z = pLsdft->delta_z * (k + pLsdft->AtomOverlap_nonlocal[at].zstart[poscnt]) -
//                                 pLsdft->range_z - pLsdft->AtomOverlap_nonlocal[at].Z0[poscnt];
//                             zz = z * z;
//                             pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m + l][k] =
//                                 (double **) malloc(sizeof(double *) * ny);
//                             assert(pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m + l][k] != NULL);

//                             for (j = 0; j < ny; j++) {
//                                 y = pLsdft->delta_y * (j + pLsdft->AtomOverlap_nonlocal[at].ystart[poscnt]) -
//                                     pLsdft->range_y - pLsdft->AtomOverlap_nonlocal[at].Y0[poscnt];
//                                 yy = y * y;
//                                 pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m + l][k][j] =
//                                     (double *) malloc(sizeof(double) * nx);
//                                 assert(pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m + l][k][j] != NULL);
//                                 for (i = 0; i < nx; i++) {
//                                     x = pLsdft->delta_x * (i + pLsdft->AtomOverlap_nonlocal[at].xstart[poscnt]) -
//                                         pLsdft->range_x - pLsdft->AtomOverlap_nonlocal[at].X0[poscnt];

// 				    //printf("SS: i=%d, X=%lf, X0=%lf,x=%lf \n",i,pLsdft->delta_x*(i + pLsdft->AtomOverlap_nonlocal[at].xstart[poscnt])-pLsdft->range_x,pLsdft->AtomOverlap_nonlocal[at].X0[poscnt],x);
//                                     r = sqrt(x * x + yy + zz);
//                                     SplineEvaluate(pLsdft->psd[at].RadialGrid, pLsdft->psd[at].UDeltaV[l],
//                                         pLsdft->psd[at].size, &r, &pUlDeltaVl, 1, pLsdft->psd[at].SplineFitUDeltaV[l]);
//                                     SpHarmonic = RealSphericalHarmonic(x, y, z, l, m);
//                                     pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m + l][k][j][i] =
//                                         pUlDeltaVl * SpHarmonic;
//                                 }
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     return 0;
// }

///////////////////////////////////////////////////////////////////////////////////////////////
//                 CalculateNonlocalProjectors: Calculate nonlocal projectors                //
///////////////////////////////////////////////////////////////////////////////////////////////
PetscScalar CalculateNonlocalProjectorsForces(LSDFT_OBJ *pLsdft)
{
    PetscInt at, poscnt, l, m, k, j, i, nz, nx, ny;
    PetscScalar x, y, z, r, zz, yy, SpHarmonic, pUlDeltaVl;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    pLsdft->Projector_nonlocalforces = (PROJECTOR_OBJ *) malloc(sizeof(PROJECTOR_OBJ) * pLsdft->Ntype);
    assert(pLsdft->Projector_nonlocalforces != NULL);
    for (at = 0; at < pLsdft->Ntype; at++) {
        pLsdft->Projector_nonlocalforces[at].Chi =
            (double ******) malloc(sizeof(double *****) * (pLsdft->AtomOverlap_nonlocalforces[at].Natoms));
        assert(pLsdft->Projector_nonlocalforces[at].Chi != NULL);
        for (poscnt = 0; poscnt < pLsdft->AtomOverlap_nonlocalforces[at].Natoms; poscnt++) {
            nz = pLsdft->AtomOverlap_nonlocalforces[at].zend[poscnt] - pLsdft->AtomOverlap_nonlocalforces[at].zstart[poscnt] + 1;
            ny = pLsdft->AtomOverlap_nonlocalforces[at].yend[poscnt] - pLsdft->AtomOverlap_nonlocalforces[at].ystart[poscnt] + 1;
            nx = pLsdft->AtomOverlap_nonlocalforces[at].xend[poscnt] - pLsdft->AtomOverlap_nonlocalforces[at].xstart[poscnt] + 1;

            pLsdft->Projector_nonlocalforces[at].Chi[poscnt] =
                (double *****) malloc(sizeof(double ****) * (pLsdft->psd[at].lmax + 1));
            assert(pLsdft->Projector_nonlocalforces[at].Chi[poscnt] != NULL);
            for (l = 0; l <= pLsdft->psd[at].lmax; l++) {
                if (l != pLsdft->localPsd[at]) {
                    pLsdft->Projector_nonlocalforces[at].Chi[poscnt][l] =
                        (double ****) malloc(sizeof(double ***) * (2 * l + 1));
                    assert(pLsdft->Projector_nonlocalforces[at].Chi[poscnt][l] != NULL);
                    for (m = -l; m <= l; m++) {
                        pLsdft->Projector_nonlocalforces[at].Chi[poscnt][l][m + l] =
                            (double ***) malloc(sizeof(double **) * nz);
                        assert(pLsdft->Projector_nonlocalforces[at].Chi[poscnt][l][m + l] != NULL);
                        for (k = 0; k < nz; k++) {
                            z = pLsdft->delta_z * (k + pLsdft->AtomOverlap_nonlocalforces[at].zstart[poscnt]) -
                                pLsdft->range_z - pLsdft->AtomOverlap_nonlocalforces[at].Z0[poscnt];
                            zz = z * z;
                            pLsdft->Projector_nonlocalforces[at].Chi[poscnt][l][m + l][k] =
                                (double **) malloc(sizeof(double *) * ny);
                            assert(pLsdft->Projector_nonlocalforces[at].Chi[poscnt][l][m + l][k] != NULL);

                            for (j = 0; j < ny; j++) {
                                y = pLsdft->delta_y * (j + pLsdft->AtomOverlap_nonlocalforces[at].ystart[poscnt]) -
                                    pLsdft->range_y - pLsdft->AtomOverlap_nonlocalforces[at].Y0[poscnt];
                                yy = y * y;
                                pLsdft->Projector_nonlocalforces[at].Chi[poscnt][l][m + l][k][j] =
                                    (double *) malloc(sizeof(double) * nx);
                                assert(pLsdft->Projector_nonlocalforces[at].Chi[poscnt][l][m + l][k][j] != NULL);
                                for (i = 0; i < nx; i++) {
                                    x = pLsdft->delta_x * (i + pLsdft->AtomOverlap_nonlocalforces[at].xstart[poscnt]) -
                                        pLsdft->range_x - pLsdft->AtomOverlap_nonlocalforces[at].X0[poscnt];

				    //printf("SS: i=%d, X=%lf, X0=%lf,x=%lf \n",i,pLsdft->delta_x*(i + pLsdft->AtomOverlap_nonlocal[at].xstart[poscnt])-pLsdft->range_x,pLsdft->AtomOverlap_nonlocal[at].X0[poscnt],x);
				    
				  

                                    r = sqrt(x * x + yy + zz);
                                    SplineEvaluate(pLsdft->psd[at].RadialGrid, pLsdft->psd[at].UDeltaV[l],
                                        pLsdft->psd[at].size, &r, &pUlDeltaVl, 1, pLsdft->psd[at].SplineFitUDeltaV[l]);
                                    SpHarmonic = RealSphericalHarmonic(x, y, z, l, m);
                                    pLsdft->Projector_nonlocalforces[at].Chi[poscnt][l][m + l][k][j][i] =
                                        pUlDeltaVl * SpHarmonic;
				    // printf("l=%d, m=%d, k=%d,j=%d,i=%d, x=%0.16lf,y=%0.16lf,z=%0.16lf, r=%0.16lf, SpHarmonic=%0.16lf \n",l,m,k,j,i,x,y,z,r,SpHarmonic);
				    //  printf("NL projector force: pLsdft->Projector_nonlocalforces[at].Chi[%d][%d][%d][%d][%d][%d]=%0.16lf \n",poscnt,l,m+l,k,j,i, pLsdft->Projector_nonlocalforces[at].Chi[poscnt][l][m + l][k][j][i]);

                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return 0;
}



///////////////////////////////////////////////////////////////////////////////////////////////
//                               Objects_Destroy: Destroys PETSc objects                     //
///////////////////////////////////////////////////////////////////////////////////////////////
void Objects_Destroy(LSDFT_OBJ *pLsdft)
{

    int at, i;
    int MAX_ITS_ANDERSON = pLsdft->MixingHistory;
    PetscInt xcor, ycor, zcor, lxdim, lydim, lzdim;


    /*
     * destroy vectors
     */
    VecDestroy(&pLsdft->elecDensRho);
    VecDestroy(&pLsdft->SuperposAtRho);
    VecDestroy(&pLsdft->chrgDensB);
    VecDestroy(&pLsdft->chrgDensB_TM);
    VecDestroy(&pLsdft->potentialPhi);
    VecDestroy(&pLsdft->twopiRhoPB);
    VecDestroy(&pLsdft->bjVj);
    VecDestroy(&pLsdft->bjVj_TM);
    VecDestroy(&pLsdft->Vxc);
    VecDestroy(&pLsdft->forces);
    VecDestroy(&pLsdft->Atompos);
    VecDestroy(&pLsdft->mvAtmConstraint);
    VecDestroy(&pLsdft->Veff);
    VecDestroy(&pLsdft->Phi_c);
    VecDestroy(&pLsdft->twopiBTMmBPS);
    VecDestroy(&pLsdft->PoissonRHSAdd);
    VecDestroy(&pLsdft->tempVec);
    VecDestroy(&pLsdft->xkprev);
    VecDestroy(&pLsdft->xk);
    VecDestroy(&pLsdft->fkprev);

    for (i = 0; i < MAX_ITS_ANDERSON; i++) {
        VecDestroy(&pLsdft->Xk[i]);
        VecDestroy(&pLsdft->Fk[i]);
        VecDestroy(&pLsdft->XpbF[i]);
    }

    /*
     * destroy matrices
     */
    MatDestroy(&pLsdft->laplaceOpr);
    MatDestroy(&pLsdft->gradient_x);
    MatDestroy(&pLsdft->gradient_y);
    MatDestroy(&pLsdft->gradient_z);

    // if(pLsdft->Nkpts==1)
    //   {
    MatDestroy(&pLsdft->XOrb);
    MatDestroy(&pLsdft->YOrb);
    MatDestroy(&pLsdft->YOrbNew);
    MatDestroy(&pLsdft->Hsub);
    MatDestroy(&pLsdft->Msub);
    // }else{
    // for(i=0;i<pLsdft->Nkpts_symGroup;i++)
    //   {
    //  MatDestroy(&pLsdft->XOrb_real[i]);
    //  MatDestroy(&pLsdft->YOrb_real[i]);
    //  MatDestroy(&pLsdft->YOrbNew_real);
    //  MatDestroy(&pLsdft->XOrb_imag[i]);
    //  MatDestroy(&pLsdft->YOrb_imag[i]);
    //  MatDestroy(&pLsdft->YOrbNew_imag);
    //  MatDestroy(&pLsdft->ZOrb1);
    //  MatDestroy(&pLsdft->Hsub);
    //  MatDestroy(&pLsdft->Msub);
    //  MatDestroy(&pLsdft->Hsub_imag);
    //  MatDestroy(&pLsdft->Msub_imag);
    //   }
    // free(pLsdft->k1);
    // free(pLsdft->k2);
    // free(pLsdft->k3);
    // free(pLsdft->lambdakpt);
    // }

    // destroy local vectors
    VecDestroy(&pLsdft->Veff_loc);
    VecDestroy(&pLsdft->Vk);
    VecDestroy(&pLsdft->Vkm1);
    VecDestroy(&pLsdft->Vkp1);

    /*
     * destroy ksp
     */
    KSPDestroy(&pLsdft->ksp);
    /*
     * vecscatter destroy
     */
    VecScatterDestroy(&pLsdft->VeffScatterCtx);

/*
   * free  memory for lssgq weights and nodes
   */
    DMDAGetCorners(pLsdft->da, &xcor, &ycor, &zcor, &lxdim, &lydim, &lzdim);


    for (i = 0; i < (lxdim * lydim * lzdim); i++) {
        PetscFree(pLsdft->lssgq_weights[i]);
        PetscFree(pLsdft->lssgq_lambda[i]);
    }
    PetscFree(pLsdft->lssgq_weights);
    PetscFree(pLsdft->lssgq_lambda);


    /*
     * free memory occupied by pseudopotentials
     */
    for (at = 0; at < pLsdft->Ntype; at++) {
        free(pLsdft->psd[at].Vloc);
        for (i = 0; i <= pLsdft->psd[at].lmax; i++) {
            free(pLsdft->psd[at].V[i]);
            free(pLsdft->psd[at].U[i]);
        }
        free(pLsdft->psd[at].U);
        free(pLsdft->psd[at].V);
        free(pLsdft->psd[at].uu);
        free(pLsdft->psd[at].RadialGrid);
    }
    free(pLsdft->psd);

    /*
     * free memory occupied by arrays
     //  */
    // PetscFree(pLsdft->nnzDArray);
    // PetscFree(pLsdft->nnzODArray);

    PetscFree(pLsdft->lambda);
    PetscFree(pLsdft->CUTOFF_x);
    PetscFree(pLsdft->CUTOFF_y);
    PetscFree(pLsdft->CUTOFF_z);
    PetscFree(pLsdft->startPos);
    PetscFree(pLsdft->endPos);
    PetscFree(pLsdft->localPsd);

    /*
     * free communicators
     */
    // MPI_Comm_free(&pLsdft->kpoint_group_comm);
    // MPI_Comm_free(&pLsdft->kpoint_intergroup_comm);


    return;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//                        Set_VecZero: Initializes vectors with zero entries                 //
///////////////////////////////////////////////////////////////////////////////////////////////
void Set_VecZero(LSDFT_OBJ *pLsdft)
{

    PetscInt at, J, l, m, nx, ny, nz, i, j, k;

    VecZeroEntries(pLsdft->chrgDensB);
    VecZeroEntries(pLsdft->chrgDensB_TM);
    VecZeroEntries(pLsdft->potentialPhi);
    VecZeroEntries(pLsdft->Phi_c);
    VecZeroEntries(pLsdft->twopiRhoPB);
    VecZeroEntries(pLsdft->twopiBTMmBPS);
    VecZeroEntries(pLsdft->bjVj);
    VecZeroEntries(pLsdft->bjVj_TM);
    VecZeroEntries(pLsdft->Veff);
    VecZeroEntries(pLsdft->Vxc);
    VecZeroEntries(pLsdft->PoissonRHSAdd);
    VecZeroEntries(pLsdft->tempVec);

    for (at = 0; at < pLsdft->Ntype; at++) {
        for (J = 0; J < pLsdft->AtomOverlap_nonlocal[at].Natoms; J++) {
            nz = pLsdft->AtomOverlap_nonlocal[at].zend[J] - pLsdft->AtomOverlap_nonlocal[at].zstart[J] + 1;
            ny = pLsdft->AtomOverlap_nonlocal[at].yend[J] - pLsdft->AtomOverlap_nonlocal[at].ystart[J] + 1;
            nx = pLsdft->AtomOverlap_nonlocal[at].xend[J] - pLsdft->AtomOverlap_nonlocal[at].xstart[J] + 1;
            for (l = 0; l <= pLsdft->psd[at].lmax; l++) {
                if (l != pLsdft->localPsd[at]) {
                    for (m = -l; m <= l; m++) {
                        for (k = 0; k < nz; k++) {
                            for (j = 0; j < ny; j++) {
                                free(pLsdft->Projector_nonlocal[at].Chi[J][l][m + l][k][j]);
                            }
                            free(pLsdft->Projector_nonlocal[at].Chi[J][l][m + l][k]);
                        }
                        free(pLsdft->Projector_nonlocal[at].Chi[J][l][m + l]);
                    }
                    free(pLsdft->Projector_nonlocal[at].Chi[J][l]);
                }
            }
            free(pLsdft->Projector_nonlocal[at].Chi[J]);
        }
        free(pLsdft->Projector_nonlocal[at].Chi);


	for (J = 0; J < pLsdft->AtomOverlap_nonlocalforces[at].Natoms; J++) {
            nz = pLsdft->AtomOverlap_nonlocalforces[at].zend[J] - pLsdft->AtomOverlap_nonlocalforces[at].zstart[J] + 1;
            ny = pLsdft->AtomOverlap_nonlocalforces[at].yend[J] - pLsdft->AtomOverlap_nonlocalforces[at].ystart[J] + 1;
            nx = pLsdft->AtomOverlap_nonlocalforces[at].xend[J] - pLsdft->AtomOverlap_nonlocalforces[at].xstart[J] + 1;
            for (l = 0; l <= pLsdft->psd[at].lmax; l++) {
                if (l != pLsdft->localPsd[at]) {
                    for (m = -l; m <= l; m++) {
                        for (k = 0; k < nz; k++) {
                            for (j = 0; j < ny; j++) {
                                free(pLsdft->Projector_nonlocalforces[at].Chi[J][l][m + l][k][j]);
                            }
                            free(pLsdft->Projector_nonlocalforces[at].Chi[J][l][m + l][k]);
                        }
                        free(pLsdft->Projector_nonlocalforces[at].Chi[J][l][m + l]);
                    }
                    free(pLsdft->Projector_nonlocalforces[at].Chi[J][l]);
                }
            }
            free(pLsdft->Projector_nonlocalforces[at].Chi[J]);
        }
        free(pLsdft->Projector_nonlocalforces[at].Chi);



        free(pLsdft->AtomOverlap_local[at].X0);
        free(pLsdft->AtomOverlap_local[at].Y0);
        free(pLsdft->AtomOverlap_local[at].Z0);
        free(pLsdft->AtomOverlap_local[at].xstart);
        free(pLsdft->AtomOverlap_local[at].ystart);
        free(pLsdft->AtomOverlap_local[at].zstart);
        free(pLsdft->AtomOverlap_local[at].xend);
        free(pLsdft->AtomOverlap_local[at].yend);
        free(pLsdft->AtomOverlap_local[at].zend);
        free(pLsdft->AtomOverlap_local[at].xindex);
        free(pLsdft->AtomOverlap_local[at].yindex);
        free(pLsdft->AtomOverlap_local[at].zindex);

        free(pLsdft->AtomOverlap_nonlocal[at].X0);
        free(pLsdft->AtomOverlap_nonlocal[at].Y0);
        free(pLsdft->AtomOverlap_nonlocal[at].Z0);
        free(pLsdft->AtomOverlap_nonlocal[at].xstart);
        free(pLsdft->AtomOverlap_nonlocal[at].ystart);
        free(pLsdft->AtomOverlap_nonlocal[at].zstart);
        free(pLsdft->AtomOverlap_nonlocal[at].xend);
        free(pLsdft->AtomOverlap_nonlocal[at].yend);
        free(pLsdft->AtomOverlap_nonlocal[at].zend);
        free(pLsdft->AtomOverlap_nonlocal[at].xindex);
        free(pLsdft->AtomOverlap_nonlocal[at].yindex);
        free(pLsdft->AtomOverlap_nonlocal[at].zindex);

	free(pLsdft->AtomOverlap_localforces[at].X0);
        free(pLsdft->AtomOverlap_localforces[at].Y0);
        free(pLsdft->AtomOverlap_localforces[at].Z0);
        free(pLsdft->AtomOverlap_localforces[at].xstart);
        free(pLsdft->AtomOverlap_localforces[at].ystart);
        free(pLsdft->AtomOverlap_localforces[at].zstart);
        free(pLsdft->AtomOverlap_localforces[at].xend);
        free(pLsdft->AtomOverlap_localforces[at].yend);
        free(pLsdft->AtomOverlap_localforces[at].zend);
        free(pLsdft->AtomOverlap_localforces[at].xindex);
        free(pLsdft->AtomOverlap_localforces[at].yindex);
        free(pLsdft->AtomOverlap_localforces[at].zindex);

        free(pLsdft->AtomOverlap_nonlocalforces[at].X0);
        free(pLsdft->AtomOverlap_nonlocalforces[at].Y0);
        free(pLsdft->AtomOverlap_nonlocalforces[at].Z0);
        free(pLsdft->AtomOverlap_nonlocalforces[at].xstart);
        free(pLsdft->AtomOverlap_nonlocalforces[at].ystart);
        free(pLsdft->AtomOverlap_nonlocalforces[at].zstart);
        free(pLsdft->AtomOverlap_nonlocalforces[at].xend);
        free(pLsdft->AtomOverlap_nonlocalforces[at].yend);
        free(pLsdft->AtomOverlap_nonlocalforces[at].zend);
        free(pLsdft->AtomOverlap_nonlocalforces[at].xindex);
        free(pLsdft->AtomOverlap_nonlocalforces[at].yindex);
        free(pLsdft->AtomOverlap_nonlocalforces[at].zindex);


    }
    free(pLsdft->AtomOverlap_local);
    free(pLsdft->AtomOverlap_nonlocal);
    free(pLsdft->Projector_nonlocal);

     free(pLsdft->AtomOverlap_localforces);
    free(pLsdft->AtomOverlap_nonlocalforces);
    free(pLsdft->Projector_nonlocalforces);


    // if(pLsdft->Nkpts>1)
    //   {
    //     free(pLsdft->Projector_nonlocal_imag);
    //   }

    return;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//    GetInfluencingAtoms: Gets the list of atoms that influence the domain of the processor  //
//                               additionaly stores details of the overlap                    //
////////////////////////////////////////////////////////////////////////////////////////////////
void GetInfluencingAtoms_Local(LSDFT_OBJ *pLsdft)
{

    int count_local = 0;
    //int count_nonlocal=0;
    int NatomDomain = 0;
    int start, end, offset_local_x, offset_nonlocal_x, offset_local_y, offset_nonlocal_y, offset_local_z,
        offset_nonlocal_z, index = 0, posindex = 0, at, poscnt;
    PetscInt xcor, ycor, zcor, lxdim, lydim, lzdim;
    PetscInt xs, ys, zs, xl, yl, zl, xs_nl, ys_nl, zs_nl, xl_nl, yl_nl, zl_nl, xi, yj, zk, xstart = -1, ystart =
        -1, zstart = -1, xend = -1, yend = -1, zend = -1, overlap_local = 0, overlap_nonlocal = 0;
    PetscInt xstart_nl = -1, ystart_nl = -1, zstart_nl = -1, xend_nl = -1, yend_nl = -1, zend_nl = -1;
    PetscScalar cutoffr_nonlocal, cutoffr_local_x, cutoffr_local_y, cutoffr_local_z, max1, max2;
    int Imax_x = 0, Imin_x = 0, Imax_y = 0, Imin_y = 0, Imax_z = 0, Imin_z = 0;
    int PP, QQ, RR;

    PetscScalar x0, y0, z0, X0, Y0, Z0;
    PetscScalar *pAtompos;

    PetscScalar R_x = pLsdft->range_x;
    PetscScalar R_y = pLsdft->range_y;
    PetscScalar R_z = pLsdft->range_z;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    DMDAGetCorners(pLsdft->da, &xcor, &ycor, &zcor, &lxdim, &lydim, &lzdim);
    //xcor=0;ycor=0;zcor=0;lxdim=pLsdft->numPoints_x;lydim=pLsdft->numPoints_y;lzdim=pLsdft->numPoints_z;

    /*
     * allocate memory for Atom overlap object in the structure
     */
    pLsdft->AtomOverlap_local = (OVERLAP_OBJ *) malloc(sizeof(OVERLAP_OBJ) * pLsdft->Ntype);
    
    assert(pLsdft->AtomOverlap_local != NULL);
   

    /*
     * first go over the list of all atoms and determine the number of atoms which have overlap
     */

   
     VecGetArray(pLsdft->AtomposSC, &pAtompos);

    /*
     * loop over different types of atoms
     */
    for (at = 0; at < pLsdft->Ntype; at++) {
        cutoffr_local_x = pLsdft->CUTOFF_x[at];
        cutoffr_local_y = pLsdft->CUTOFF_y[at];
        cutoffr_local_z = pLsdft->CUTOFF_z[at];

        offset_local_x = (int) ceil(cutoffr_local_x / pLsdft->delta_x+0.5);
        offset_local_y = (int) ceil(cutoffr_local_y / pLsdft->delta_y+0.5);
        offset_local_z = (int) ceil(cutoffr_local_z / pLsdft->delta_z+0.5);

        // max1=  pLsdft->psd[at].rc[0] > pLsdft->psd[at].rc[1] ? pLsdft->psd[at].rc[0]:pLsdft->psd[at].rc[1];
        // max2=  pLsdft->psd[at].rc[2] > pLsdft->psd[at].rc[3] ? pLsdft->psd[at].rc[2]:pLsdft->psd[at].rc[3];
        // cutoffr_nonlocal =  max1>max2 ? max1:max2;
        // offset_nonlocal_x = ceil(cutoffr_nonlocal/pLsdft->delta_x);
        // offset_nonlocal_y = ceil(cutoffr_nonlocal/pLsdft->delta_y);
        // offset_nonlocal_z = ceil(cutoffr_nonlocal/pLsdft->delta_z);

	//    start = (int) floor(pLsdft->startPos[at] / 3);
	// end = (int) floor(pLsdft->endPos[at] / 3);


	start = (int) floor(pLsdft->startPosSC[at] / 3);
        end = (int) floor(pLsdft->endPosSC[at] / 3);

	printf("at=%d, start=%d, end=%d \n",at,start,end);

        count_local = 0;
        // count_nonlocal=0;

        Imax_x = 0;
        Imin_x = 0;
        Imax_y = 0;
        Imin_y = 0;
        Imax_z = 0;
        Imin_z = 0;

        if (pLsdft->BC == 1)    // 3D periodic
        {
            Imax_x = ceil(cutoffr_local_x / R_x);
            Imin_x = -ceil(cutoffr_local_x / R_x);

            Imax_y = ceil(cutoffr_local_y / R_y);
            Imin_y = -ceil(cutoffr_local_y / R_y);

            Imax_z = ceil(cutoffr_local_z / R_z);
            Imin_z = -ceil(cutoffr_local_z / R_z);
        } else if (pLsdft->BC == 3) // 2D periodic
        {
            Imax_x = ceil(cutoffr_local_x / R_x);
            Imin_x = -ceil(cutoffr_local_x / R_x);

            Imax_y = ceil(cutoffr_local_y / R_y);
            Imin_y = -ceil(cutoffr_local_y / R_y);
        } else if (pLsdft->BC == 4) // 1D periodic
        {
            Imax_x = ceil(cutoffr_local_x / R_x);
            Imin_x = -ceil(cutoffr_local_x / R_x);
        }


        /*
         * loop over every atom of a given type
         */
        for (poscnt = start; poscnt <= end; poscnt++) {
            X0 = pAtompos[index++];
            Y0 = pAtompos[index++];
            Z0 = pAtompos[index++];

            for (PP = Imin_x; PP <= Imax_x; PP++)
                for (QQ = Imin_y; QQ <= Imax_y; QQ++)
                    for (RR = Imin_z; RR <= Imax_z; RR++) {

                        /*
                         * periodic map of atomic position
                         */
                        x0 = X0 + PP * 2.0 * R_x;
                        y0 = Y0 + QQ * 2.0 * R_y;
                        z0 = Z0 + RR * 2.0 * R_z;

                        xi = roundf((x0 + R_x) / pLsdft->delta_x);
                        yj = roundf((y0 + R_y) / pLsdft->delta_y);
                        zk = roundf((z0 + R_z) / pLsdft->delta_z);

                        xs = xi - offset_local_x;
                        xl = xi + offset_local_x;
                        ys = yj - offset_local_y;
                        yl = yj + offset_local_y;
                        zs = zk - offset_local_z;
                        zl = zk + offset_local_z;

                        /*
                         * find if domain of influence of pseudocharge overlaps with the domain stored
                         * by processor
                         */
                        overlap_local =
                            Check_Overlap(xcor, ycor, zcor, lxdim, lydim, lzdim, xs, ys, zs, xl, yl, zl, &xstart,
                            &ystart, &zstart, &xend, &yend, &zend);
                        /*
                         * increment counter for pseudocharge
                         */
                        count_local += overlap_local;
			// if(xs==lxdim || ys==lydim || zs==lzdim)
			//   {
			//     printf("SSS xs=%d,ys=%d,zs=%d, xl=%d,yl=%d,zl=%d \n", xs, ys, zs, xl, yl, zl);
			//   }
                        /*
                         * check overlap for nonlocal
                         */
                        // if(overlap_local==1 && pLsdft->psd[at].lmax!=0)
                        //   {
                        //     xs_nl = xi-offset_nonlocal_x; xl_nl = xi+offset_nonlocal_x;
                        //     ys_nl = yj-offset_nonlocal_y; yl_nl = yj+offset_nonlocal_y;
                        //     zs_nl = zk-offset_nonlocal_z; zl_nl = zk+offset_nonlocal_z;

                        //     overlap_nonlocal = Check_Overlap(xcor,ycor,zcor,lxdim,lydim,lzdim,xs_nl,ys_nl,zs_nl,xl_nl,yl_nl,zl_nl,&xstart_nl,&ystart_nl,&zstart_nl,&xend_nl,&yend_nl,&zend_nl);
                        //     count_nonlocal +=overlap_nonlocal;
                        //   }
                    }
        }

        pLsdft->AtomOverlap_local[at].Natoms = count_local;
       

	printf("rank=%d, pLsdft->AtomOverlap_local[%d].Natoms=%d \n", rank,at, pLsdft->AtomOverlap_local[at].Natoms);


        /*
         * once the number of local and nonlocal influencing atoms are determined,
         * we allocate memory for the datastructure and store the positions, starting and ending indices of local regions
         */
        pLsdft->AtomOverlap_local[at].X0 = (double *) malloc(sizeof(double) * count_local);
        pLsdft->AtomOverlap_local[at].Y0 = (double *) malloc(sizeof(double) * count_local);
        pLsdft->AtomOverlap_local[at].Z0 = (double *) malloc(sizeof(double) * count_local);
        pLsdft->AtomOverlap_local[at].xstart = (int *) malloc(sizeof(int) * count_local);
        pLsdft->AtomOverlap_local[at].ystart = (int *) malloc(sizeof(int) * count_local);
        pLsdft->AtomOverlap_local[at].zstart = (int *) malloc(sizeof(int) * count_local);
        pLsdft->AtomOverlap_local[at].xend = (int *) malloc(sizeof(int) * count_local);
        pLsdft->AtomOverlap_local[at].yend = (int *) malloc(sizeof(int) * count_local);
        pLsdft->AtomOverlap_local[at].zend = (int *) malloc(sizeof(int) * count_local);

        pLsdft->AtomOverlap_local[at].xindex = (int *) malloc(sizeof(int) * count_local);
        pLsdft->AtomOverlap_local[at].yindex = (int *) malloc(sizeof(int) * count_local);
        pLsdft->AtomOverlap_local[at].zindex = (int *) malloc(sizeof(int) * count_local);

        count_local = 0;
  

        /*
         * loop over every atom of a given type
         */
        for (poscnt = start; poscnt <= end; poscnt++) {
            X0 = pAtompos[posindex++];
            Y0 = pAtompos[posindex++];
            Z0 = pAtompos[posindex++];

            for (PP = Imin_x; PP <= Imax_x; PP++)
                for (QQ = Imin_y; QQ <= Imax_y; QQ++)
                    for (RR = Imin_z; RR <= Imax_z; RR++) {
                        /*
                         * periodic map of atomic position
                         */
                        x0 = X0 + PP * 2.0 * R_x;
                        y0 = Y0 + QQ * 2.0 * R_y;
                        z0 = Z0 + RR * 2.0 * R_z;

                        xi = roundf((x0 + R_x) / pLsdft->delta_x);
                        yj = roundf((y0 + R_y) / pLsdft->delta_y);
                        zk = roundf((z0 + R_z) / pLsdft->delta_z);

                        xs = xi - offset_local_x;
                        xl = xi + offset_local_x;
                        ys = yj - offset_local_y;
                        yl = yj + offset_local_y;
                        zs = zk - offset_local_z;
                        zl = zk + offset_local_z;

                        overlap_local =
                            Check_Overlap(xcor, ycor, zcor, lxdim, lydim, lzdim, xs, ys, zs, xl, yl, zl, &xstart,
                            &ystart, &zstart, &xend, &yend, &zend);

                        if (overlap_local == 1) {
                            /*
                             * store data in structure
                             */
                            pLsdft->AtomOverlap_local[at].X0[count_local] = x0;
                            pLsdft->AtomOverlap_local[at].Y0[count_local] = y0;
                            pLsdft->AtomOverlap_local[at].Z0[count_local] = z0;

                            pLsdft->AtomOverlap_local[at].xstart[count_local] = xstart;
                            pLsdft->AtomOverlap_local[at].ystart[count_local] = ystart;
                            pLsdft->AtomOverlap_local[at].zstart[count_local] = zstart;

                            pLsdft->AtomOverlap_local[at].xend[count_local] = xend;
                            pLsdft->AtomOverlap_local[at].yend[count_local] = yend;
                            pLsdft->AtomOverlap_local[at].zend[count_local] = zend;

                            pLsdft->AtomOverlap_local[at].xindex[count_local] = posindex - 3;
                            pLsdft->AtomOverlap_local[at].yindex[count_local] = posindex - 2;
                            pLsdft->AtomOverlap_local[at].zindex[count_local] = posindex - 1;

                            count_local++;
                            

                        }
                    }
        }
    }

   
     VecRestoreArray(pLsdft->AtomposSC, &pAtompos);

}

////////////////////////////////////////////////////////////////////////////////////////////////
//    GetInfluencingAtoms: Gets the list of atoms that influence the domain of the processor  //
//                               additionaly stores details of the overlap                    //
////////////////////////////////////////////////////////////////////////////////////////////////
void GetInfluencingAtoms_Nonlocal(LSDFT_OBJ *pLsdft)
{

    //int count_local=0;
    int count_nonlocal = 0;
    int NatomDomain = 0;
    int start, end, offset_nonlocal_x, offset_nonlocal_y, offset_nonlocal_z, index = 0, posindex = 0, at, poscnt;
    PetscInt xcor, ycor, zcor, lxdim, lydim, lzdim;
    PetscInt xs, ys, zs, xl, yl, zl, xs_nl, ys_nl, zs_nl, xl_nl, yl_nl, zl_nl, xi, yj, zk, xstart = -1, ystart =
        -1, zstart = -1, xend = -1, yend = -1, zend = -1, overlap_local = 0, overlap_nonlocal = 0;
    PetscInt xstart_nl = -1, ystart_nl = -1, zstart_nl = -1, xend_nl = -1, yend_nl = -1, zend_nl = -1;
    PetscScalar cutoffr_nonlocal, max1, max2, cutoffr_local_x, cutoffr_local_y, cutoffr_local_z;
    int Imax_x = 0, Imin_x = 0, Imax_y = 0, Imin_y = 0, Imax_z = 0, Imin_z = 0;
    int PP, QQ, RR;

    PetscScalar x0, y0, z0, X0, Y0, Z0;
    PetscScalar *pAtompos;

    PetscScalar R_x = pLsdft->range_x;
    PetscScalar R_y = pLsdft->range_y;
    PetscScalar R_z = pLsdft->range_z;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // DMDAGetCorners(pLsdft->daloccell, &xcor, &ycor, &zcor, &lxdim, &lydim, &lzdim);
    // printf("rank=%d, daloccell: xcor=%d, ycor=%d, zcor=%d, lxdim=%d, lydim=%d, lzdim=%d \n", rank, xcor, ycor, zcor, lxdim,
    //     lydim, lzdim);

    DMDAGetCorners(pLsdft->da, &xcor, &ycor, &zcor, &lxdim, &lydim, &lzdim);

    xcor=xcor-ceil(pLsdft->Rcut/pLsdft->delta_x);
    ycor=ycor-ceil(pLsdft->Rcut/pLsdft->delta_y);
    zcor=zcor-ceil(pLsdft->Rcut/pLsdft->delta_z);
    
    lxdim=lxdim+2*ceil(pLsdft->Rcut/pLsdft->delta_x);
    lydim=lydim+2*ceil(pLsdft->Rcut/pLsdft->delta_y);
    lzdim=lzdim+2*ceil(pLsdft->Rcut/pLsdft->delta_z);
    
   
    //    printf("NL: xcor=%d,ycor=%d,zcor=%d, lxdim=%d,lydim=%d,lzdim=%d \n",xcor,ycor,zcor,lxdim,lydim,lzdim);

    /*
     * allocate memory for Atom overlap object in the structure
     */
    
    pLsdft->AtomOverlap_nonlocal = (OVERLAP_OBJ *) malloc(sizeof(OVERLAP_OBJ) * pLsdft->Ntype);
   
    assert(pLsdft->AtomOverlap_nonlocal != NULL);

    /*
     * first go over the list of all atoms and determine the number of atoms which have overlap
     */

    
     VecGetArray(pLsdft->AtomposSC, &pAtompos);

    /*
     * loop over different types of atoms
     */
    for (at = 0; at < pLsdft->Ntype; at++) {
        cutoffr_local_x = pLsdft->CUTOFF_x[at];
        cutoffr_local_y = pLsdft->CUTOFF_y[at];
        cutoffr_local_z = pLsdft->CUTOFF_z[at];

        //  offset_local_x = (int)ceil(cutoffr_local_x/pLsdft->delta_x+0.5);
        //offset_local_y = (int)ceil(cutoffr_local_y/pLsdft->delta_y+0.5);
        //offset_local_z = (int)ceil(cutoffr_local_z/pLsdft->delta_z+0.5);

        max1 = pLsdft->psd[at].rc[0] > pLsdft->psd[at].rc[1] ? pLsdft->psd[at].rc[0] : pLsdft->psd[at].rc[1];
        max2 = pLsdft->psd[at].rc[2] > pLsdft->psd[at].rc[3] ? pLsdft->psd[at].rc[2] : pLsdft->psd[at].rc[3];
        cutoffr_nonlocal = max1 > max2 ? max1 : max2;
        offset_nonlocal_x = ceil(cutoffr_nonlocal / pLsdft->delta_x);
        offset_nonlocal_y = ceil(cutoffr_nonlocal / pLsdft->delta_y);
        offset_nonlocal_z = ceil(cutoffr_nonlocal / pLsdft->delta_z);

	//	printf("offset_nonlocal_x=%d, offset_nonlocal_y=%d, offset_nonlocal_z=%d\n",offset_nonlocal_x,offset_nonlocal_y,offset_nonlocal_z);

	//  start = (int) floor(pLsdft->startPos[at] / 3); 
	// end = (int) floor(pLsdft->endPos[at] / 3);
	  start = (int) floor(pLsdft->startPosSC[at] / 3);
        end = (int) floor(pLsdft->endPosSC[at] / 3);

        // count_local=0;
        count_nonlocal = 0;

        Imax_x = 0;
        Imin_x = 0;
        Imax_y = 0;
        Imin_y = 0;
        Imax_z = 0;
        Imin_z = 0;

        if (pLsdft->BC == 1)    // 3D periodic
        {
            Imax_x = ceil(cutoffr_local_x / R_x);
            Imin_x = -ceil(cutoffr_local_x / R_x);

            Imax_y = ceil(cutoffr_local_y / R_y);
            Imin_y = -ceil(cutoffr_local_y / R_y);

            Imax_z = ceil(cutoffr_local_z / R_z);
            Imin_z = -ceil(cutoffr_local_z / R_z);
        } else if (pLsdft->BC == 3) // 2D periodic
        {
            Imax_x = ceil(cutoffr_local_x / R_x);
            Imin_x = -ceil(cutoffr_local_x / R_x);

            Imax_y = ceil(cutoffr_local_y / R_y);
            Imin_y = -ceil(cutoffr_local_y / R_y);
        } else if (pLsdft->BC == 4) // 1D periodic
        {
            Imax_x = ceil(cutoffr_local_x / R_x);
            Imin_x = -ceil(cutoffr_local_x / R_x);
        }


        /*
         * loop over every atom of a given type
         */
        for (poscnt = start; poscnt <= end; poscnt++) {
            X0 = pAtompos[index++];
            Y0 = pAtompos[index++];
            Z0 = pAtompos[index++];

            for (PP = Imin_x; PP <= Imax_x; PP++)
                for (QQ = Imin_y; QQ <= Imax_y; QQ++)
                    for (RR = Imin_z; RR <= Imax_z; RR++) {

                        /*
                         * periodic map of atomic position
                         */
                        x0 = X0 + PP * 2.0 * R_x;
                        y0 = Y0 + QQ * 2.0 * R_y;
                        z0 = Z0 + RR * 2.0 * R_z;

                        // xi = roundf((x0 + R_x) / pLsdft->delta_x);
                        // yj = roundf((y0 + R_y) / pLsdft->delta_y);
                        // zk = roundf((z0 + R_z) / pLsdft->delta_z);

		        xi = roundf((x0 + R_x + 0*pLsdft->Rcut) / pLsdft->delta_x);
                        yj = roundf((y0 + R_y + 0*pLsdft->Rcut) / pLsdft->delta_y);
                        zk = roundf((z0 + R_z + 0*pLsdft->Rcut) / pLsdft->delta_z);
			

                        // xs = xi-offset_local_x; xl = xi+offset_local_x;
                        // ys = yj-offset_local_y; yl = yj+offset_local_y;
                        // zs = zk-offset_local_z; zl = zk+offset_local_z;

                        // /*
                        //  * find if domain of influence of pseudocharge overlaps with the domain stored
                        //  * by processor
                        //  */
                        // overlap_local = Check_Overlap(xcor,ycor,zcor,lxdim,lydim,lzdim,xs,ys,zs,xl,yl,zl,&xstart,&ystart,&zstart,&xend,&yend,&zend);
                        // /*
                        //  * increment counter for pseudocharge
                        //  */
                        // count_local +=overlap_local;
                        /*
                         * check overlap for nonlocal
                         */
                        if (pLsdft->psd[at].lmax != 0) {
                            xs_nl = xi - offset_nonlocal_x;
                            xl_nl = xi + offset_nonlocal_x;
                            ys_nl = yj - offset_nonlocal_y;
                            yl_nl = yj + offset_nonlocal_y;
                            zs_nl = zk - offset_nonlocal_z;
                            zl_nl = zk + offset_nonlocal_z;
                            overlap_nonlocal =
                                Check_Overlap(xcor, ycor, zcor, lxdim, lydim, lzdim, xs_nl, ys_nl, zs_nl, xl_nl, yl_nl,
                                zl_nl, &xstart_nl, &ystart_nl, &zstart_nl, &xend_nl, &yend_nl, &zend_nl);
                            count_nonlocal += overlap_nonlocal;
                        }
                    }
        }

        //  pLsdft->AtomOverlap_local[at].Natoms = count_local;
        pLsdft->AtomOverlap_nonlocal[at].Natoms = count_nonlocal;
	//    printf("rank=%d, pLsdft->AtomOverlap_nonlocal[at].Natoms=%d \n", rank, pLsdft->AtomOverlap_nonlocal[at].Natoms);


       

        pLsdft->AtomOverlap_nonlocal[at].X0 = (double *) malloc(sizeof(double) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocal[at].Y0 = (double *) malloc(sizeof(double) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocal[at].Z0 = (double *) malloc(sizeof(double) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocal[at].xstart = (int *) malloc(sizeof(int) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocal[at].ystart = (int *) malloc(sizeof(int) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocal[at].zstart = (int *) malloc(sizeof(int) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocal[at].xend = (int *) malloc(sizeof(int) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocal[at].yend = (int *) malloc(sizeof(int) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocal[at].zend = (int *) malloc(sizeof(int) * count_nonlocal);

        pLsdft->AtomOverlap_nonlocal[at].xindex = (int *) malloc(sizeof(int) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocal[at].yindex = (int *) malloc(sizeof(int) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocal[at].zindex = (int *) malloc(sizeof(int) * count_nonlocal);

      
        count_nonlocal = 0;

        /*
         * loop over every atom of a given type
         */
        for (poscnt = start; poscnt <= end; poscnt++) {
            X0 = pAtompos[posindex++];
            Y0 = pAtompos[posindex++];
            Z0 = pAtompos[posindex++];

            for (PP = Imin_x; PP <= Imax_x; PP++)
                for (QQ = Imin_y; QQ <= Imax_y; QQ++)
                    for (RR = Imin_z; RR <= Imax_z; RR++) {
                        /*
                         * periodic map of atomic position
                         */
                        x0 = X0 + PP * 2.0 * R_x;
                        y0 = Y0 + QQ * 2.0 * R_y;
                        z0 = Z0 + RR * 2.0 * R_z;

                       
                        /*
                         * check overlap for nonlocal
                         */
                        if (pLsdft->psd[at].lmax != 0) {
                       

			xi = roundf((x0 + R_x + 0*pLsdft->Rcut) / pLsdft->delta_x);
                        yj = roundf((y0 + R_y + 0*pLsdft->Rcut) / pLsdft->delta_y);
                        zk = roundf((z0 + R_z + 0*pLsdft->Rcut) / pLsdft->delta_z);
			    

                            xs_nl = xi - offset_nonlocal_x;
                            xl_nl = xi + offset_nonlocal_x;
                            ys_nl = yj - offset_nonlocal_y;
                            yl_nl = yj + offset_nonlocal_y;
                            zs_nl = zk - offset_nonlocal_z;
                            zl_nl = zk + offset_nonlocal_z;

                            //  printf("offset_nonlocal_x=%d \n",offset_nonlocal_x);

                            overlap_nonlocal =
                                Check_Overlap(xcor, ycor, zcor, lxdim, lydim, lzdim, xs_nl, ys_nl, zs_nl, xl_nl, yl_nl,
                                zl_nl, &xstart_nl, &ystart_nl, &zstart_nl, &xend_nl, &yend_nl, &zend_nl);
                            // printf("*rank=%d, x0=%lf, y0=%lf, z0=%lf, X0=%lf, Y0=%lf, Z0=%lf, xstart_nl=%d, ystart_nl=%d, zstart_nl=%d, xend_nl=%d, yend_nl=%d, zend_nl=%d",rank,x0,y0,z0,X0,Y0,Z0,xstart_nl,ystart_nl,zstart_nl,xend_nl,yend_nl,zend_nl);
                            if (overlap_nonlocal == 1) {
                                /*
                                 * store data in structure
                                 */
                                pLsdft->AtomOverlap_nonlocal[at].X0[count_nonlocal] = x0;
                                pLsdft->AtomOverlap_nonlocal[at].Y0[count_nonlocal] = y0;
                                pLsdft->AtomOverlap_nonlocal[at].Z0[count_nonlocal] = z0;

                                pLsdft->AtomOverlap_nonlocal[at].xstart[count_nonlocal] = xstart_nl;
                                pLsdft->AtomOverlap_nonlocal[at].ystart[count_nonlocal] = ystart_nl;
                                pLsdft->AtomOverlap_nonlocal[at].zstart[count_nonlocal] = zstart_nl;

                                pLsdft->AtomOverlap_nonlocal[at].xend[count_nonlocal] = xend_nl;
                                pLsdft->AtomOverlap_nonlocal[at].yend[count_nonlocal] = yend_nl;
                                pLsdft->AtomOverlap_nonlocal[at].zend[count_nonlocal] = zend_nl;

                                pLsdft->AtomOverlap_nonlocal[at].xindex[count_nonlocal] = posindex - 3;
                                pLsdft->AtomOverlap_nonlocal[at].yindex[count_nonlocal] = posindex - 2;
                                pLsdft->AtomOverlap_nonlocal[at].zindex[count_nonlocal] = posindex - 1;

				// printf("rank=%d, x0=%lf, y0=%lf, z0=%lf, xi=%d, yk=%d, zk=%d, xstart_nl=%d, ystart_nl=%d, zstart_nl=%d, xend_nl=%d, yend_nl=%d, zend_nl=%d \n",
				//       rank,x0,y0,z0,xi,yj,zk,xstart_nl,ystart_nl,zstart_nl,xend_nl,yend_nl,zend_nl);

                                count_nonlocal++;
                            }
                        }
                       
                    }
        }
    }

   
     VecRestoreArray(pLsdft->AtomposSC, &pAtompos);
}

////////////////////////////////////////////////////////////////////////////////////////////////
//    GetInfluencingAtoms: Gets the list of atoms that influence the domain of the processor  //
//                               additionaly stores details of the overlap                    //
////////////////////////////////////////////////////////////////////////////////////////////////
void GetInfluencingAtoms_LocalForces(LSDFT_OBJ *pLsdft)
{

    int count_local = 0;
    //int count_nonlocal=0;
    int NatomDomain = 0;
    int start, end, offset_local_x, offset_nonlocal_x, offset_local_y, offset_nonlocal_y, offset_local_z,
        offset_nonlocal_z, index = 0, posindex = 0, at, poscnt;
    PetscInt xcor, ycor, zcor, lxdim, lydim, lzdim;
    PetscInt xs, ys, zs, xl, yl, zl, xs_nl, ys_nl, zs_nl, xl_nl, yl_nl, zl_nl, xi, yj, zk, xstart = -1, ystart =
        -1, zstart = -1, xend = -1, yend = -1, zend = -1, overlap_local = 0, overlap_nonlocal = 0;
    PetscInt xstart_nl = -1, ystart_nl = -1, zstart_nl = -1, xend_nl = -1, yend_nl = -1, zend_nl = -1;
    PetscScalar cutoffr_nonlocal, cutoffr_local_x, cutoffr_local_y, cutoffr_local_z, max1, max2;
    int Imax_x = 0, Imin_x = 0, Imax_y = 0, Imin_y = 0, Imax_z = 0, Imin_z = 0;
    int PP, QQ, RR;

    PetscScalar x0, y0, z0, X0, Y0, Z0;
    PetscScalar *pAtompos;
    PetscInt XS,XE,YS,YE,ZS,ZE;

    PetscScalar R_x = pLsdft->range_x;
    PetscScalar R_y = pLsdft->range_y;
    PetscScalar R_z = pLsdft->range_z;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    DMDAGetCorners(pLsdft->da, &xcor, &ycor, &zcor, &lxdim, &lydim, &lzdim);
    //xcor=0;ycor=0;zcor=0;lxdim=pLsdft->numPoints_x;lydim=pLsdft->numPoints_y;lzdim=pLsdft->numPoints_z;
    XS=xcor; XE=xcor+lxdim; YS=ycor; YE=ycor+lydim; ZS=zcor; ZE=zcor+lzdim;


    if(xcor==0)
      {
    	// xcor=xcor-ceil(pLsdft->Rcut/pLsdft->delta_x);
    	XS=XS-ceil(pLsdft->Rcut/pLsdft->delta_x);
	
      }
    if(ycor==0)
      {
    	// ycor=ycor-ceil(pLsdft->Rcut/pLsdft->delta_y);
    	YS=YS-ceil(pLsdft->Rcut/pLsdft->delta_y);

      }
    if(zcor==0)
      {
    	// zcor=zcor-ceil(pLsdft->Rcut/pLsdft->delta_z);
    	ZS=ZS-ceil(pLsdft->Rcut/pLsdft->delta_z);
      }
    if(xcor+lxdim==pLsdft->numPoints_x)
      {
    // lxdim=lxdim+ceil(pLsdft->Rcut/pLsdft->delta_x);
    	XE=XE+ceil(pLsdft->Rcut/pLsdft->delta_x);
      }
    if(ycor+lydim==pLsdft->numPoints_y)
      {
    // lydim=lydim+ceil(pLsdft->Rcut/pLsdft->delta_y);
       YE=YE+ceil(pLsdft->Rcut/pLsdft->delta_y);
      }
     if(zcor+lzdim==pLsdft->numPoints_z)
      {
    // lzdim=lzdim+ceil(pLsdft->Rcut/pLsdft->delta_z);
       ZE=ZE+ceil(pLsdft->Rcut/pLsdft->delta_z);
      }

    
     // printf("Localforce, rank=%d,XS=%d,YS=%d,ZS=%d,XE=%d,YE=%d,ZE=%d \n",rank,XS,YS,ZS,XE,YE,ZE);
     // printf("lxdim=%d, lydim=%d, lzdim=%d \n",lxdim,lydim,lzdim);
     xcor=XS;lxdim=XE-XS;
     ycor=YS;lydim=YE-YS; 
     zcor=ZS;lzdim=ZE-ZS; 

     // printf("Localforce, rank=%d, lxdim=%d, lydim=%d, lzdim=%d \n",rank,lxdim,lydim,lzdim);


    /*
     * allocate memory for Atom overlap object in the structure
     */
    pLsdft->AtomOverlap_localforces = (OVERLAP_OBJ *) malloc(sizeof(OVERLAP_OBJ) * pLsdft->Ntype);
   
    assert(pLsdft->AtomOverlap_localforces != NULL);
    

    /*
     * first go over the list of all atoms and determine the number of atoms which have overlap
     */

     VecGetArray(pLsdft->Atompos, &pAtompos); 
     // VecGetArray(pLsdft->AtomposSC, &pAtompos); // changed

    /*
     * loop over different types of atoms
     */
    for (at = 0; at < pLsdft->Ntype; at++) {
        cutoffr_local_x = pLsdft->CUTOFF_x[at];
        cutoffr_local_y = pLsdft->CUTOFF_y[at];
        cutoffr_local_z = pLsdft->CUTOFF_z[at];

        offset_local_x = (int) ceil(cutoffr_local_x / pLsdft->delta_x + 0.5);
        offset_local_y = (int) ceil(cutoffr_local_y / pLsdft->delta_y + 0.5);
        offset_local_z = (int) ceil(cutoffr_local_z / pLsdft->delta_z + 0.5);

        // max1=  pLsdft->psd[at].rc[0] > pLsdft->psd[at].rc[1] ? pLsdft->psd[at].rc[0]:pLsdft->psd[at].rc[1];
        // max2=  pLsdft->psd[at].rc[2] > pLsdft->psd[at].rc[3] ? pLsdft->psd[at].rc[2]:pLsdft->psd[at].rc[3];
        // cutoffr_nonlocal =  max1>max2 ? max1:max2;
        // offset_nonlocal_x = ceil(cutoffr_nonlocal/pLsdft->delta_x);
        // offset_nonlocal_y = ceil(cutoffr_nonlocal/pLsdft->delta_y);
        // offset_nonlocal_z = ceil(cutoffr_nonlocal/pLsdft->delta_z);

	    start = (int) floor(pLsdft->startPos[at] / 3);
	 end = (int) floor(pLsdft->endPos[at] / 3);


	//	start = (int) floor(pLsdft->startPosSC[at] / 3);
	// end = (int) floor(pLsdft->endPosSC[at] / 3);

        count_local = 0;
        // count_nonlocal=0;

        Imax_x = 0;
        Imin_x = 0;
        Imax_y = 0;
        Imin_y = 0;
        Imax_z = 0;
        Imin_z = 0;

        if (pLsdft->BC == 2)    // 3D periodic
        {
            Imax_x = ceil(cutoffr_local_x / R_x);
            Imin_x = -ceil(cutoffr_local_x / R_x);

            Imax_y = ceil(cutoffr_local_y / R_y);
            Imin_y = -ceil(cutoffr_local_y / R_y);

            Imax_z = ceil(cutoffr_local_z / R_z);
            Imin_z = -ceil(cutoffr_local_z / R_z);
        } else if (pLsdft->BC == 3) // 2D periodic
        {
            Imax_x = ceil(cutoffr_local_x / R_x);
            Imin_x = -ceil(cutoffr_local_x / R_x);

            Imax_y = ceil(cutoffr_local_y / R_y);
            Imin_y = -ceil(cutoffr_local_y / R_y);
        } else if (pLsdft->BC == 4) // 1D periodic
        {
            Imax_x = ceil(cutoffr_local_x / R_x);
            Imin_x = -ceil(cutoffr_local_x / R_x);
        }


        /*
         * loop over every atom of a given type
         */
        for (poscnt = start; poscnt <= end; poscnt++) {
            X0 = pAtompos[index++];
            Y0 = pAtompos[index++];
            Z0 = pAtompos[index++];

            for (PP = Imin_x; PP <= Imax_x; PP++)
                for (QQ = Imin_y; QQ <= Imax_y; QQ++)
                    for (RR = Imin_z; RR <= Imax_z; RR++) {

                        /*
                         * periodic map of atomic position
                         */
                        x0 = X0 + PP * 2.0 * R_x;
                        y0 = Y0 + QQ * 2.0 * R_y;
                        z0 = Z0 + RR * 2.0 * R_z;

                        xi = roundf((x0 + R_x) / pLsdft->delta_x);
                        yj = roundf((y0 + R_y) / pLsdft->delta_y);
                        zk = roundf((z0 + R_z) / pLsdft->delta_z);

                        xs = xi - offset_local_x;
                        xl = xi + offset_local_x;
                        ys = yj - offset_local_y;
                        yl = yj + offset_local_y;
                        zs = zk - offset_local_z;
                        zl = zk + offset_local_z;

                        /*
                         * find if domain of influence of pseudocharge overlaps with the domain stored
                         * by processor
                         */
                        overlap_local =
                            Check_Overlap(xcor, ycor, zcor, lxdim, lydim, lzdim, xs, ys, zs, xl, yl, zl, &xstart,
                            &ystart, &zstart, &xend, &yend, &zend);
                        /*
                         * increment counter for pseudocharge
                         */
                        count_local += overlap_local;
                        /*
                         * check overlap for nonlocal
                         */
                        // if(overlap_local==1 && pLsdft->psd[at].lmax!=0)
                        //   {
                        //     xs_nl = xi-offset_nonlocal_x; xl_nl = xi+offset_nonlocal_x;
                        //     ys_nl = yj-offset_nonlocal_y; yl_nl = yj+offset_nonlocal_y;
                        //     zs_nl = zk-offset_nonlocal_z; zl_nl = zk+offset_nonlocal_z;

                        //     overlap_nonlocal = Check_Overlap(xcor,ycor,zcor,lxdim,lydim,lzdim,xs_nl,ys_nl,zs_nl,xl_nl,yl_nl,zl_nl,&xstart_nl,&ystart_nl,&zstart_nl,&xend_nl,&yend_nl,&zend_nl);
                        //     count_nonlocal +=overlap_nonlocal;
                        //   }
                    }
        }

        pLsdft->AtomOverlap_localforces[at].Natoms = count_local;
        //pLsdft->AtomOverlap_nonlocal[at].Natoms = count_nonlocal;

	//	  printf("LForces rank=%d, pLsdft->AtomOverlap_localforces[at].Natoms=%d \n", rank, pLsdft->AtomOverlap_localforces[at].Natoms);


        /*
         * once the number of local and nonlocal influencing atoms are determined,
         * we allocate memory for the datastructure and store the positions, starting and ending indices of local regions
         */
        pLsdft->AtomOverlap_localforces[at].X0 = (double *) malloc(sizeof(double) * count_local);
        pLsdft->AtomOverlap_localforces[at].Y0 = (double *) malloc(sizeof(double) * count_local);
        pLsdft->AtomOverlap_localforces[at].Z0 = (double *) malloc(sizeof(double) * count_local);
        pLsdft->AtomOverlap_localforces[at].xstart = (int *) malloc(sizeof(int) * count_local);
        pLsdft->AtomOverlap_localforces[at].ystart = (int *) malloc(sizeof(int) * count_local);
        pLsdft->AtomOverlap_localforces[at].zstart = (int *) malloc(sizeof(int) * count_local);
        pLsdft->AtomOverlap_localforces[at].xend = (int *) malloc(sizeof(int) * count_local);
        pLsdft->AtomOverlap_localforces[at].yend = (int *) malloc(sizeof(int) * count_local);
        pLsdft->AtomOverlap_localforces[at].zend = (int *) malloc(sizeof(int) * count_local);

        pLsdft->AtomOverlap_localforces[at].xindex = (int *) malloc(sizeof(int) * count_local);
        pLsdft->AtomOverlap_localforces[at].yindex = (int *) malloc(sizeof(int) * count_local);
        pLsdft->AtomOverlap_localforces[at].zindex = (int *) malloc(sizeof(int) * count_local);

        count_local = 0;
        // count_nonlocal=0;

        /*
         * loop over every atom of a given type
         */
        for (poscnt = start; poscnt <= end; poscnt++) {
            X0 = pAtompos[posindex++];
            Y0 = pAtompos[posindex++];
            Z0 = pAtompos[posindex++];

            for (PP = Imin_x; PP <= Imax_x; PP++)
                for (QQ = Imin_y; QQ <= Imax_y; QQ++)
                    for (RR = Imin_z; RR <= Imax_z; RR++) {
                        /*
                         * periodic map of atomic position
                         */
                        x0 = X0 + PP * 2.0 * R_x;
                        y0 = Y0 + QQ * 2.0 * R_y;
                        z0 = Z0 + RR * 2.0 * R_z;

                        xi = roundf((x0 + R_x) / pLsdft->delta_x);
                        yj = roundf((y0 + R_y) / pLsdft->delta_y);
                        zk = roundf((z0 + R_z) / pLsdft->delta_z);
	    
                        xs = xi - offset_local_x;
                        xl = xi + offset_local_x;
                        ys = yj - offset_local_y;
                        yl = yj + offset_local_y;
                        zs = zk - offset_local_z;
                        zl = zk + offset_local_z;

                        overlap_local =
                            Check_Overlap(xcor, ycor, zcor, lxdim, lydim, lzdim, xs, ys, zs, xl, yl, zl, &xstart,
                            &ystart, &zstart, &xend, &yend, &zend);

			// printf("offset_local_x=%d \n",offset_local_x);
			// printf("xcor=%d,lxdim=%d,xs=%d,xl=%d,xstart=%d,xend=%d \n",xcor,lxdim,xs,xl,xstart,xend);
			// printf("ycor=%d,lydim=%d,ys=%d,yl=%d,ystart=%d,yend=%d \n",ycor,lydim,ys,yl,ystart,yend);
			// printf("zcor=%d,lzdim=%d,zs=%d,zl=%d,zstart=%d,zend=%d \n",zcor,lzdim,zs,zl,zstart,zend);


                        if (overlap_local == 1) {
                            /*
                             * store data in structure
                             */
                            pLsdft->AtomOverlap_localforces[at].X0[count_local] = x0;
                            pLsdft->AtomOverlap_localforces[at].Y0[count_local] = y0;
                            pLsdft->AtomOverlap_localforces[at].Z0[count_local] = z0;

                            pLsdft->AtomOverlap_localforces[at].xstart[count_local] = xstart;
                            pLsdft->AtomOverlap_localforces[at].ystart[count_local] = ystart;
                            pLsdft->AtomOverlap_localforces[at].zstart[count_local] = zstart;

                            pLsdft->AtomOverlap_localforces[at].xend[count_local] = xend;
                            pLsdft->AtomOverlap_localforces[at].yend[count_local] = yend;
                            pLsdft->AtomOverlap_localforces[at].zend[count_local] = zend;

                            pLsdft->AtomOverlap_localforces[at].xindex[count_local] = posindex - 3;
                            pLsdft->AtomOverlap_localforces[at].yindex[count_local] = posindex - 2;
                            pLsdft->AtomOverlap_localforces[at].zindex[count_local] = posindex - 1;

                            count_local++;
                      

                        }
                    }
        }
    }

    VecRestoreArray(pLsdft->Atompos, &pAtompos);
  

}

////////////////////////////////////////////////////////////////////////////////////////////////
//    GetInfluencingAtoms: Gets the list of atoms that influence the domain of the processor  //
//                               additionaly stores details of the overlap                    //
////////////////////////////////////////////////////////////////////////////////////////////////
void GetInfluencingAtoms_NonlocalForces(LSDFT_OBJ *pLsdft)
{

    //int count_local=0;
    int count_nonlocal = 0;
    int NatomDomain = 0;
    int start, end, offset_nonlocal_x, offset_nonlocal_y, offset_nonlocal_z, index = 0, posindex = 0, at, poscnt;
    PetscInt xcor, ycor, zcor, lxdim, lydim, lzdim;
    PetscInt xs, ys, zs, xl, yl, zl, xs_nl, ys_nl, zs_nl, xl_nl, yl_nl, zl_nl, xi, yj, zk, xstart = -1, ystart =
        -1, zstart = -1, xend = -1, yend = -1, zend = -1, overlap_local = 0, overlap_nonlocal = 0;
    PetscInt xstart_nl = -1, ystart_nl = -1, zstart_nl = -1, xend_nl = -1, yend_nl = -1, zend_nl = -1;
    PetscScalar cutoffr_nonlocal, max1, max2, cutoffr_local_x, cutoffr_local_y, cutoffr_local_z;
    int Imax_x = 0, Imin_x = 0, Imax_y = 0, Imin_y = 0, Imax_z = 0, Imin_z = 0;
    int PP, QQ, RR;

    PetscScalar x0, y0, z0, X0, Y0, Z0;
    PetscScalar *pAtompos;
      PetscInt XS,XE,YS,YE,ZS,ZE;

    PetscScalar R_x = pLsdft->range_x;
    PetscScalar R_y = pLsdft->range_y;
    PetscScalar R_z = pLsdft->range_z;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // DMDAGetCorners(pLsdft->daloccell, &xcor, &ycor, &zcor, &lxdim, &lydim, &lzdim);
    // printf("rank=%d, daloccell: xcor=%d, ycor=%d, zcor=%d, lxdim=%d, lydim=%d, lzdim=%d \n", rank, xcor, ycor, zcor, lxdim,
    //     lydim, lzdim);

    DMDAGetCorners(pLsdft->da, &xcor, &ycor, &zcor, &lxdim, &lydim, &lzdim);
   
     // if(xcor==0)
    //   {
    // 	xcor=xcor-ceil(pLsdft->Rcut/pLsdft->delta_x);
    //   }
    // if(ycor=0)
    //   {
    // 	ycor=ycor-ceil(pLsdft->Rcut/pLsdft->delta_y);
    //   }
    // if(zcor==0)
    //   {
    // 	zcor=zcor-ceil(pLsdft->Rcut/pLsdft->delta_z);
    //   }
    // if(lxdim==pLsdft->numPoints_x)
    //   {
    // lxdim=lxdim+ceil(pLsdft->Rcut/pLsdft->delta_x);
    //   }
    // if(lydim==pLsdft->numPoints_y)
    //   {
    // lydim=lydim+ceil(pLsdft->Rcut/pLsdft->delta_y);
    //   }
    //  if(lzdim==pLsdft->numPoints_z)
    //   {
    // lzdim=lzdim+ceil(pLsdft->Rcut/pLsdft->delta_z);
    //   }


    //  XS=xcor; XE=xcor+lxdim; YS=ycor; YE=ycor+lydim; ZS=zcor; ZE=zcor+lzdim;

    // if(xcor==0)
    //   {
    // 	xcor=xcor-ceil(pLsdft->Rcut/pLsdft->delta_x);
    // 	XS=XS-ceil(pLsdft->Rcut/pLsdft->delta_x);
	
    //   }
    // if(ycor==0)
    //   {
    // 	ycor=ycor-ceil(pLsdft->Rcut/pLsdft->delta_y);
    // 	YS=YS-ceil(pLsdft->Rcut/pLsdft->delta_y);

    //   }
    // if(zcor==0)
    //   {
    // 	zcor=zcor-ceil(pLsdft->Rcut/pLsdft->delta_z);
    // 	ZS=ZS-ceil(pLsdft->Rcut/pLsdft->delta_z);
    //   }
    // if(lxdim==pLsdft->numPoints_x)
    //   {
    // lxdim=lxdim+ceil(pLsdft->Rcut/pLsdft->delta_x);
    // 	XE=XE+ceil(pLsdft->Rcut/pLsdft->delta_x);
    //   }
    // if(lydim==pLsdft->numPoints_y)
    //   {
    // lydim=lydim+ceil(pLsdft->Rcut/pLsdft->delta_y);
    //    YE=YE+ceil(pLsdft->Rcut/pLsdft->delta_y);
    //   }
    //  if(lzdim==pLsdft->numPoints_z)
    //   {
    // lzdim=lzdim+ceil(pLsdft->Rcut/pLsdft->delta_z);
    //    ZE=ZE+ceil(pLsdft->Rcut/pLsdft->delta_z);
    //   }

      XS=xcor; XE=xcor+lxdim; YS=ycor; YE=ycor+lydim; ZS=zcor; ZE=zcor+lzdim;

     	XS=XS-ceil(pLsdft->Rcut/pLsdft->delta_x);
	YS=YS-ceil(pLsdft->Rcut/pLsdft->delta_y);
	ZS=ZS-ceil(pLsdft->Rcut/pLsdft->delta_z);

	XE=XE+ceil(pLsdft->Rcut/pLsdft->delta_x);
	YE=YE+ceil(pLsdft->Rcut/pLsdft->delta_y);
	ZE=ZE+ceil(pLsdft->Rcut/pLsdft->delta_z);

     // printf("NL: XS=%d,YS=%d,ZS=%d,XE=%d,YE=%d,ZE=%d \n",XS,YS,ZS,XE,YE,ZE);
     // printf("NL: lxdim=%d, lydim=%d, lzdim=%d \n",lxdim,lydim,lzdim);
     xcor=XS;lxdim=XE-XS;
     ycor=YS;lydim=YE-YS; 
     zcor=ZS;lzdim=ZE-ZS; 

     // printf("NL: lxdim=%d, lydim=%d, lzdim=%d \n",lxdim,lydim,lzdim);


     // printf("NLforces: xcor=%d,ycor=%d,zcor=%d, lxdim=%d,lydim=%d,lzdim=%d \n",xcor,ycor,zcor,lxdim,lydim,lzdim);

    /*
     * allocate memory for Atom overlap object in the structure
     */
   
    pLsdft->AtomOverlap_nonlocalforces = (OVERLAP_OBJ *) malloc(sizeof(OVERLAP_OBJ) * pLsdft->Ntype);
   
    assert(pLsdft->AtomOverlap_nonlocalforces != NULL);

    /*
     * first go over the list of all atoms and determine the number of atoms which have overlap
     */

     VecGetArray(pLsdft->Atompos, &pAtompos);
 

    /*
     * loop over different types of atoms
     */
    for (at = 0; at < pLsdft->Ntype; at++) {
        cutoffr_local_x = pLsdft->CUTOFF_x[at];
        cutoffr_local_y = pLsdft->CUTOFF_y[at];
        cutoffr_local_z = pLsdft->CUTOFF_z[at];

        //  offset_local_x = (int)ceil(cutoffr_local_x/pLsdft->delta_x+0.5);
        //offset_local_y = (int)ceil(cutoffr_local_y/pLsdft->delta_y+0.5);
        //offset_local_z = (int)ceil(cutoffr_local_z/pLsdft->delta_z+0.5);

        max1 = pLsdft->psd[at].rc[0] > pLsdft->psd[at].rc[1] ? pLsdft->psd[at].rc[0] : pLsdft->psd[at].rc[1];
        max2 = pLsdft->psd[at].rc[2] > pLsdft->psd[at].rc[3] ? pLsdft->psd[at].rc[2] : pLsdft->psd[at].rc[3];
        cutoffr_nonlocal = max1 > max2 ? max1 : max2;
        offset_nonlocal_x = ceil(cutoffr_nonlocal / pLsdft->delta_x);
        offset_nonlocal_y = ceil(cutoffr_nonlocal / pLsdft->delta_y);
        offset_nonlocal_z = ceil(cutoffr_nonlocal / pLsdft->delta_z);

	// printf("offset_nonlocal_x=%d, offset_nonlocal_y=%d, offset_nonlocal_z=%d\n",offset_nonlocal_x,offset_nonlocal_y,offset_nonlocal_z);

	  start = (int) floor(pLsdft->startPos[at] / 3);
	 end = (int) floor(pLsdft->endPos[at] / 3);
	//	  start = (int) floor(pLsdft->startPosSC[at] / 3);
	// end = (int) floor(pLsdft->endPosSC[at] / 3);

        // count_local=0;
        count_nonlocal = 0;

        Imax_x = 0;
        Imin_x = 0;
        Imax_y = 0;
        Imin_y = 0;
        Imax_z = 0;
        Imin_z = 0;

        if (pLsdft->BC == 2)    // 3D periodic
        {
            Imax_x = ceil(cutoffr_local_x / R_x);
            Imin_x = -ceil(cutoffr_local_x / R_x);

            Imax_y = ceil(cutoffr_local_y / R_y);
            Imin_y = -ceil(cutoffr_local_y / R_y);

            Imax_z = ceil(cutoffr_local_z / R_z);
            Imin_z = -ceil(cutoffr_local_z / R_z);
        } else if (pLsdft->BC == 3) // 2D periodic
        {
            Imax_x = ceil(cutoffr_local_x / R_x);
            Imin_x = -ceil(cutoffr_local_x / R_x);

            Imax_y = ceil(cutoffr_local_y / R_y);
            Imin_y = -ceil(cutoffr_local_y / R_y);
        } else if (pLsdft->BC == 4) // 1D periodic
        {
            Imax_x = ceil(cutoffr_local_x / R_x);
            Imin_x = -ceil(cutoffr_local_x / R_x);
        }


        /*
         * loop over every atom of a given type
         */
        for (poscnt = start; poscnt <= end; poscnt++) {
            X0 = pAtompos[index++];
            Y0 = pAtompos[index++];
            Z0 = pAtompos[index++];

            for (PP = Imin_x; PP <= Imax_x; PP++)
                for (QQ = Imin_y; QQ <= Imax_y; QQ++)
                    for (RR = Imin_z; RR <= Imax_z; RR++) {

                        /*
                         * periodic map of atomic position
                         */
                        x0 = X0 + PP * 2.0 * R_x;
                        y0 = Y0 + QQ * 2.0 * R_y;
                        z0 = Z0 + RR * 2.0 * R_z;

                        // xi = roundf((x0 + R_x) / pLsdft->delta_x);
                        // yj = roundf((y0 + R_y) / pLsdft->delta_y);
                        // zk = roundf((z0 + R_z) / pLsdft->delta_z);

		        xi = roundf((x0 + R_x) / pLsdft->delta_x);
                        yj = roundf((y0 + R_y) / pLsdft->delta_y);
                        zk = roundf((z0 + R_z) / pLsdft->delta_z);
			

                        // xs = xi-offset_local_x; xl = xi+offset_local_x;
                        // ys = yj-offset_local_y; yl = yj+offset_local_y;
                        // zs = zk-offset_local_z; zl = zk+offset_local_z;

                        // /*
                        //  * find if domain of influence of pseudocharge overlaps with the domain stored
                        //  * by processor
                        //  */
                        // overlap_local = Check_Overlap(xcor,ycor,zcor,lxdim,lydim,lzdim,xs,ys,zs,xl,yl,zl,&xstart,&ystart,&zstart,&xend,&yend,&zend);
                        // /*
                        //  * increment counter for pseudocharge
                        //  */
                        // count_local +=overlap_local;
                        /*
                         * check overlap for nonlocal
                         */
                        if (pLsdft->psd[at].lmax != 0) {
                            xs_nl = xi - offset_nonlocal_x;
                            xl_nl = xi + offset_nonlocal_x;
                            ys_nl = yj - offset_nonlocal_y;
                            yl_nl = yj + offset_nonlocal_y;
                            zs_nl = zk - offset_nonlocal_z;
                            zl_nl = zk + offset_nonlocal_z;
                            overlap_nonlocal =
                                Check_Overlap(xcor, ycor, zcor, lxdim, lydim, lzdim, xs_nl, ys_nl, zs_nl, xl_nl, yl_nl,
                                zl_nl, &xstart_nl, &ystart_nl, &zstart_nl, &xend_nl, &yend_nl, &zend_nl);
                            count_nonlocal += overlap_nonlocal;
                        }
                    }
        }

        //  pLsdft->AtomOverlap_local[at].Natoms = count_local;
        pLsdft->AtomOverlap_nonlocalforces[at].Natoms = count_nonlocal;
	//   printf("rank=%d, pLsdft->AtomOverlap_nonlocalforces[at].Natoms=%d \n", rank, pLsdft->AtomOverlap_nonlocalforces[at].Natoms);


      
        pLsdft->AtomOverlap_nonlocalforces[at].X0 = (double *) malloc(sizeof(double) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocalforces[at].Y0 = (double *) malloc(sizeof(double) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocalforces[at].Z0 = (double *) malloc(sizeof(double) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocalforces[at].xstart = (int *) malloc(sizeof(int) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocalforces[at].ystart = (int *) malloc(sizeof(int) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocalforces[at].zstart = (int *) malloc(sizeof(int) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocalforces[at].xend = (int *) malloc(sizeof(int) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocalforces[at].yend = (int *) malloc(sizeof(int) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocalforces[at].zend = (int *) malloc(sizeof(int) * count_nonlocal);

        pLsdft->AtomOverlap_nonlocalforces[at].xindex = (int *) malloc(sizeof(int) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocalforces[at].yindex = (int *) malloc(sizeof(int) * count_nonlocal);
        pLsdft->AtomOverlap_nonlocalforces[at].zindex = (int *) malloc(sizeof(int) * count_nonlocal);

        
        count_nonlocal = 0;

        /*
         * loop over every atom of a given type
         */
        for (poscnt = start; poscnt <= end; poscnt++) {
            X0 = pAtompos[posindex++];
            Y0 = pAtompos[posindex++];
            Z0 = pAtompos[posindex++];

            for (PP = Imin_x; PP <= Imax_x; PP++)
                for (QQ = Imin_y; QQ <= Imax_y; QQ++)
                    for (RR = Imin_z; RR <= Imax_z; RR++) {
                        /*
                         * periodic map of atomic position
                         */
                        x0 = X0 + PP * 2.0 * R_x;
                        y0 = Y0 + QQ * 2.0 * R_y;
                        z0 = Z0 + RR * 2.0 * R_z;

                        

                        /*
                         * check overlap for nonlocal
                         */
                        if (pLsdft->psd[at].lmax != 0) {
                           
			xi = roundf((x0 + R_x) / pLsdft->delta_x);
                        yj = roundf((y0 + R_y) / pLsdft->delta_y);
                        zk = roundf((z0 + R_z) / pLsdft->delta_z);
			    

                            xs_nl = xi - offset_nonlocal_x;
                            xl_nl = xi + offset_nonlocal_x;
                            ys_nl = yj - offset_nonlocal_y;
                            yl_nl = yj + offset_nonlocal_y;
                            zs_nl = zk - offset_nonlocal_z;
                            zl_nl = zk + offset_nonlocal_z;

                            //  printf("offset_nonlocal_x=%d \n",offset_nonlocal_x);

                            overlap_nonlocal =
                                Check_Overlap(xcor, ycor, zcor, lxdim, lydim, lzdim, xs_nl, ys_nl, zs_nl, xl_nl, yl_nl,
                                zl_nl, &xstart_nl, &ystart_nl, &zstart_nl, &xend_nl, &yend_nl, &zend_nl);
                            // printf("*rank=%d, x0=%lf, y0=%lf, z0=%lf, X0=%lf, Y0=%lf, Z0=%lf, xstart_nl=%d, ystart_nl=%d, zstart_nl=%d, xend_nl=%d, yend_nl=%d, zend_nl=%d",rank,x0,y0,z0,X0,Y0,Z0,xstart_nl,ystart_nl,zstart_nl,xend_nl,yend_nl,zend_nl);
                            if (overlap_nonlocal == 1) {
                                /*
                                 * store data in structure
                                 */
                                pLsdft->AtomOverlap_nonlocalforces[at].X0[count_nonlocal] = x0;
                                pLsdft->AtomOverlap_nonlocalforces[at].Y0[count_nonlocal] = y0;
                                pLsdft->AtomOverlap_nonlocalforces[at].Z0[count_nonlocal] = z0;

                                pLsdft->AtomOverlap_nonlocalforces[at].xstart[count_nonlocal] = xstart_nl;
                                pLsdft->AtomOverlap_nonlocalforces[at].ystart[count_nonlocal] = ystart_nl;
                                pLsdft->AtomOverlap_nonlocalforces[at].zstart[count_nonlocal] = zstart_nl;

                                pLsdft->AtomOverlap_nonlocalforces[at].xend[count_nonlocal] = xend_nl;
                                pLsdft->AtomOverlap_nonlocalforces[at].yend[count_nonlocal] = yend_nl;
                                pLsdft->AtomOverlap_nonlocalforces[at].zend[count_nonlocal] = zend_nl;

                                pLsdft->AtomOverlap_nonlocalforces[at].xindex[count_nonlocal] = posindex - 3;
                                pLsdft->AtomOverlap_nonlocalforces[at].yindex[count_nonlocal] = posindex - 2;
                                pLsdft->AtomOverlap_nonlocalforces[at].zindex[count_nonlocal] = posindex - 1;

				       // printf("NLFORCES: rank=%d, x0=%lf, y0=%lf, z0=%lf, xi=%d, yj=%d, zk=%d, xstart_nl=%d, ystart_nl=%d, zstart_nl=%d, xend_nl=%d, yend_nl=%d, zend_nl=%d \n",rank,x0,y0,z0,xi,yj,zk,xstart_nl,ystart_nl,zstart_nl,xend_nl,yend_nl,zend_nl);

                                count_nonlocal++;
                            }
                        }
                      
                    }
        }
    }

    VecRestoreArray(pLsdft->Atompos, &pAtompos);
    
}



////////////////////////////////////////////////////////////////////////////////////////////////
//Check_Overlap: checks if the influence region of an atom overlaps with the processor region //
////////////////////////////////////////////////////////////////////////////////////////////////
int Check_Overlap(int xcor, int ycor, int zcor, int lxdim, int lydim, int lzdim, int xs, int ys, int zs, int xl, int yl,
    int zl, int *xstart, int *ystart, int *zstart, int *xend, int *yend, int *zend)
{
    /*
     * find if domain of influence of pseudocharge overlaps with the domain stored
     * by processor
     */
    *xstart = -1000;
    *ystart = -1000;
    *zstart = -1000;
    *xend = -1000;
    *yend = -1000;
    *zend = -1000;
    int overlap = 0;
    if (xs >= xcor && xs <= xcor + lxdim - 1)
        *xstart = xs;
    else if (xcor >= xs && xcor <= xl)
        *xstart = xcor;

    if (xl >= xcor && xl <= xcor + lxdim - 1)
        *xend = xl;
    else if (xcor + lxdim - 1 >= xs && xcor + lxdim - 1 <= xl)
        *xend = xcor + lxdim - 1;

    if (ys >= ycor && ys <= ycor + lydim - 1)
        *ystart = ys;
    else if (ycor >= ys && ycor <= yl)
        *ystart = ycor;

    if (yl >= ycor && yl <= ycor + lydim - 1)
        *yend = yl;
    else if (ycor + lydim - 1 >= ys && ycor + lydim - 1 <= yl)
        *yend = ycor + lydim - 1;

    if (zs >= zcor && zs <= zcor + lzdim - 1)
        *zstart = zs;
    else if (zcor >= zs && zcor <= zl)
        *zstart = zcor;

    if (zl >= zcor && zl <= zcor + lzdim - 1)
        *zend = zl;
    else if (zcor + lzdim - 1 >= zs && zcor + lzdim - 1 <= zl)
        *zend = zcor + lzdim - 1;

    if ((*xstart != -1000) && (*xend != -1000) && (*ystart != -1000) && (*yend != -1000) && (*zstart != -1000)
        && (*zend != -1000)) {
        overlap = 1;
    } else {
        overlap = 0;
    }
    return overlap;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//                 Create_Veffcommunicator:        //
// written with help from Steve
////////////////////////////////////////////////////////////////////////////////////////////////
PetscErrorCode Create_VeffCommunicaor(LSDFT_OBJ *pLsdft)
{
  PetscInt xcor, ycor, zcor, lxdim, lydim, lzdim;
  //PetscInt *neigh;
  PetscInt Xs,Ys,Zs,Lxdim,Lydim,Lzdim,xstart,ystart,zstart,xend,yend,zend,Lx,Ly,Lz,rowstart,rowend,LxLyLz;
  // PetscInt *pLsdft->ProcDomains;
  PetscMPIInt comm_size;
  PetscErrorCode ierr;
  int rank,count_overlap,overlap,reorder=0,ctr,p;
  MPI_Comm lssgq_comm_dist_graph;
  MPI_Request request;

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &comm_size);
  DMDAGetCorners(pLsdft->da, &xcor, &ycor, &zcor, &lxdim, &lydim, &lzdim);
  VecGetOwnershipRange(pLsdft->Veff,&rowstart,&rowend);
  
  // // allocate memory
  // pLsdft->ProcDomains = (PetscInt *) calloc(comm_size*7,sizeof(PetscInt));
  //pLsdft->OwnershipRange = (PetscInt *) calloc(comm_size*2,sizeof(PetscInt));
  
  // set values
  pLsdft->ProcDomains[rank*7]=rank;
  pLsdft->ProcDomains[rank*7+1]=xcor;
  pLsdft->ProcDomains[rank*7+2]=ycor;
  pLsdft->ProcDomains[rank*7+3]=zcor;
  pLsdft->ProcDomains[rank*7+4]=xcor+lxdim-1;
  pLsdft->ProcDomains[rank*7+5]=ycor+lydim-1;
  pLsdft->ProcDomains[rank*7+6]=zcor+lzdim-1;

  pLsdft->OwnershipRange[rank*2]=rowstart;
  pLsdft->OwnershipRange[rank*2+1]=rowend;

  // Lx=0;Ly=0,Lz=0;
  LxLyLz=0;

  // call all reduce
  ierr = MPI_Allreduce(MPI_IN_PLACE, pLsdft->ProcDomains, comm_size*7,MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
    CHKERRQ(ierr);
     // call all reduce
  ierr = MPI_Allreduce(MPI_IN_PLACE, pLsdft->OwnershipRange, comm_size*2,MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
    CHKERRQ(ierr);

    // define supercell (processor domain + rcut region around each supercell)
    Xs=xcor-ceil(pLsdft->Rcut/pLsdft->delta_x);
    Ys=ycor-ceil(pLsdft->Rcut/pLsdft->delta_y);
    Zs=zcor-ceil(pLsdft->Rcut/pLsdft->delta_z);
    Lxdim=lxdim+2*ceil(pLsdft->Rcut/pLsdft->delta_x);
    Lydim=lydim+2*ceil(pLsdft->Rcut/pLsdft->delta_y);
    Lzdim=lzdim+2*ceil(pLsdft->Rcut/pLsdft->delta_z);
    // Xe=xcor+lxdim+ceil(pLsdft->Rcut/pLsdft->delta_x)-1;
    // Ye=ycor+lydim+ceil(pLsdft->Rcut/pLsdft->delta_y)-1;
    // Ze=zcor+lzdim+ceil(pLsdft->Rcut/pLsdft->delta_z)-1;

    // loop over all the processor domains and check for overlap
    count_overlap=0;ctr=0;
    for(p=0;p<comm_size;p++)
      {
	// note that this part will change when periodicity is introduced
	overlap=0;
	overlap=Check_Overlap(Xs,Ys,Zs,Lxdim,Lydim,Lzdim,pLsdft->ProcDomains[p*7+1],pLsdft->ProcDomains[p*7+2],pLsdft->ProcDomains[p*7+3],pLsdft->ProcDomains[p*7+4],pLsdft->ProcDomains[p*7+5],pLsdft->ProcDomains[p*7+6],&xstart,&ystart,&zstart,&xend,&yend,&zend);

	if(pLsdft->ProcDomains[p*7]==rank)
	  {
	    overlap=0;
	  }
	if(overlap==1)
	  {
	    //Lx=Lx + (pLsdft->ProcDomains[p*7+4]-pLsdft->ProcDomains[p*7+1]+1);
	    //Ly=Ly + (pLsdft->ProcDomains[p*7+5]-pLsdft->ProcDomains[p*7+2]+1);
	    //Lz=Lz + (pLsdft->ProcDomains[p*7+6]-pLsdft->ProcDomains[p*7+3]+1);
	    LxLyLz+=(pLsdft->ProcDomains[p*7+4]-pLsdft->ProcDomains[p*7+1]+1)*(pLsdft->ProcDomains[p*7+5]-pLsdft->ProcDomains[p*7+2]+1)*(pLsdft->ProcDomains[p*7+6]-pLsdft->ProcDomains[p*7+3]+1);
	  }
	count_overlap +=overlap;
      }
    // allocate memory to store neighbors
    pLsdft->neigh = (PetscInt *)calloc(count_overlap,sizeof(PetscInt));
    // allocate memory to store the effective potential coming from neighbors
    printf("rank=%d, LxLyLz=%d \n",rank,LxLyLz);

    pLsdft->Veff_neigh = (PetscScalar *)calloc(LxLyLz,sizeof(PetscScalar));
    // // allocate memory to store the effective potential coming from neighbors
    // pLsdft->VeffNodal = (PetscScalar *)calloc(pLsdft->nx_loc*pLsdft->ny_loc*pLsdft->nz_loc,sizeof(PetscScalar));
    // also allocate memory to store indices
    pLsdft->sendcounts = (PetscInt *)calloc(count_overlap,sizeof(PetscInt));
    pLsdft->sdispls = (PetscInt *)calloc(count_overlap,sizeof(PetscInt));
    pLsdft->recvcounts = (PetscInt *)calloc(count_overlap,sizeof(PetscInt));
    pLsdft->rdispls = (PetscInt *)calloc(count_overlap,sizeof(PetscInt));

    // loop over all the processor domains and store the ranks of overlapping domains
      for(p=0;p<comm_size;p++)
      {
	// note that this part will change when periodicity is introduced
	overlap=0;
	overlap=Check_Overlap(Xs,Ys,Zs,Lxdim,Lydim,Lzdim,pLsdft->ProcDomains[p*7+1],pLsdft->ProcDomains[p*7+2],pLsdft->ProcDomains[p*7+3],pLsdft->ProcDomains[p*7+4],pLsdft->ProcDomains[p*7+5],pLsdft->ProcDomains[p*7+6],&xstart,&ystart,&zstart,&xend,&yend,&zend);

	if(pLsdft->ProcDomains[p*7]==rank)
	  {
	    overlap=0;
	  }

	if(overlap==1)
	  {
	    pLsdft->neigh[ctr]=p;
	    pLsdft->sendcounts[ctr]=lxdim*lydim*lzdim;
	    pLsdft->sdispls[ctr]=0;
	    pLsdft->recvcounts[ctr]=(pLsdft->ProcDomains[p*7+4]-pLsdft->ProcDomains[p*7+1]+1)*(pLsdft->ProcDomains[p*7+5]-pLsdft->ProcDomains[p*7+2]+1)*(pLsdft->ProcDomains[p*7+6]-pLsdft->ProcDomains[p*7+3]+1);
	    if(ctr==0)
	      {
		pLsdft->rdispls[ctr]=0;
	      }else
	      {
		//pLsdft->rdispls[ctr]=pLsdft->recvcounts[ctr-1];
		pLsdft->rdispls[ctr]=pLsdft->rdispls[ctr-1]+pLsdft->recvcounts[ctr-1];
	      }
	    //  printf("@@@ rank=%d,pLsdft->neigh[ctr]=%d,pLsdft->sendcounts[ctr]=%d,pLsdft->sdispls[ctr]=%d,pLsdft->recvcounts[ctr]=%d,pLsdft->rdispls[ctr]=%d \n",rank,pLsdft->neigh[ctr],pLsdft->sendcounts[ctr],pLsdft->sdispls[ctr],pLsdft->recvcounts[ctr],pLsdft->rdispls[ctr]);
	    ctr++;
	   
	  }
      }
      // now create the mpi communicator
       	MPI_Dist_graph_create_adjacent(PETSC_COMM_WORLD,count_overlap,pLsdft->neigh,MPI_UNWEIGHTED,count_overlap,pLsdft->neigh,MPI_UNWEIGHTED,MPI_INFO_NULL,reorder,&lssgq_comm_dist_graph); // creates a distributed graph topology (adjacent, cartesian cubical)

        pLsdft->lssgq_comm_effpot = lssgq_comm_dist_graph;
	pLsdft->nneigh=count_overlap;

	// // print to check
	//   for(p=0;p<count_overlap;p++)
	//     {
	//      
	//     }

	// //////////////////////////////////////////////////////////////////////////
	// // check the communicator by communicating the neighbor ranks
	// PetscMPIInt *neighranks;
	// PetscInt *rr,one=1;
	// PetscInt *sendcounts, *sdispls, *recvcounts, *rdispls;
	

	// neighranks = (PetscMPIInt *)calloc(count_overlap,sizeof(PetscMPIInt));
	// rr = (PetscInt *)calloc(count_overlap,sizeof(PetscInt));
	// sendcounts = (PetscInt *)calloc(count_overlap,sizeof(PetscInt));
	// sdispls = (PetscInt *)calloc(count_overlap,sizeof(PetscInt));
	// recvcounts = (PetscInt *)calloc(count_overlap,sizeof(PetscInt));
	// rdispls = (PetscInt *)calloc(count_overlap,sizeof(PetscInt));

	// /////////////////////////////////////////////////////////////
	// for(p=0;p<count_overlap;p++)
	//   {
	//     rr[p]=rank;
	//     sendcounts[p]=1;
	//     sdispls[p]=0;
	//     recvcounts[p]=1;
	//     rdispls[p]=p;
	    
	//   }
	  
	// //ierr=MPI_Ineighbor_alltoall(rr,1,MPI_INT,neighranks,1,MPI_INT,pLsdft->lssgq_comm_effpot,&request);
	
	// // vector operation
	// ierr=MPI_Ineighbor_alltoallv(&rank,sendcounts,sdispls,MPI_INT,neighranks,recvcounts,rdispls,MPI_INT,pLsdft->lssgq_comm_effpot,&request);

	// MPI_Wait(&request,MPI_STATUS_IGNORE); // ensure above communication has finished
	
	// printf("rank=%d, count_overlap=%d \n",rank,count_overlap);
	//   for(p=0;p<count_overlap;p++)
	//     {
	//       printf("rank=%d,  neigh[%d]=%d,  neighranks[%d]=%d \n",rank,p,neigh[p],p,neighranks[p]);
	//     }
	//	free(neighranks);

	// free arrays
	//	free(neigh);
		
 
  return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//Create_VeffNodalBoundary: boundary Veff values                                              //
////////////////////////////////////////////////////////////////////////////////////////////////
void Create_VeffNodalBoundary(LSDFT_OBJ *pLsdft)
{
  // first go over the Rcut regions and then store the value from the structure
  PetscScalar ***arrVeffBoundary,*arrBVeff;
  PetscInt xcor, ycor, zcor, lxdim, lydim, lzdim;
  PetscInt xcorda, ycorda, zcorda, lxdimda, lydimda, lzdimda;
  PetscInt SCn_x,SCn_y,SCn_z;
  
  int i0,j0,k0,i1,j1,k1,LI,k,j,i,LI2,ii,jj,kk;
  PetscScalar x,y,z,X,Y,Z,dx,dy,dz,x0,y0,z0,x1,y1,z1,hx,hy,hz;
  PetscScalar Lx,Ly,Lz,shift_x,shift_y,shift_z;
  PetscScalar c000,c001,c010,c011,c100,c101,c110,c111;

  Lx=2.0*pLsdft->Brange_x[0];
  Ly=2.0*pLsdft->Brange_y[0];
  Lz=2.0*pLsdft->Brange_z[0];
  
  hx=Lx/pLsdft->BNx[0];
  hy=Ly/pLsdft->BNy[0];
  hz=Lz/pLsdft->BNz[0];

  VecZeroEntries(pLsdft->VeffNodalBoundary);

  DMDAGetCorners(pLsdft->da,&xcorda,&ycorda,&zcorda,&lxdimda,&lydimda,&lzdimda);
  // SCn_x=lxdimda+ 2*ceil(pLsdft->Rcut/pLsdft->delta_x);
  // SCn_y=lydimda+ 2*ceil(pLsdft->Rcut/pLsdft->delta_y);
  // SCn_z=lzdimda+ 2*ceil(pLsdft->Rcut/pLsdft->delta_z);

  shift_x=(ceil(pLsdft->Rcut/pLsdft->delta_x)-xcorda)*pLsdft->delta_x+pLsdft->range_x;
  shift_y=(ceil(pLsdft->Rcut/pLsdft->delta_y)-ycorda)*pLsdft->delta_y+pLsdft->range_y;
  shift_z=(ceil(pLsdft->Rcut/pLsdft->delta_z)-zcorda)*pLsdft->delta_z+pLsdft->range_z;

  // shift_x=(ceil(pLsdft->Rcut/pLsdft->delta_x)+xcorda)*pLsdft->delta_x+0*pLsdft->range_x;
  // shift_y=(ceil(pLsdft->Rcut/pLsdft->delta_y)+ycorda)*pLsdft->delta_y+0*pLsdft->range_y;
  // shift_z=(ceil(pLsdft->Rcut/pLsdft->delta_z)+zcorda)*pLsdft->delta_z+0*pLsdft->range_z;


  //	printf("shift_x=%lf",shift_x);

  DMDAGetCorners(pLsdft->daloccell, &xcor, &ycor, &zcor, &lxdim, &lydim, &lzdim);
  DMDAVecGetArray(pLsdft->daloccell,pLsdft->VeffNodalBoundary,&arrVeffBoundary);
  VecGetArray(pLsdft->BVeff[0],&arrBVeff);

  // now go over all the nodes and check if it is outside the computational domain
  //  printf("lxdim=%d, pLsdft->nx_loc=%d \n",lxdim,pLsdft->nx_loc);
  for(k=zcor;k<zcor+lzdim;k++)
    {
      for(j=ycor;j<ycor+lydim;j++)
	{
	  for(i=xcor;i<xcor+lxdim;i++)
	    {

	      // outside point, now calculate the (x,y,z) coordinate of this file
	      // note the origin is at the center of the computational domain

	      // x=pLsdft->delta_x*i-pLsdft->range_x-pLsdft->Rcut;
	      // y=pLsdft->delta_y*j-pLsdft->range_y-pLsdft->Rcut;
	      // z=pLsdft->delta_z*k-pLsdft->range_z-pLsdft->Rcut;

	      x=pLsdft->delta_x*i-shift_x;
	      y=pLsdft->delta_y*j-shift_y;
	      z=pLsdft->delta_z*k-shift_z;
		  

	      LI2= k*lxdim*lydim+j*lxdim+i;

	      // coordinate of the point with origin at the corner of daloccell 
	      // x=pLsdft->delta_x*i;
	      // y=pLsdft->delta_y*j;
	      // z=pLsdft->delta_z*k;

	      //determine if this node is outside
	      //if((z<-pLsdft->range_z || z>pLsdft->range_z) || (y<-pLsdft->range_y || y>pLsdft->range_y) || (x<-pLsdft->range_x || x>pLsdft->range_x))
	      // if(fabs(z+pLsdft->range_z)<1e-2 || fabs(z-pLsdft->range_z)>1e-2 || fabs(y+pLsdft->range_y)<1e-2 || fabs(y-pLsdft->range_y)>1e-2 || fabs(x+pLsdft->range_x)<1e-2 || fabs(x-pLsdft->range_x)>1e-2)
	      x=round(x*100000000)/100000000;
	      y=round(y*100000000)/100000000;
	      z=round(z*100000000)/100000000;
	      if(z<-pLsdft->range_z || z>pLsdft->range_z || y<-pLsdft->range_y || y>pLsdft->range_y || x<-pLsdft->range_x || x>pLsdft->range_x)
		{
		 
		  // use fmod to calculate the respective point in the boundary info cell
		  // cell is of length Lx,Ly,Lz
		  // probably buggy!!
		  		  
		  if(x<0)
		    dx=x+pLsdft->range_x;
		  else if(x>0)
		    dx=x-pLsdft->range_x;
		  
		  if(y<0)
		    dy=y+pLsdft->range_y;
		  else if(y>0)
		    dy=y-pLsdft->range_y;

		  if(z<0)
		    dz=z+pLsdft->range_z;
		  else if(z>0)
		    dz=z-pLsdft->range_z;


		  X=fmod(dx,Lx);
		  Y=fmod(dy,Ly);
		  Z=fmod(dz,Lz);


		  // probably buggy!!
		  if(X<0)
		    X=X+2.0*pLsdft->Brange_x[0];
		  if(Y<0)
		    Y=Y+2.0*pLsdft->Brange_y[0];
		  if(Z<0)
		    Z=Z+2.0*pLsdft->Brange_z[0];

		  ii=roundf((X+0*pLsdft->Brange_x[0])/pLsdft->delta_x);
		  jj=roundf((Y+0*pLsdft->Brange_y[0])/pLsdft->delta_y);
		  kk=roundf((Z+0*pLsdft->Brange_z[0])/pLsdft->delta_z);
		  
		  if(ii>=pLsdft->BNx[0])
		    ii=ii-pLsdft->BNx[0];
		  if(jj>=pLsdft->BNy[0])
		    jj=jj-pLsdft->BNy[0];
		  if(kk>=pLsdft->BNz[0])
		    kk=kk-pLsdft->BNz[0];

		 
		  LI= kk*pLsdft->BNx[0]*pLsdft->BNy[0]+jj*pLsdft->BNx[0]+ii;
		  arrVeffBoundary[k][j][i]=arrBVeff[LI];

		  // comment from here
		  /*
		  // todo:need to check pointer issue
		  // first find mesh corners
		  FindMeshCorners(x,y,z,pLsdft->range_x,pLsdft->range_y,pLsdft->range_z,pLsdft->Brange_x[0],pLsdft->Brange_y[0],pLsdft->Brange_z[0],pLsdft->BNx[0],pLsdft->BNy[0],pLsdft->BNz[0],&i0,&j0,&k0,&i1,&j1,&k1,&X,&Y,&Z);
		  // next evaluate function values at corners
		  EvaluateCorners(arrBVeff,pLsdft->BNx[0],pLsdft->BNy[0],pLsdft->BNz[0],i0,j0,k0,i1,j1,k1,&c000,&c001,&c010,&c011,&c100,&c101,&c110,&c111);
		  if(pLsdft->LinInterpEFieldsFlag==1)
		    {
		  // need to calculate x0,..., x1,...
		  x0=i0*hx-pLsdft->Brange_x[0];
		  y0=j0*hy-pLsdft->Brange_y[0];
		  z0=k0*hz-pLsdft->Brange_z[0];
		  
		  x1=i1*hx-pLsdft->Brange_x[0];
		  y1=j1*hy-pLsdft->Brange_y[0];
		  z1=k1*hz-pLsdft->Brange_z[0];
		  
		  X=X-pLsdft->Brange_x[0];
		  Y=Y-pLsdft->Brange_y[0];
		  Z=Z-pLsdft->Brange_z[0];

		  // evaluate the function value using interpolation 
		  arrVeffBoundary[k][j][i]=TrilinearInterpolation(x0,y0,z0,x1,y1,z1,c000,c001,c010,c011,c100,c101,c110,c111,X,Y,Z);
		    }else
		    {
		      // nearest neighbor interpolation (co-incident point)
		      arrVeffBoundary[k][j][i]=c000;
		    }
		  
		  //
		  */
		  // till here
		}
	      else{
		arrVeffBoundary[k][j][i]=0.0;
	      }
	    }
	}
    }
  
  DMDAVecRestoreArray(pLsdft->daloccell,pLsdft->VeffNodalBoundary,&arrVeffBoundary);
  VecRestoreArray(pLsdft->BVeff[0],&arrBVeff);
  //printf("printing VeffNodalBoundary \n");

  //  VecView(pLsdft->VeffNodalBoundary,PETSC_VIEWER_STDOUT_SELF);


  return;
}

void Create_SupercellAtomList(LSDFT_OBJ *pLsdft)
{

  VecDuplicate(pLsdft->Atompos,&pLsdft->AtomposSC);
  VecCopy(pLsdft->Atompos,pLsdft->AtomposSC);
  pLsdft->nAtomsSC = pLsdft->nAtoms;
  pLsdft->startPosSC[0]=0;
  pLsdft->endPosSC[0]=3*pLsdft->nAtoms;

  
}

////////////////////////////////////////////////////////////////////////////////////////////////
//Create_VeffNodalBoundary: boundary Veff values                                              //
////////////////////////////////////////////////////////////////////////////////////////////////
void Create_SupercellAtomList_OLD(LSDFT_OBJ *pLsdft)
{
  // first go over the region Rc+rc(=4) from the edges of the computational domain

  int Ix,Iy,Iz;
  PetscScalar *pAtompos,*pAtomposSC,*pBAtompos;
  PetscScalar rc=4.0;
  PetscScalar Xlim=pLsdft->range_x+pLsdft->Rcut+rc;
  PetscScalar Ylim=pLsdft->range_y+pLsdft->Rcut+rc;
  PetscScalar Zlim=pLsdft->range_z+pLsdft->Rcut+rc;

  PetscScalar Shiftx,Shifty,Shiftz;
   

  // PetscScalar Xlim=pLsdft->range_x;
  //PetscScalar Ylim=pLsdft->range_y;
  //PetscScalar Zlim=pLsdft->range_z;

  // go over all types of atoms
  // first find the maximum of Ntype and BNtype[0] (later iterate on all 6 sides)
  // int NtypeMax=pLsdft->BNtype[0];
  // if(pLsdft->BNtype[0]<pLsdft->Ntype)
  //   {
  // 	 maxNtype=pLsdft->Ntype;
  //   }

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);



  int indexcell;

  PetscScalar Lx,Ly,Lz;
  Lx=2.0*pLsdft->Brange_x[0];
  Ly=2.0*pLsdft->Brange_y[0];
  Lz=2.0*pLsdft->Brange_z[0];
  
  Shiftx=pLsdft->range_x-pLsdft->Brange_x[0];
  Shifty=pLsdft->range_y-pLsdft->Brange_y[0];
  Shiftz=pLsdft->range_z-pLsdft->Brange_z[0];

  // Ix=ceil(2.0*Xlim/Lx);
  // Iy=ceil(2.0*Ylim/Ly);
  // Iz=ceil(2.0*Zlim/Lz);

  Ix=ceil(2.0*(Xlim+Shiftx)/Lx);
  Iy=ceil(2.0*(Ylim+Shifty)/Ly);
  Iz=ceil(2.0*(Zlim+Shiftz)/Lz);


  int index=0,at,count=0,indexSC,start,end,poscnt,PP,QQ,RR;
  PetscScalar x0,y0,z0,X0,Y0,Z0;
  pLsdft->noetotSC = pLsdft->noetot;

  printf("Lx=%lf, Xlim=%lf Ix=%d\n",Lx,Xlim,Ix);
  printf("Ly=%lf, Ylim=%lf Iy=%d\n",Ly,Ylim,Iy);
  printf("Lz=%lf, Zlim=%lf Iz=%d\n",Lz,Zlim,Iz);

  // next go over these atoms and check if they are in the supercell but outside the computational domain
  // do getarray
  
  VecGetArray(pLsdft->BAtompos[0],&pBAtompos);
  for(at=0;at<pLsdft->BNtype[0]; at++)
    {
      start = (int) floor(pLsdft->BstartPos[0][at]/3);
      end = (int) floor(pLsdft->BendPos[0][at]/3);
     
      for(poscnt = start; poscnt <= end; poscnt++)
	{
	  X0 = pBAtompos[index++]; 
	  Y0 = pBAtompos[index++];
	  Z0 = pBAtompos[index++];

	  for(RR=-Iz;RR<=Iz;RR++)
	    {
	      for(QQ=-Iy;QQ<=Iy;QQ++)
		{
		  for(PP=-Ix;PP<=Ix;PP++)
		    {
		      /*
		       * periodic map of atomic position
		       */
		      x0 = X0-Shiftx + 2.0*PP*pLsdft->Brange_x[0];
		      y0 = Y0-Shifty + 2.0*QQ*pLsdft->Brange_y[0];
		      z0 = Z0-Shiftz + 2.0*RR*pLsdft->Brange_z[0];
		      // if(((x0>=-Xlim && x0<-pLsdft->range_x) || (x0<=Xlim && x0>pLsdft->range_x)) || ((y0>=-Ylim && y0<-pLsdft->range_y) || (y0<=Ylim && y0>pLsdft->range_y)) || ((z0>=-Zlim && z0<-pLsdft->range_z) || (z0<=Zlim && z0>pLsdft->range_z)))
		      // x0=round(x0*1000)/1000;                                                                                                                                                 
		      // y0=round(y0*1000)/1000;                                                                                                                                                  
		      // z0=round(z0*1000)/1000;  
		      if(x0>=-Xlim && x0<=Xlim && y0>=-Ylim && y0<=Ylim && z0>=-Zlim && z0<=Zlim)
			{
			  if (x0<-pLsdft->range_x || x0>=pLsdft->range_x || y0<-pLsdft->range_y || y0>=pLsdft->range_y || z0<-pLsdft->range_z || z0>=pLsdft->range_z )
			    {
			      // REMEMBER to change the above condition >= to > later when the computational domain also includes the right face atoms!
			      count++;
			      pLsdft->noetotSC += pLsdft->Bnoe[0][at];
			    }
			}
			
		    }
		}
      
	    }
	}
    }

  //  printf(" pLsdft->noetotSC=%d \n", pLsdft->noetotSC);
  // now allocate memory for storing all atoms in the supercell
  pLsdft->nAtomsSC = count+pLsdft->nAtoms;

  printf("pLsdft->nAtomsSC=%d,count=%d \n", pLsdft->nAtomsSC,count);

  VecCreate(PETSC_COMM_SELF, &pLsdft->AtomposSC);
  VecSetSizes(pLsdft->AtomposSC,PETSC_DECIDE,3*pLsdft->nAtomsSC);
  VecSetFromOptions(pLsdft->AtomposSC);
  //    printf("pLsdft->Ntype=%d\n",pLsdft->Ntype);


  // vecgetarray
  VecGetArray(pLsdft->AtomposSC,&pAtomposSC);
  VecGetArray(pLsdft->Atompos,&pAtompos);
  index=0,indexSC=0,indexcell=0;

  for(at=0;at<pLsdft->Ntype;at++)
    {
      //first copy all atoms of this type which are in the computational domain 
      // NOTE: Ntype>=BNtype for now
      // NOTE: something might be missing here for Ntype neq BNtype
      
      start = (int) floor(pLsdft->startPos[at] / 3);
      end = (int) floor(pLsdft->endPos[at] / 3);
      // store startpos for this type of atom
      //	printf("aaa");
      pLsdft->startPosSC[at]=indexSC;
      //	printf("pLsdft->startPosSC[at]=%d \n",	pLsdft->startPosSC[at]);

      for (poscnt = start; poscnt <= end; poscnt++) 
	{
	  pAtomposSC[indexSC++]=pAtompos[indexcell++]; //x
	  pAtomposSC[indexSC++]=pAtompos[indexcell++]; //y
	  pAtomposSC[indexSC++]=pAtompos[indexcell++]; //z
	}

      //	 printf("indexSC=%d \n",indexSC);

      if(at<pLsdft->BNtype[0])
	{
	  start = (int) floor(pLsdft->BstartPos[0][at]/3);
	  end = (int) floor(pLsdft->BendPos[0][at]/3);
      
	  for(poscnt = start; poscnt <= end; poscnt++)
	    {
	      X0 = pBAtompos[index++]; // to change
	      Y0 = pBAtompos[index++];
	      Z0 = pBAtompos[index++];

	      for(RR=-Iz;RR<=Iz;RR++)
		{
		  for(QQ=-Iy;QQ<=Iy;QQ++)
		    {
		      for(PP=-Ix;PP<=Ix;PP++)
			{
			  /*
			   * periodic map of atomic position
			   */
			  x0 = X0-Shiftx + PP*Lx;
			  y0 = Y0-Shifty + QQ*Ly;
			  z0 = Z0-Shiftz + RR*Lz;
			  // if(((x0>=-Xlim && x0<-pLsdft->range_x) || (x0<=Xlim && x0>pLsdft->range_x)) || ((y0>=-Ylim && y0<-pLsdft->range_y) || (y0<=Ylim && y0>pLsdft->range_y)) || ((z0>=-Zlim && z0<-pLsdft->range_z) || (z0<=Zlim && z0>pLsdft->range_z)))

			  // x0=round(x0*1000)/1000;                                                                                                                                                   
			  // y0=round(y0*1000)/1000;                                                                                                                                                   
			  // z0=round(z0*1000)/1000; 

			  if(x0>=-Xlim && x0<=Xlim && y0>=-Ylim && y0<=Ylim && z0>=-Zlim && z0<=Zlim)
			    {
			      if (x0<-pLsdft->range_x || x0>=pLsdft->range_x || y0<-pLsdft->range_y || y0>=pLsdft->range_y || z0<-pLsdft->range_z || z0>=pLsdft->range_z )
				{
				  // REMEMBER to change the above condition >= to > later when the computational domain also includes the right face atoms!
				  pAtomposSC[indexSC++]=x0;
				  pAtomposSC[indexSC++]=y0;
				  pAtomposSC[indexSC++]=z0;
			    
				}
			    }
			
			}
		    }
      
		}
	    }
	}
	  pLsdft->endPosSC[at]=indexSC-1;
	  //    printf("pLsdft->endPosSC[at]=%d \n",pLsdft->endPosSC[at]);
	  pLsdft->noaSC[at]=floor(pLsdft->endPosSC[at]/3)-floor(pLsdft->startPosSC[at]/3)+1;
	  pLsdft->noeSC[at]= pLsdft->noaSC[at]*pLsdft->noe[at];
	  printf("pLsdft->noaSC[at]=%d, pLsdft->noeSC[at]=%d,  pLsdft->endPosSC[at]=%d",pLsdft->noaSC[at],pLsdft->noeSC[at],pLsdft->endPosSC[at]);
    }


  VecRestoreArray(pLsdft->AtomposSC,&pAtomposSC);
  VecRestoreArray(pLsdft->Atompos,&pAtompos);
  VecRestoreArray(pLsdft->BAtompos[0],&pBAtompos);

  if(rank==0)
    {
      printf("SClist \n");
      VecView(pLsdft->AtomposSC,PETSC_VIEWER_STDOUT_SELF);

    }


  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//                    Calculate_IntWts: Calculates and stores the integration weights        //
///////////////////////////////////////////////////////////////////////////////////////////////
void Calculate_IntWts(LSDFT_OBJ* pLsdft)
{

  PetscInt i,j,k,l,colidx,gxdim,gydim,gzdim,xcor,ycor,zcor,lxdim,lydim,lzdim;
  PetscScalar r,dr,dV,***WtArray;
 
  DMDAGetCorners(pLsdft->da,&xcor,&ycor,&zcor,&lxdim,&lydim,&lzdim);
  dV = pLsdft->delVol;

  DMDAVecGetArray(pLsdft->da,pLsdft->IntWts,&WtArray);
   for(k=zcor; k<zcor+lzdim; k++)
    for(j=ycor; j<ycor+lydim; j++)
      for(i=xcor; i<xcor+lxdim; i++)
	{
	  WtArray[k][j][i] = 1.0;
	  
	  if(i==0)
	    WtArray[k][j][i] = 0.5*WtArray[k][j][i];
	  if(i==pLsdft->numPoints_x-1)
	    WtArray[k][j][i] = 0.5*WtArray[k][j][i];
	   if(j==0)
	     WtArray[k][j][i] = 0.5*WtArray[k][j][i];
	  if(j==pLsdft->numPoints_y-1)
	     WtArray[k][j][i] = 0.5*WtArray[k][j][i];
	  if(k==0)
	     WtArray[k][j][i] = 0.5*WtArray[k][j][i];
	  if(k==pLsdft->numPoints_z-1)
	     WtArray[k][j][i] = 0.5*WtArray[k][j][i];
	  
	}
  DMDAVecRestoreArray(pLsdft->da,pLsdft->IntWts,&WtArray);

}

