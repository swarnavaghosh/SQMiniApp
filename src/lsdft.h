/*
  |
  | file name: sddft.h
  |
  | Description: This file contains the function declarations
  |
  | Authors: Swarnava Ghosh
  |
  |
  |-------------------------------------------------------------------------------------------*/
#ifndef _LSDFT_
#define _LSDFT_
#include "petsc.h"
#include "petscksp.h"
#include "ilsdft.h"
#include "assert.h"
#undef _DEBUG
//#define KSP_TYPE KSPCG
#define KSP_TYPE KSPGMRES
#define MAX_TABLE_SIZE 11000
#define ITMAXBRENTS 100
#define EPSILON 1e-16

/*
 * initObjs.cc
 */
void Get_Input(LSDFT_OBJ *pLsdft);
void Initialize(LSDFT_OBJ *pLsdft);
PetscErrorCode Create_Objects(LSDFT_OBJ *pLsdft);
PetscErrorCode Create_Laplacian(LSDFT_OBJ *pLsdft);
PetscErrorCode Create_LaplacianLocal(LSDFT_OBJ *pLsdft);
PetscErrorCode Create_Gradient(LSDFT_OBJ *pLsdft);
void Calculate_SplineFit(LSDFT_OBJ *pLsdft);
PetscErrorCode Create_Orbitals(LSDFT_OBJ *pLsdft);
//PetscErrorCode Create_Orbitals_Complex(LSDFT_OBJ* pLsdft);
void Calculate_PseudochargeCutoff(LSDFT_OBJ *pLsdft);
void Calculate_NonlocalIndex(LSDFT_OBJ *pLsdft);
PetscScalar Calculate_NonlocalCoefficients(LSDFT_OBJ *pLsdft);
PetscScalar CalculateNonlocalProjectors(LSDFT_OBJ *pLsdft);
PetscScalar CalculateNonlocalProjectorsForces(LSDFT_OBJ *pLsdft);
//PetscScalar CalculateNonlocalProjectors_Complex(LSDFT_OBJ* pLsdft);
void Objects_Destroy(LSDFT_OBJ *pOfdft);
void Set_VecZero(LSDFT_OBJ *pLsdft);
void GetInfluencingAtoms_Local(LSDFT_OBJ *pLsdft);
void GetInfluencingAtoms_Nonlocal(LSDFT_OBJ *pLsdft);
void GetInfluencingAtoms_LocalForces(LSDFT_OBJ *pLsdft);
void GetInfluencingAtoms_NonlocalForces(LSDFT_OBJ *pLsdft);
int Check_Overlap(int xcor, int ycor, int zcor, int lxdim, int lydim, int lzdim, int xs, int ys, int zs, int xl, int yl,
    int zl, int *xstart, int *ystart, int *zstart, int *xend, int *yend, int *zend);
//PetscScalar kpointWeight(PetscScalar kx,PetscScalar ky,PetscScalar kz);
//PetscScalar Calculate_kpoints(LSDFT_OBJ* pLsdft);
//PetscErrorCode Create_subcommunicator(LSDFT_OBJ* pLsdft);
int InitializeHamiltonian_Serial(LSDFT_OBJ *pLsdft);
PetscErrorCode Create_VeffCommunicaor(LSDFT_OBJ *pLsdft);
void Create_VeffNodalBoundary(LSDFT_OBJ *pLsdft);
void Create_SupercellAtomList(LSDFT_OBJ *pLsdft);
void Calculate_IntWts(LSDFT_OBJ* pLsdft);


/*
 * relaxatoms.cc
 */
PetscErrorCode Calculate_Electrondensity(LSDFT_OBJ *pLsdft);
void Get_AtomicForces(LSDFT_OBJ *pLsdft);
void NLCG(LSDFT_OBJ *pLsdft);
void Display_Atompos(LSDFT_OBJ *pLsdft);
void Display_Relax(LSDFT_OBJ *pLsdft);
PetscErrorCode Final_Output(LSDFT_OBJ *pLsdft);
void Periodic_MapAtoms(LSDFT_OBJ *pLsdft);
void NVE_MD(LSDFT_OBJ *pLsdft);
void Print_ElecDens(LSDFT_OBJ *pLsdft);
void Read_ElecDens(LSDFT_OBJ *pLsdft);
void Print_MinMaxLambda(LSDFT_OBJ *pLsdft);
void Read_MinMaxLambda(LSDFT_OBJ *pLsdft);


/*
 * readfiles.cc
 */
void Read_Parameters(LSDFT_OBJ *pLsdft);
void Read_Ion(LSDFT_OBJ *pLsdft);
void Read_Relax(LSDFT_OBJ *pLsdft);
void Read_Pseudopotential(LSDFT_OBJ *pLsdft);
void TransferFields(LSDFT_OBJ *pLsdft);
void DisplayPotential(LSDFT_OBJ *pLsdft);

/*
 * poisson.cc
 */
PetscErrorCode SolvePoisson(LSDFT_OBJ *pLsdft);

/*
 * occupation.cc
 */
PetscScalar constraint_lssgq(LSDFT_OBJ *pLsdft, PetscScalar lambdaf);
PetscScalar FermiDirac(PetscScalar bet, PetscScalar lambda, PetscScalar lambdaf);
PetscErrorCode CalculateDensity_lssgq(LSDFT_OBJ *pLsdft);
//PetscErrorCode CalculateDensity_Complex(LSDFT_OBJ* pLsdft,Mat* U1,Mat* U2,int k);
PetscScalar findRootBrent(LSDFT_OBJ *pLsdft, PetscScalar x1, PetscScalar x2, PetscScalar tol);
//PetscErrorCode RotatePsi(LSDFT_OBJ* pLsdft, Mat* Psi, Mat* Q, Mat* PsiQ);
//PetscErrorCode RotatePsi_Complex(LSDFT_OBJ* pLsdft, Mat* U1,Mat* U2, Mat* Q1,Mat* Q2, Mat* UQ1, Mat* UQ2);
//void Calculate_nonlocalindex(LSDFT_OBJ* pLsdft);
//void Display_EigenValues(LSDFT_OBJ* pLsdft);
//PetscErrorCode SumDensity_Groups(LSDFT_OBJ* pLsdft);
//void Scale_Density(LSDFT_OBJ* pLsdft);

/*
 * Chebyshev filtering
 */
//void ChebyshevFiltering(LSDFT_OBJ* pLsdft,int m,PetscScalar a,PetscScalar b,PetscScalar a0,int k);
//PetscErrorCode Mult_HamiltonianOrbital(LSDFT_OBJ* pLsdft,Mat *Psi1,Mat *Psi2);
PetscErrorCode Mult_HamiltonianVector(LSDFT_OBJ *pLsdft, Vec *V1, Vec *V2);
//PetscErrorCode Mult_HamiltonianOrbital_Complex(LSDFT_OBJ* pLsdft,PetscScalar k1,PetscScalar k2,PetscScalar k3,Mat *Psi1,Mat *Psi1_imag,Mat *Psi2,Mat *Psi2_imag);
//PetscErrorCode Mult_HamiltonianVector_Complex(LSDFT_OBJ* pLsdft,PetscScalar k1,PetscScalar k2,PetscScalar k3,Vec *V1,Vec *V1_imag,Vec *V2,Vec *V2_imag);
PetscErrorCode Mult_HamiltonianVector(LSDFT_OBJ *pLsdft, Vec *V1, Vec *V2, int k,int j, int i);


/*
 * scf.cc
 */
PetscScalar SCF_LSSGQ(LSDFT_OBJ *pLsdft);
//PetscErrorCode ProjectMatrices(LSDFT_OBJ* pLsdft, Mat* Psi,Mat *Hsub,Mat* Msub);
//PetscErrorCode ProjectMatrices_Complex(LSDFT_OBJ* pLsdft, Mat* U1,Mat* U2,Mat* Hsub1,Mat* Hsub2,Mat* Msub1,Mat* Msub2,int k);
//PetscErrorCode SolveGeneralizedEigen_Real(LSDFT_OBJ* pLsdft,Mat* Hsub, Mat* Msub);
//PetscErrorCode SolveGeneralizedEigen_Complex(LSDFT_OBJ* pLsdft,Mat* Hsub1, Mat* Hsub2, Mat* Msub1, Mat* Msub2,int k);
//void  SetZero_Eigen(LSDFT_OBJ* pLsdft);

PetscScalar CalculateSpectralNodesAndWeights(LSDFT_OBJ *pLsdft, int p, int LIp,int k,int j, int i);
PetscScalar NodesAndWeights_LSSGQ(LSDFT_OBJ *pLsdft);
int VeffGlobalToLocal(LSDFT_OBJ *pLsdft);
int InitializeHamiltonian_Serial(LSDFT_OBJ *pLsdft);
int ReInitializeHamiltonian_Serial(LSDFT_OBJ *pLsdft);
PetscScalar CalculateSpectralNodesAndWeightsTwoPoint(LSDFT_OBJ *pLsdft, int p, int q, int LIp, int flag);
int VeffTransfer(LSDFT_OBJ *pLsdft);
int FormVeffNodal(LSDFT_OBJ *pLsdft,int K,int J,int I);

/*
 * tools.cc
 */
PetscScalar fract(PetscInt n, PetscInt k);
PetscErrorCode Mult_ComplexScalar(PetscScalar A1, PetscScalar A2, PetscScalar B1, PetscScalar B2, PetscScalar *C1,
    PetscScalar *C2);
PetscErrorCode Exp_ComplexScalar(PetscScalar k1, PetscScalar k2, PetscScalar k3, PetscScalar v1, PetscScalar v2,
    PetscScalar v3, PetscScalar *C1, PetscScalar *C2);
PetscErrorCode Mult_ComplexMatMat(LSDFT_OBJ *pLsdft, Mat *A1, Mat *A2, Mat *B1, Mat *B2, Mat *C1, Mat *C2);
PetscErrorCode Mult_ComplexMatVec(LSDFT_OBJ *pLsdft, Mat *A1, Mat *A2, Vec *V1, Vec *V2, Vec *C1, Vec *C2);
PetscErrorCode VecNorm_Complex(Vec *Vr, Vec *Vi, PetscScalar *norm);
PetscErrorCode VecDot_Complex(Vec *V1r, Vec *V1i, Vec *V2r, Vec *V2i, PetscScalar *rvecdot, PetscScalar *ivecdot);
PetscErrorCode VecScale_Complex(Vec *Vr, Vec *Vi, PetscScalar scale);
void SplineEvaluate(PetscScalar *X1, PetscScalar *Y1, int len1, PetscScalar *X2, PetscScalar *Y2, int len2,
    PetscScalar *YD);
PetscScalar RealSphericalHarmonic(PetscScalar x, PetscScalar y, PetscScalar z, int l, int m);
PetscScalar ComplexSphericalHarmonic(PetscScalar x, PetscScalar y, PetscScalar z, int l, int m, PetscScalar *SHreal,
    PetscScalar *SHimag);
void getYD_gen(PetscScalar *X, PetscScalar *Y, PetscScalar *YD, int len);
void tridiag_gen(PetscScalar *A, PetscScalar *B, PetscScalar *C, PetscScalar *D, int len);
void Lanczos(LSDFT_OBJ *pLsdft, PetscScalar *EigenMin, PetscScalar *EigenMax, int kpt);
void TridiagEigenSolve(PetscScalar diag[], PetscScalar subdiag[], int n, PetscScalar *EigenMin,
    PetscScalar *EigenMax);
void TridiagEigenVecSolve_NodesAndWeights(LSDFT_OBJ *pLsdft, PetscScalar diag[], PetscScalar subdiag[], int n,
    int LIp);


PetscScalar ScalarNorm_Complex(PetscScalar ar, PetscScalar ai);
//void kPointLanczos(LSDFT_OBJ* pLsdft,PetscScalar* EigenMin,PetscScalar* EigenMax,int kpt);

PetscScalar TrilinearInterpolation(PetscScalar x0,PetscScalar y0,PetscScalar z0,PetscScalar x1,PetscScalar y1,PetscScalar z1,PetscScalar c000,PetscScalar c001,PetscScalar c010,PetscScalar c011,PetscScalar c100,PetscScalar c101,PetscScalar c110,PetscScalar c111,PetscScalar x,PetscScalar y,PetscScalar z);
void FindMeshCorners(PetscScalar x, PetscScalar y,PetscScalar z,PetscScalar Rx,PetscScalar Ry,PetscScalar Rz, PetscScalar rx,PetscScalar ry,PetscScalar rz,int nx,int ny,int nz,int *i0,int *j0,int *k0,int *i1,int *j1,int *k1,PetscScalar *X, PetscScalar *Y,PetscScalar *Z);
void EvaluateCorners(PetscScalar *arrVals,int nx,int ny,int nz,int i0,int j0,int k0,int i1,int j1,int k1,PetscScalar *c000,PetscScalar *c001,PetscScalar *c010,PetscScalar *c011,PetscScalar *c100,PetscScalar *c101,PetscScalar *c110,PetscScalar *c111);

/*
 * multipole.cc
 */
PetscErrorCode MultipoleExpansion_Phi(LSDFT_OBJ *pLsdft, Vec *RhopBvec);
PetscErrorCode BulkBoundary_Phi(LSDFT_OBJ *pLsdft, Vec *RhopBvec);

/*
 * mix.cc
 */
PetscScalar mix(LSDFT_OBJ *pLsdft, PetscInt its, Vec fk);
 
/*
 * forces.cc
 */
void Calculate_localforce(LSDFT_OBJ *pLsdft);
PetscErrorCode Calculate_nonlocalforce(LSDFT_OBJ *pLsdft, Mat *Psi);
void Symmetrysize_force(LSDFT_OBJ *pLsdft);
void Display_force(LSDFT_OBJ *pLsdft);
//PetscErrorCode Calculate_nonlocalforce_Complex(LSDFT_OBJ* pLsdft);

/*
 * ExchangeCorrelation.cc
 */
PetscErrorCode Vxc_Calc_CA(LSDFT_OBJ *pLsdft);
PetscErrorCode Exc_Calc_CA(LSDFT_OBJ *pLsdft);
PetscErrorCode Vxc_Calc_CA_PZ(LSDFT_OBJ *pLsdft);
PetscErrorCode Exc_Calc_CA_PZ(LSDFT_OBJ *pLsdft);
PetscErrorCode Vxc_Calc_CA_PW(LSDFT_OBJ *pLsdft);
PetscErrorCode Exc_Calc_CA_PW(LSDFT_OBJ *pLsdft);
PetscErrorCode Vxc_Calc(LSDFT_OBJ *pLsdft);
PetscErrorCode Exc_Calc(LSDFT_OBJ *pLsdft);

/*
 * energy.cc
 */
void SystemEnergy_Calc_lssgq(LSDFT_OBJ *pLsdft);
void Display_Energy(LSDFT_OBJ *pLsdft);
void Display_FinalEnergy(LSDFT_OBJ *pLsdft);
void MD_energies(LSDFT_OBJ *pLsdft);

/*
 * electrostatics.cc
 */
void Generate_pseudocharge(LSDFT_OBJ *pLsdft);
PetscScalar PseudopotReference(PetscScalar r, PetscScalar rcut, PetscScalar Znucl);

/*
 * lsdft_sq.cc
 */
PetscErrorCode Calculate_nonlocalforce(LSDFT_OBJ *pLsdft);
PetscScalar CalculateGradientDensityMatrixColumn(LSDFT_OBJ *pLsdft, int p, int LIp,int Kp, int Jp, int Ip, double Zeta, double Chi);
PetscScalar ChebyshevCoeff(LSDFT_OBJ *pLsdft, double *d, double lambda_fhat, double bethat);
PetscErrorCode CalculateNonlocalForceContribution(LSDFT_OBJ *pLsdft,int Ip,int Jp,int Kp,PetscScalar weights);
PetscScalar LanczosMaxMinEigen(LSDFT_OBJ *pLsdft, int p, int LIp,int Kp, int Jp, int Ip);


/*
 * preconditioner
 */
int Create_HelmholtzOperator(LSDFT_OBJ *pLsdft);
int Initialize_HelmholtzSolver(LSDFT_OBJ *pLsdft);
int Apply_Preconditioner(LSDFT_OBJ *pLsdft,int t);

#endif
