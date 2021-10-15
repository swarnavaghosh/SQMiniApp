/* 
  |
  | file name: isddft.h          
  |
  | Description: This file contains the variables required
  |
  | Authors: Swarnava Ghosh
  |
  |    
  |-------------------------------------------------------------------------------------------*/
#ifndef _IODFD_
#define _IODFD_
#include "petsc.h"
#include "petscdmda.h"
#define MAX_ORDER 10
#define EPSNEV 1
/*
 * structure storing the atomic positions of the overlap atoms
 */
typedef struct {
    /*
     * this structure is created for every atom type.
     * the size of the pointer is the number of overlapping atoms
     * two copies of this are created. One for rb and other for rc
     */
    int Natoms;                 // number of influence atoms in simulation domain + image atoms 
    //  int Natoms_domain; // number of influence atoms in simulation domain
    double *X0;                 // x coordinates of simulation domain atoms
    double *Y0;                 // y coordinates of simulation domain atoms
    double *Z0;                 // z coordinates of simulation domain atoms
    //double  ******Chi; // nonlocal projectors, Chi_{Jlm}(i,j,k)
    int *xstart;                // starting node number in x direction of overlap region
    int *ystart;                // starting node number in y direction of overlap region
    int *zstart;                // starting node number in z direction of overlap region
    int *xend;                  // ending node number in x direction of overlap region
    int *yend;                  // ending node number in y direction of overlap region
    int *zend;                  // ending node number in z direction of overlap region

    int *xindex;
    int *yindex;
    int *zindex;

} OVERLAP_OBJ;

typedef struct {
    double ******Chi;           // nonlocal projectors, Chi_{Jlm}(i,j,k) 
} PROJECTOR_OBJ;

/*
 * structure storing the pseudopotential information 
 */
typedef struct {
    PetscScalar *Vloc;          ///< stores local part of pseudopotential 
    PetscScalar **V;            ///< stores the component of pseudopotential
    PetscScalar **U;            ///< stores the component of pseudowavefunction
    PetscScalar **UDeltaV;
    PetscScalar **SplineFitUDeltaV;
    PetscScalar **SplineFitU;
    PetscScalar *SplineFitVloc;
    PetscScalar *SplineFitIsoAtomDen;

    PetscScalar *uu;            ///< stores isolated atom electron density
    PetscScalar *RadialGrid;    ///< stores the radial grid
    PetscScalar *rc;            ///< component pseudopotential cutoff
    PetscInt lmax;              ///< maximum pseudopotential component
    PetscInt size;              ///< size of the arrays storing the pseudopotentials   
    PetscScalar *Gamma;

} PSD_OBJ;

/*
 * structure storing the variables required by the functions of SPARC 
 */
typedef struct {
    PetscInt numPoints_x;       ///< stores the number of nodes in x-direction
    PetscInt numPoints_y;       ///< stores the number of nodes in y-direction
    PetscInt numPoints_z;       ///< stores the number of nodes in z-direction

    PetscInt SCnumPoints_x;     ///< stores the number of supercell nodes in x-direction
    PetscInt SCnumPoints_y;     ///< stores the number of supercell nodes in y-direction
    PetscInt SCnumPoints_z;     ///< stores the number of supercell nodes in z-direction

    PetscInt nx_loc;
    PetscInt ny_loc;
    PetscInt nz_loc;

    AO aodmda1;

    PetscInt order;             ///< stores half of finite difference order
    PetscInt nAtoms;            ///< stores the total number of atoms in the simulation
    PetscInt nAtoms_UnitCell;   ///< total number of atoms in a unit cell for periodic calculations
  PetscInt nAtomsSC;
    PetscInt *noa, *noe;        ///< For each type of atom, noa stores the number of atoms of that type, noe stores the number of valence electrons of that type
    PetscInt noetot;            ///< stores the total number of valence electrons in the simulation
  PetscInt noetotSC;
    PetscInt *noaSC, *noeSC; 
    PetscInt ChebDegree; ///< stores the degree of Chebychev polynomial used for the Chebychev filter in predictor and boundary condition calculation
    PetscInt N_qp;              // number of quadrature  points for gauss quadrature
    PetscInt NplCC;             // degree of polynomial for Clenshaw Curtis quadrature
    PetscInt MAXITSCF;          ///< Maximum number of SCF iterations  
    PetscInt MAXIT_NLCG;        ///< Maximum number of NLCG iterations
    PetscInt T_indexlength;
    PetscInt T_length;
    PetscInt *T_index;
    PetscInt IP_length;
    PetscInt *IP_index;
    PetscInt *ProcDomains;
    PetscInt *sendcounts;
    PetscInt *sdispls;
    PetscInt *recvcounts;
    PetscInt *rdispls;
    PetscInt *neigh;
    PetscInt *OwnershipRange;
    PetscInt nneigh;
    
  PetscInt Bnoetot[6]; //TODO allocate
  PetscInt BnAtoms[6]; // TODO allocate
  PetscInt *Bnoa[6], *Bnoe[6]; // TODO allocate
  PetscScalar Brange_x[6];        ///< stores half of domain length in x direction
    PetscScalar Brange_y[6];        ///< stores half of domain length in y direction
    PetscScalar Brange_z[6];        ///< stores half of domain length in z direction
   PetscInt BNx[6];
   PetscInt BNy[6];
   PetscInt BNz[6];
  

    VecScatter VeffScatterCtx;

    PetscScalar *Veff_neigh; // effective potential from neighboring processors
  //PetscScalar *VeffNodal;  // effective potential around a node
    PetscScalar range_x;        ///< stores half of domain length in x direction
    PetscScalar range_y;        ///< stores half of domain length in y direction
    PetscScalar range_z;        ///< stores half of domain length in z direction
    PetscScalar delta_x;        ///< mesh spacing
    PetscScalar delta_y;        ///< mesh spacing
    PetscScalar delta_z;        ///< mesh spacing
    PetscScalar elecN;          ///< stores integral of total pseudocharge density 
    PetscScalar Eatom;          ///< stores total energy per atom  
    PetscScalar Eself;          ///< stores self energy of the nuclei
    PetscScalar Eself_TM;       ///< stores self energy of the nuclei calculated using reference pseudopotential
    PetscScalar Exc;            ///< stores exchange correlation energy
    PetscScalar Ecorrection;    ///< stores electrostatic correction energy
    PetscScalar Entropy;        ///< stores the electronic entropy
    PetscScalar Eband;          ///< stores the band structure energy
    PetscScalar TOLSCF;         ///< stores the scf tolerence
    PetscScalar KSPTOL;         ///< stores the tolerence for solving the Poisson equation
    PetscScalar NLCGTOL;        ///< stores the tolerance for the non-linear conjugate gradient method for atomic relaxation
    PetscScalar LANCZOSTOL;     ///< stores the tolerence for the lanczos iteration
    PetscScalar PseudochargeRadiusTOL;  ///< stores the tolerance for calculation of pseudocharge density radius
    PetscScalar REFERENCE_CUTOFF;   ///< stores the cutoff for the reference pseudopotential
    PetscScalar ConstShift;     /// stores the shift for correction electrostatic potential under periodic boundary conditions
    PetscScalar MaxForce;       /// maximum component of forces
    PetscScalar NormForce;      /// average l2 norm of forces
    PetscScalar delVol;
    PetscScalar HalfdelVol;
    PetscScalar OrbScalingFactor;
    PetscScalar DMx;
    PetscScalar DMy;
    PetscScalar DMz;
    PetscScalar kB;
    PetscScalar Ceh;
    PetscScalar Rcut;

    // lssgq weights and nodes
    PetscScalar **lssgq_weights;
    PetscScalar **lssgq_lambda;
    PetscScalar *lambda_max;
    PetscScalar *lambda_min;
  PetscScalar *ChebCoeff;

    PetscScalar coeffs_x[MAX_ORDER + 1];    ///< stores the weights of the finite-difference laplacian
    PetscScalar coeffs_y[MAX_ORDER + 1];    ///< stores the weights of the finite-difference laplacian
    PetscScalar coeffs_z[MAX_ORDER + 1];    ///< stores the weights of the finite-difference laplacian
    PetscScalar stencil_coeffs_x[MAX_ORDER + 1];    ///< stores the weights of the finite-difference laplacian
    PetscScalar stencil_coeffs_y[MAX_ORDER + 1];    ///< stores the weights of the finite-difference laplacian
    PetscScalar stencil_coeffs_z[MAX_ORDER + 1];    ///< stores the weights of the finite-difference laplacian
    PetscScalar coeffs_grad_x[MAX_ORDER + 1];   ///< stores the weights of the finite difference gradient
    PetscScalar coeffs_grad_y[MAX_ORDER + 1];   ///< stores the weights of the finite difference gradient
    PetscScalar coeffs_grad_z[MAX_ORDER + 1];   ///< stores the weights of the finite difference gradient
    PetscScalar Beta;           ///< stores the electronic smearing (1/(k_B * T)) (units: 1/Ha)
    PetscScalar Eelec;
    PetscInt NstatesBC;
    PetscInt Nstates;           ///< stores total number of electronic states used in the simulation  
    PetscInt Ntype;             ///< stores total number of types of atoms used for the simulation
  PetscInt BNtype[6];
    PetscInt SCFNewRhoCalcCtr;  ///< stores the number of times Chebychev filtering is done in the first SCF iteration before updation of electron density

    int RelaxCount;             ///< stores the number of relaxation steps done. This is incremented for every relaxation step.
    PetscScalar MixingParameter;    ///< stores mixing parameter for Anderson Mixing
    int MixingHistory;          ///< stores mixing history for Anderson Mixing 
    PetscInt BC;                ///< Boundary Condition
    MatNullSpace nullspace;     ///< nullspace for poisson solve under periodic boundary conditions
    int PulayFrequency;
    int PulayRestartFlag;
    // int NkptsGroup;
    //int kptParalFlag;
    int PrintForceFlag;
    int PrintAtomPosFlag;
    int PrintElecDensFlag;
    int PrintEigenFlag;
    int ReadElecDensFlag;
    int DipoleMomentFlag;
  int LinInterpEFieldsFlag;
  
    //int NprocGroup;
    //int kptGroupNumber;
    //int Nkpts_symGroup;
    //int Nkpts_start;
    //int Nkpts_end;

    // MPI_Comm kpoint_group_comm;  // dont forget to free communicator
    //MPI_Comm kpoint_intergroup_comm;  // dont forget to free communicator
  MPI_Comm lssgq_comm_effpot;

    PetscScalar lambda_f;       ///< stores the Fermi energy
    PetscScalar *lambda;        ///< stores the eigenvalues in the subspace eigen problem. Size of the array is Nstates
    PetscScalar *CUTOFF_x;      ///< stores the pseudocharge cutoff \f$ r_J^b \f$. Size of array is nAtoms. 
    PetscScalar *CUTOFF_y;      ///< stores the pseudocharge cutoff \f$ r_J^b \f$. Size of array is nAtoms. 
    PetscScalar *CUTOFF_z;      ///< stores the pseudocharge cutoff \f$ r_J^b \f$. Size of array is nAtoms. 
    PetscInt *startPos;         ///< stores the starting position of a particular type of atom in the list of atom positions (stored in Vec Atompos)
    PetscInt *endPos;           ///< stores the ending position of a particular type of atom in the list of atom positions (stored in Vec Atompos)
     PetscInt *startPosSC;         ///< stores the starting position of a particular type of atom in the list of atom positions (stored in Vec Atompos)
    PetscInt *endPosSC;           ///< stores the ending position of a particular type of atom in the list of atom positions (stored in Vec Atompos)
    char **atomType;            ///< stores the string containing the atom type name for every type of atom. example H for Hydrogen, He for Helium, Li for Lithium, etc. This is read from the .ion file 
    char **psdName;             ///< stores the string containing the pseudopotential name for every type of atom. This is read from the .inpt file    
    int *localPsd;              ///< stores the respective local component of pseudopotential for each type of atom. 0 for s, 1 for p, 2 for d, 3 for f.
    char file[50];              ///< stores the input filename 
    char XC[30];                ///< stores the exchange correlation name
    int RelaxFlag;              ///< flag for relaxation of atoms. 1=relaxation, any other number=NO relaxation
    //int ChebyshevCallCounter; ///< variable for storing the number of times chebyshev filter has been called 

  
  char *BElectronicFieldFilename[6];
   int *BlocalPsd[6];  
  int BReadBoundaryEfield[6];
  char **BpsdName[6];   
  char **BatomType[6]; // TODO Allocate
  PetscInt *BstartPos[6];  // TODO allocate
  PetscInt *BendPos[6]; // TODO allocate 
  Vec BAtompos[6]; // TODO allocate 
  Vec AtomposSC; // all atom positions including the ones in supercell

    //PetscInt *nnzDArray;  ///< stores the number of nonzeros in the diagonal block of the nonlocal matrix
    //PetscInt *nnzODArray; ///< stores the number of nonzeros in the off diagonal block of the nonlocal matrix

    DM da;                      ///< DMDA for finite difference laplacian. see PETSc DMDA for more details
    DM daloc;                   //local DMDA
    DM daloccell;               //local DMDA for supercell
  //DM daVeff;
    // DM da_grad; ///< DMDA for finite difference laplacian. see PETSc DMDA for more details  
  
  Vec BelecDensRho[6]; // boundary electron density
  Vec BpotentialPhi[6]; // boundary electrostatic potential
  Vec BVeff[6]; // boundary effective potential


    Vec elecDensRho;            ///< PETSc vector for storing the electron density
    Vec SuperposAtRho;          ///< PETSc vector for storing the guess electron density as obtained from superposition of atoms 
    Vec chrgDensB;              ///< PETSc vector for storing the pseudocharge density. This is calculated from the pseudopotential files given as user input
    Vec chrgDensB_TM;           ///< PETSc vector for storing the reference pseudocharge density. This is used for calculation of energy and force corrections
    Vec potentialPhi;           ///< PETSc vector for storing the electrostatic potential
    Vec twopiRhoPB;             ///< PETSc vector for storing \f$ 2*\pi*(\rho+b) \f$ (the right hand side of the poisson equation solves for the elctrostatic potential)
    Vec Veff;                   ///< PETSc vector for storing effective potential \f$ V_{eff} = \phi + V_{xc} \f$

    Vec Veff_loc;               // local vector storing Veff
    Vec Vk;                     // local vector for lanczos
    Vec Vkm1;                   // local vector for lanczos
    Vec Vkp1;                   // local vector for lanczos
  Vec VeffVktemp;               // temporary vector
    Vec VeffNodal;
  Vec VeffNodalBoundary; // nodal veff if the nodes are outside
  Vec DMcolgrad_x;
  Vec DMcolgrad_y;
  Vec DMcolgrad_z;
  Vec Tkp1;
  Vec Tkm1;
  Vec Tk;
  Vec DMcol;
  Vec IntWts;

    Vec Vxc;                    ///< PETSc vector for storing exchange-correlation potential \f$ V_{xc} \f$
    Vec bjVj;                   ///< PETSc vector for storing the nodal contribution of the repulsive energy
    Vec bjVj_TM;                ///< PETSc vector for storing the nodal contribution of the reference repulsive energy
    Vec Phi_c;                  ///< PETSc vector for storing the correction electrostatic potential
    Vec twopiBTMmBPS;           ///< PETSc vector for storing the difference between the reference pseudocharge density and the pseudocharge density
    Vec PoissonRHSAdd;          ///< PETSc vector for storing the correction term in multipole expansion
    PetscScalar *pForces_corr;  ///< stores the correction in forces on atoms
    Vec forces;                 ///< PETSc vector storing the forces on the atoms. Size of this vector is 3*nAtoms. The first, second and third entries are the x,y and z coordinates of the atomic forces on first atom, the fourth, fifth and sixth entries are the x,y,z coordinates on the atomic positions of second atom and so on.
   Vec nonlocalforces;
    Vec Atompos;                ///< PETSc vector storing the atomic positions. Size of this vector is 3*nAtoms. The first, second and third entries are the x,y and z coordinates of the atomic positions of first atom, the fourth, fifth and sixth entries are the x,y,z coordinates of the atomic positions of second atom and so on.
    Vec mvAtmConstraint;        ///< PETSc vector storing the constraints for movement of atoms. Size of this vector is 3*nAtoms. An entry of 1 means the atom is "movable", hence "net" force on it is non-zero, an entry of 0 means the atom is "fixed" and the "net" force on it is zero. 
    
    Mat laplaceOpr;             ///< PETSc distributed sparse matrix (row major format) for storing -1/2 laplacian. This is a PETSc DMDA type matrix. See PETSc manuals for details on DMDA matrices
    // Mat laplaceOprloc; // local laplacian
  // Mat HamiltonianOpr;         ///< PETSc distributed sparse matrix (row major format) for storing the Hamiltonian operator
    Mat LaplacianOprloc;      // local Hamiltonian operator
  Mat LapPlusVeffOprloc;
  //  Mat HamiltonianOpr_imag;    ///< PETSc distributed sparse matrix (row major format) for storing the imaginary part of Hamiltonian operator
    Mat gradient_x;             ///< PETSc distributed sparse matrix (row major format) for storing x component of the gradient operator. This is a PETSc DMDA type matrix. See PETSc manuals for details on DMDA matrices
    Mat gradient_y;             ///< PETSc distributed sparse matrix (row major format) for storing y component of the gradient operator. This is a PETSc DMDA type matrix. See PETSc manuals for details on DMDA matrices
    Mat gradient_z;             ///< PETSc distributed sparse matrix (row major format) for storing z component of the gradient operator. This is a PETSc DMDA type matrix. See PETSc manuals for details on DMDA matrices

    /* 
     * matrices for storing orbitals
     */
    Mat XOrb;                   ///< PETSc distributed dense matrix (row major format) for storing the orbitals. Stores input to the chebychev filter
    Mat YOrb;                   ///< PETSc distributed dense matrix (row major format) for storing the orbitals. Stores the filtered orbitals
    Mat YOrbNew;                ///< PETSc distributed dense matrix (row major format) for storing the orbitals. Stores a temporary copy of the orbitals in chebychev filtering 
    Mat Hsub;
    Mat Msub;
    Mat Hsub_imag;
    Mat Msub_imag;

    /*
       Mat *XOrb_real; ///< Real part of input to chebyshev filter orbitals for k-point code
       Mat *XOrb_imag; ///< Imaginary part of input to chebyshev filter orbitals for k-point code
       Mat *YOrb_real; ///< Real part of chebyshev filtered orbitals for k-point code
       Mat *YOrb_imag; ///< Imaginary part of chebyshev filtered orbitals for k-point code
       Mat YOrbNew_real; ///< Stores a temporary copy of the real part of orbitals in chebyshev filtering
       Mat YOrbNew_imag; ///< Stores a temporary copy of the imaginary part of orbitals in chebyshev filtering
       Mat ZOrb1; ///< Stores a temporary copy of the real part of orbitals in chebyshev filtering 
     */

    //int Nkpts; ///< number of k-points
    //   int Nkpts_sym;///< number of k-points after symmetry reduction
       int Kx; ///< number of kpoints in x direction
       int Ky; ///< number of kpoints in y direction
       int Kz; ///< number of kpoints in z direction
  //   int kctr; ///< counter for k point. 
  //     PetscScalar *lambdakpt; ///< stores the eigenvalues in k point code. size: Nstates*Nkpts_sym
  //     PetscScalar *kptWts; ///< stores the weights of the k points
  //     PetscScalar *k1;
  //     PetscScalar *k2;
  //     PetscScalar *k3;
    

    /* /\* */
    /*  * variables for MD */
    /*  *\/ */
    /* PetscScalar *vels; */
    /* PetscScalar *vels_old; */
    /* PetscScalar *accels; */
    /* PetscScalar *Mass; */
    /* PetscScalar T_MD;           // md temperature */
    /* int MaxMDSteps; */
    /* int MDStepCount; */
    /* int MDFlag; */
    /* PetscScalar time_step; */
    /* PetscScalar TE; */
    /* PetscScalar KE; */
    /* PetscScalar PE; */
    /* PetscScalar TE_mean; */
    /* PetscScalar PE_mean; */
    /* PetscScalar KE_mean; */
    /* PetscScalar T_mean; */
    /* PetscScalar TE_std; */
    /* PetscScalar PE_std; */
    /* PetscScalar KE_std; */
    /* PetscScalar T_std; */

    /*
     * Vectors for Anderson mixing
     */
    Vec xkprev;
    Vec xk;
    Vec fkprev;
    Vec *Xk;
    Vec *Fk;
    Vec *XpbF;

    Vec tempVec;                ///< temporary vector used for storing copy of vector during computation 
    KSP ksp;                    ///< PETSc (Krylov subspace) ksp object for solving the Poisson equation. See PETSc manual for more details on ksp

    PSD_OBJ *psd;               ///< datatype for the pseudopotentials. This is defined for each type of atom.
    OVERLAP_OBJ *AtomOverlap_local;
    OVERLAP_OBJ *AtomOverlap_localforces;
    OVERLAP_OBJ *AtomOverlap_nonlocal;
    OVERLAP_OBJ *AtomOverlap_nonlocalforces;
    PROJECTOR_OBJ *Projector_nonlocal;
    PROJECTOR_OBJ *Projector_nonlocalforces;

  // data structures for preconditioner
  Mat HelmholtzOpr;
  KSP kspHelm;
  Vec rhs;
  Vec Pfk;
  Vec hk;
  int precond;
  
   
    /*
     * filenames
     */
    char OutFilename[100];
    char ForceFilename[100];
    char AtomFilename[100];
    char EigenFilename[100];
    char MDFilename[100];
    char DenFilename[100];
    char ReadDenFilename[100];

    /*
     * times
     */
    PetscLogDouble timeInitialization;  //X
    PetscLogDouble timePseudocharge;    //X
    PetscLogDouble timePoisson; //X
    PetscLogDouble timeDensity; //X
    PetscLogDouble timeEnergy;  //X
    PetscLogDouble timeForces;  //X 
    PetscLogDouble timeTotal;   //X
    PetscLogDouble timeLSSGQ;
    //PetscLogDouble timeSubspace;
    //PetscLogDouble timeRotation;
    //PetscLogDouble timeEigen;
    PetscLogDouble timeSCF;
    PetscLogDouble timeLanczos;
    PetscLogDouble timeFermienergy;

    PetscLogDouble DtimeLSSGQ;  // cheb filter time per SCF
  // PetscLogDouble DtimeSQ;     // subspace time per SCF
    //PetscLogDouble DtimeEigen; // eigen time per SCF
    //PetscLogDouble DtimeRotation; // rotation time per SCF
    PetscLogDouble DtimeDensity;    // density calculation time per SCF
    PetscLogDouble DtimeEnergy;
    PetscLogDouble DtimeForces;
    PetscLogDouble DtimePoisson;    // poisson time per SCF  
    PetscLogDouble DtimeLanczos;
    PetscLogDouble DtimeFermienergy;

    PetscLogDouble timeInitialization_start;
    PetscLogDouble timePseudocharge_start;
    PetscLogDouble timePoisson_start;
    PetscLogDouble timeDensity_start;
    PetscLogDouble timeEnergy_start;
    PetscLogDouble timeForces_start;
    PetscLogDouble timeTotal_start;
     PetscLogDouble timeLSSGQ_start;
    //PetscLogDouble timeSubspace_start;
    //PetscLogDouble timeRotation_start;
    //PetscLogDouble timeEigen_start;
    PetscLogDouble timeSCF_start;
    PetscLogDouble timeLanczos_start;
    PetscLogDouble timeFermienergy_start;

} LSDFT_OBJ;
#endif
