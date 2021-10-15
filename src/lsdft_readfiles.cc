/*=============================================================================================
  |
  | file name: readfiles.cc
  |
  | Description: This file contains the functions required for reading the input files
  |
  | Authors: Swarnava Ghosh
  |
  | Last Modified: 11/28/2016
  |-------------------------------------------------------------------------------------------*/
#include "lsdft.h"
#include "petscsys.h"
#include <iostream>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace std;
///////////////////////////////////////////////////////////////////////////////////////////////
//                              Read_Parameters: reads the .inpt file                        //
///////////////////////////////////////////////////////////////////////////////////////////////
void Read_Parameters(LSDFT_OBJ *pLsdft)
{

    /*
     * time stamp
     */
    char *c_time_str;
    time_t current_time = time(NULL);
    c_time_str = ctime(&current_time);
    PetscTime(&pLsdft->timeTotal_start);

    FILE *fConfFile;
    FILE *fOutFile;
    PetscInt p, a;
    int rank, size;
    char ConfigFile[100] = "./";
    strcat(ConfigFile, pLsdft->file);
    strcat(ConfigFile, ".inpt");

    pLsdft->kB = 8.617343e-5;   // Boltzmann constant in eV/K
    pLsdft->Ceh = 27.211384523; // conversion from eV to Ha, 1 Ha=27.211384523 eV, so 1eV=1/Ceh Ha

    /*
     * initialize time variables
     */
    pLsdft->timeInitialization = 0.0;
    pLsdft->timePoisson = 0.0;
    pLsdft->timeDensity = 0.0;
    pLsdft->timeForces = 0.0;
    pLsdft->timeEnergy = 0.0;
    pLsdft->timePseudocharge = 0.0;
    pLsdft->timeTotal = 0.0;
    pLsdft->timeLSSGQ = 0.0;
    // pLsdft->timeSubspace=0.0;
    //pLsdft->timeRotation=0.0;
    //pLsdft->timeEigen=0.0;
    pLsdft->timeDensity = 0.0;
    pLsdft->timeSCF = 0.0;
    pLsdft->timeLanczos = 0.0;
    pLsdft->timeFermienergy = 0.0;

    pLsdft->Ntype = -1;
    pLsdft->SCFNewRhoCalcCtr = 3;   // default value unless overwritten
    pLsdft->RelaxCount = 0;     // incremented as NLCG progresses
    //pLsdft->ChebyshevCallCounter=0; // number of times chebyshev filtering has been called
    pLsdft->MixingParameter = 0.3;  // default mixing parameter
    pLsdft->MixingHistory = 7;  // default mixing history
    pLsdft->RelaxFlag = 0;      // default: no relaxation
    // pLsdft->MDFlag = 0;         // default: no MD
    pLsdft->TOLSCF = 1e-6;      // default tolerance for SCF
    pLsdft->KSPTOL = 1e-8;      // default tolerance for KSP
    pLsdft->NLCGTOL = 1e-10;    // default tolerance for NLCG
    pLsdft->LANCZOSTOL = 1e-6;  // default tolerance for Lanczos
    pLsdft->PseudochargeRadiusTOL = 1e-6;   // default tolerance for pseudocharge density radius
    pLsdft->MAXITSCF = 100;     // default maximum number of SCF iterations
    pLsdft->MAXIT_NLCG = 300;   // default maximum number of iterations in NLCG
    pLsdft->order = 6;          // default finite-difference order is 12
    pLsdft->Beta = 1000;        // default smearing
    pLsdft->ChebDegree = 20; // default chebychev filter polynomial degree
    pLsdft->REFERENCE_CUTOFF = 0.5; // default cutoff for reference pseudopotential
    strcpy(pLsdft->XC, "LDA");  // default exchange correlation=lda
    pLsdft->BC = 1;             // default is Dirichlet boundary condition
    //pLsdft->Nkpts=1; // default number of k points of periodic calculation
    //pLsdft->Nkpts_sym=1;
    pLsdft->PulayFrequency = 1; // default mixing is pulay
    pLsdft->PulayRestartFlag = 0;   // default mixing is pulay, 1 for r-pulay
    //pLsdft->NkptsGroup=1; // number of k-point groups for parallelization
    //pLsdft->kptParalFlag=0; // k point parallelization is zero by default
    pLsdft->PrintForceFlag = 1; // flag for printing forces
    pLsdft->PrintAtomPosFlag = 0;   // flag for printing atomic positions
    pLsdft->PrintElecDensFlag = 0;  // flag for printing final electron density
    pLsdft->ReadElecDensFlag = 0;   // flag for reading final electron density
    pLsdft->PrintEigenFlag = 0; // Flag for printing final eigenvalues and occupation
    //pLsdft->Nkpts_symGroup=1; // number of kpoints (symmetricized) in each k-point group
    //pLsdft->NprocGroup=1; // number of processors in each group
    // pLsdft->time_step = 1.0;    // default time step is 1 femtosecond
    //  pLsdft->MaxMDSteps = 1000;  // default number of timesteps for MD
    pLsdft->N_qp = 50;          // number of lssgq quadrature points
    pLsdft->NplCC = 100;         // number of clenshaw curtis quadrature points
    pLsdft->Rcut = 12.00;       // Radius of truncation in Bohr
    pLsdft->LinInterpEFieldsFlag=0; // DO not perform linear interpolation of electronic fields
    pLsdft->precond=1; // preconditioning=1

    /*
     * default Monkhorst-Pack grid is 1x1x1
     */
    //pLsdft->Kx=1;
    // pLsdft->Ky=1;
    //pLsdft->Kz=1;

    /*
     * default file names (unless otherwise supplied)
     */
    strcpy(pLsdft->OutFilename, "./");
    strcpy(pLsdft->ForceFilename, "./");
    strcpy(pLsdft->AtomFilename, "./");
    strcpy(pLsdft->EigenFilename, "./");
    //  strcpy(pLsdft->MDFilename, "./");
    strcpy(pLsdft->DenFilename, "./");

    strcat(pLsdft->OutFilename, pLsdft->file);
    strcat(pLsdft->OutFilename, ".out");

    strcat(pLsdft->ForceFilename, pLsdft->file);
    strcat(pLsdft->ForceFilename, ".force");

    strcat(pLsdft->AtomFilename, pLsdft->file);
    strcat(pLsdft->AtomFilename, ".atom");

    strcat(pLsdft->EigenFilename, pLsdft->file);
    strcat(pLsdft->EigenFilename, ".eigen");

    // strcat(pLsdft->MDFilename, pLsdft->file);
    // strcat(pLsdft->MDFilename, ".md");

    strcat(pLsdft->DenFilename, pLsdft->file);
    char fname[100], denfname[100]; // output file name
    strcpy(fname, pLsdft->file);

    int i;
    if ((fConfFile = fopen((const char *) ConfigFile, "rb")) == NULL) {
        cout << "Couldn't open config file " << ConfigFile << '\n';
        cout << "Exiting for this configuration...\n";
        exit(1);
    }
    char str[60];

    do {
        fscanf(fConfFile, "%s", str);
        if (strcmp(str, "FD_ORDER:") == 0) {
            fscanf(fConfFile, "%d", &pLsdft->order);
            pLsdft->order = (pLsdft->order) / 2;    // store half order
        } else if (strcmp(str, "BOUNDARY_CONDITION:") == 0) {
            fscanf(fConfFile, "%d \n", &pLsdft->BC);
        } else if (strcmp(str, "CELL:") == 0) {
            fscanf(fConfFile, "%lf", &pLsdft->range_x);
            fscanf(fConfFile, "%lf", &pLsdft->range_y);
            fscanf(fConfFile, "%lf", &pLsdft->range_z);
        } else if (strcmp(str, "FD_GRID:") == 0) {
            fscanf(fConfFile, "%d", &pLsdft->numPoints_x);
            fscanf(fConfFile, "%d", &pLsdft->numPoints_y);
            fscanf(fConfFile, "%d", &pLsdft->numPoints_z);

            // default value
            pLsdft->SCnumPoints_x = pLsdft->numPoints_x;
            pLsdft->SCnumPoints_y = pLsdft->numPoints_y;
            pLsdft->SCnumPoints_z = pLsdft->numPoints_z;

        }else if(strcmp(str,"KPOINT_GRID:")==0)
          {
            fscanf(fConfFile,"%d",&pLsdft->Kx);
            fscanf(fConfFile,"%d",&pLsdft->Ky);
            fscanf(fConfFile,"%d",&pLsdft->Kz);
            //pLsdft->Nkpts = pLsdft->Kx*pLsdft->Ky*pLsdft->Kz;
            /*
             * currently only takes into account the time reversal symmetry
             */
            //pLsdft->Nkpts_sym = ceil(pLsdft->Kx/2.0)*pLsdft->Ky*pLsdft->Kz;
          }
        // else if(strcmp(str,"KPOINT_PARAL_GROUP:")==0)
        //   {
        //     fscanf(fConfFile,"%d",&pLsdft->NkptsGroup);
        //     if(pLsdft->NkptsGroup > 1)
        //       {
        //         pLsdft->kptParalFlag=1;
        //       }
        //   }
        else if (strcmp(str, "BETA:") == 0) {
            fscanf(fConfFile, "%lf", &pLsdft->Beta);
        } else if (strcmp(str, "RADIUS_TRUNC:") == 0) {
            fscanf(fConfFile, "%lf", &pLsdft->Rcut);
        }
        else if(strcmp(str,"CHEB_DEGREE:")==0)
           {
             fscanf(fConfFile,"%d", &pLsdft->NplCC);
           }
	else if(strcmp(str,"CHEB_FILTER_ORDER:")==0)
           {
             fscanf(fConfFile,"%d", &pLsdft->ChebDegree);
           }
        else if (strcmp(str, "LANCZOS_DEGREE:") == 0) {
            fscanf(fConfFile, "%d \n", &pLsdft->N_qp);
        }
        else if (strcmp(str, "EXCHANGE_CORRELATION:") == 0) {
            fscanf(fConfFile, "%s", pLsdft->XC);
        }
        else if(strcmp(str,"NSTATES:")==0)
           {
             fscanf(fConfFile,"%d",&pLsdft->NstatesBC);
             /*
              * allocate memory for storing eigenvalues
              */
             
           }
        else if (strcmp(str, "ION_RELAX:") == 0) {
            fscanf(fConfFile, "%d", &pLsdft->RelaxFlag);
        } else if (strcmp(str, "SCF_RHO_CALC_COUNT:") == 0) {
            fscanf(fConfFile, "%d", &pLsdft->SCFNewRhoCalcCtr);
        } else if (strcmp(str, "MIXING_PARAMETER:") == 0) {
            fscanf(fConfFile, "%lf", &pLsdft->MixingParameter);
        } else if (strcmp(str, "PULAY_FREQUENCY:") == 0) {
            fscanf(fConfFile, "%d", &pLsdft->PulayFrequency);
        } else if (strcmp(str, "PULAY_RESTART:") == 0) {
            fscanf(fConfFile, "%d", &pLsdft->PulayRestartFlag);
        } else if (strcmp(str, "REFERENCE_CUTOFF:") == 0) {
            fscanf(fConfFile, "%lf", &pLsdft->REFERENCE_CUTOFF);
        } else if (strcmp(str, "MIXING_HISTORY:") == 0) {
            fscanf(fConfFile, "%d", &pLsdft->MixingHistory);
        } else if (strcmp(str, "TOL_SCF:") == 0) {
            fscanf(fConfFile, "%lf", &pLsdft->TOLSCF);
        } else if (strcmp(str, "TOL_POISSON:") == 0) {
            fscanf(fConfFile, "%lf", &pLsdft->KSPTOL);
        } else if (strcmp(str, "TOL_NLCG:") == 0) {
            fscanf(fConfFile, "%lf", &pLsdft->NLCGTOL);
        } else if (strcmp(str, "TOL_LANCZOS:") == 0) {
            fscanf(fConfFile, "%lf", &pLsdft->LANCZOSTOL);
        } else if (strcmp(str, "TOL_PSEUDOCHARGE:") == 0) {
            fscanf(fConfFile, "%lf", &pLsdft->PseudochargeRadiusTOL);
        } else if (strcmp(str, "MAXIT_SCF:") == 0) {
            fscanf(fConfFile, "%d", &pLsdft->MAXITSCF);
        } else if (strcmp(str, "MAXIT_NLCG:") == 0) {
            fscanf(fConfFile, "%d", &pLsdft->MAXIT_NLCG);
        } else if (strcmp(str, "PRINT_FORCES:") == 0) {
            fscanf(fConfFile, "%d", &pLsdft->PrintForceFlag);
        } else if (strcmp(str, "PRINT_ATOMS:") == 0) {
            fscanf(fConfFile, "%d", &pLsdft->PrintAtomPosFlag);
        } else if (strcmp(str, "PRINT_EIGEN:") == 0) {
            fscanf(fConfFile, "%d", &pLsdft->PrintEigenFlag);
        } else if (strcmp(str, "PRINT_DENSITY:") == 0) {
            fscanf(fConfFile, "%d", &pLsdft->PrintElecDensFlag);
        } else if (strcmp(str, "CALCULATE_DIPOLE:") == 0) {
            fscanf(fConfFile, "%d", &pLsdft->DipoleMomentFlag);
        } else if (strcmp(str, "NTYPE:") == 0) {
            fscanf(fConfFile, "%d", &pLsdft->Ntype);
            /*
             * allocate memory
             */
            PetscMalloc(sizeof(PetscScalar) * pLsdft->Ntype, &pLsdft->CUTOFF_x);
            PetscMalloc(sizeof(PetscScalar) * pLsdft->Ntype, &pLsdft->CUTOFF_y);
            PetscMalloc(sizeof(PetscScalar) * pLsdft->Ntype, &pLsdft->CUTOFF_z);
            PetscMalloc(sizeof(PetscInt) * pLsdft->Ntype, &pLsdft->startPos);
            PetscMalloc(sizeof(PetscInt) * pLsdft->Ntype, &pLsdft->endPos);
	    PetscMalloc(sizeof(PetscInt) * pLsdft->Ntype, &pLsdft->startPosSC);
            PetscMalloc(sizeof(PetscInt) * pLsdft->Ntype, &pLsdft->endPosSC);
            PetscMalloc(sizeof(PetscInt) * pLsdft->Ntype, &pLsdft->noa);
            PetscMalloc(sizeof(PetscInt) * pLsdft->Ntype, &pLsdft->noe);
	    PetscMalloc(sizeof(PetscInt) * pLsdft->Ntype, &pLsdft->noaSC);
            PetscMalloc(sizeof(PetscInt) * pLsdft->Ntype, &pLsdft->noeSC);
            PetscMalloc(sizeof(PetscInt) * pLsdft->Ntype, &pLsdft->localPsd);
           
            pLsdft->psdName = (char **) malloc(pLsdft->Ntype * sizeof(char *));
            for (i = 0; i < pLsdft->Ntype; i++) {
                pLsdft->psdName[i] = (char *) malloc(25 + 1);
            }
            pLsdft->atomType = (char **) malloc(pLsdft->Ntype * sizeof(char *));
            for (i = 0; i < pLsdft->Ntype; i++) {
                pLsdft->atomType[i] = (char *) malloc(10 + 1);
            }
        } else if (strcmp(str, "PSEUDOPOTENTIAL_LOCAL:") == 0) {
            if (pLsdft->Ntype <= 0) {
                cout << "Input number of types of elements first\n";
                exit(1);
            }
            for (i = 0; i < pLsdft->Ntype; i++) {
                fscanf(fConfFile, "%d", &pLsdft->localPsd[i]);
            }
        } else if (strcmp(str, "PSEUDOPOTENTIAL_FILE:") == 0) {
            if (pLsdft->Ntype <= 0) {
                cout << "Input number of types of elements first\n";
                exit(1);
            }
            for (i = 0; i < pLsdft->Ntype; i++) {
                fscanf(fConfFile, "%s", pLsdft->psdName[i]);
            }
        } else if (strcmp(str, "OUTPUT_FILE:") == 0) {
            fscanf(fConfFile, "%s", fname);

            strcpy(pLsdft->OutFilename, "./");
            strcpy(pLsdft->ForceFilename, "./");
            strcpy(pLsdft->AtomFilename, "./");
            strcpy(pLsdft->EigenFilename, "./");
            strcpy(pLsdft->MDFilename, "./");
            strcpy(pLsdft->DenFilename, "./");


            strcat(pLsdft->OutFilename, fname);
            strcat(pLsdft->OutFilename, ".out");

            strcat(pLsdft->ForceFilename, fname);
            strcat(pLsdft->ForceFilename, ".force");

            strcat(pLsdft->AtomFilename, fname);
            strcat(pLsdft->AtomFilename, ".atom");

            strcat(pLsdft->EigenFilename, fname);
            strcat(pLsdft->EigenFilename, ".eigen");

            strcat(pLsdft->MDFilename, fname);
            strcat(pLsdft->MDFilename, ".md");

            strcat(pLsdft->DenFilename, fname);

        } else if (strcmp(str, "READ_DEN_FILE:") == 0) {
            fscanf(fConfFile, "%s", denfname);
            strcpy(pLsdft->ReadDenFilename, "./");
            strcat(pLsdft->ReadDenFilename, denfname);
            pLsdft->ReadElecDensFlag = 1;
        }


    } while (!feof(fConfFile));

    fclose(fConfFile);

   
   
    /*
     * initialization starts here
     */

    PetscTime(&pLsdft->timeInitialization_start);

    if (pLsdft->RelaxFlag == 0) {
        pLsdft->MAXIT_NLCG = 0;
    }

    // first update Rcut such that 2(Rcut_x+range_x) is divisible by the mesh
    // pLsdft->Rcut_x=pLsdft->delta_x*ceil(2.0*(pLsdft->Rcut+pLsdft->range_x)/pLsdft->delta_x) - 2.0*pLsdft->range_x;
    // pLsdft->Rcut_y=pLsdft->delta_y*ceil(2.0*(pLsdft->Rcut+pLsdft->range_y)/pLsdft->delta_y) - 2.0*pLsdft->range_y;
    // pLsdft->Rcut_z=pLsdft->delta_z*ceil(2.0*(pLsdft->Rcut+pLsdft->range_z)/pLsdft->delta_z) - 2.0*pLsdft->range_z;


    /*
     * now that the .inpt file is read, we calculate the mesh spacing
     */
    if (pLsdft->BC == 1)        // non-periodic boundary condition
    {
        pLsdft->delta_x = 2.0 * pLsdft->range_x / (pLsdft->numPoints_x - 1);
        pLsdft->delta_y = 2.0 * pLsdft->range_y / (pLsdft->numPoints_y - 1);
        pLsdft->delta_z = 2.0 * pLsdft->range_z / (pLsdft->numPoints_z - 1);

        // pLsdft->delta_x = 2.0 * pLsdft->range_x / (pLsdft->numPoints_x);
        // pLsdft->delta_y = 2.0 * pLsdft->range_y / (pLsdft->numPoints_y);
        // pLsdft->delta_z = 2.0 * pLsdft->range_z / (pLsdft->numPoints_z);

    } else if (pLsdft->BC == 2) // 3D periodic boundary condition
    {
        pLsdft->delta_x = 2.0 * pLsdft->range_x / (pLsdft->numPoints_x);
        pLsdft->delta_y = 2.0 * pLsdft->range_y / (pLsdft->numPoints_y);
        pLsdft->delta_z = 2.0 * pLsdft->range_z / (pLsdft->numPoints_z);

        // // additional number of points for supercell
        // pLsdft->SCnumPoints_x = pLsdft->numPoints_x + 2.0 * pLsdft->delta_x * ceil(pLsdft->Rcut / pLsdft->delta_x);
        // pLsdft->SCnumPoints_y = pLsdft->numPoints_y + 2.0 * pLsdft->delta_y * ceil(pLsdft->Rcut / pLsdft->delta_y);
        // pLsdft->SCnumPoints_z = pLsdft->numPoints_z + 2.0 * pLsdft->delta_z * ceil(pLsdft->Rcut / pLsdft->delta_z);

        // if(pLsdft->Nkpts > 1) // allocate memory for storing eigenvalues of all kpoints
        //    {
        //      PetscMalloc(sizeof(PetscScalar*)*(pLsdft->Nkpts_sym),&pLsdft->kptWts);
        //      PetscMalloc(sizeof(PetscScalar*)*(pLsdft->Nkpts_sym),&pLsdft->k1);
        //      PetscMalloc(sizeof(PetscScalar*)*(pLsdft->Nkpts_sym),&pLsdft->k2);
        //      PetscMalloc(sizeof(PetscScalar*)*(pLsdft->Nkpts_sym),&pLsdft->k3);
        //      PetscMalloc(sizeof(PetscScalar*)*(pLsdft->Nkpts_sym*pLsdft->Nstates),&pLsdft->lambdakpt);
        //      /*
        //       * calculate and store the k point weights
        //       */
        //      Calculate_kpoints(pLsdft);
        //    }
    } else if (pLsdft->BC == 3) // 2D periodic boundary condition
    {
        pLsdft->delta_x = 2.0 * pLsdft->range_x / (pLsdft->numPoints_x);
        pLsdft->delta_y = 2.0 * pLsdft->range_y / (pLsdft->numPoints_y);
        pLsdft->delta_z = 2.0 * pLsdft->range_z / (pLsdft->numPoints_z - 1);

        // // additional number of points for supercell
        // pLsdft->SCnumPoints_x = pLsdft->numPoints_x + 2.0 * pLsdft->delta_x * ceil(pLsdft->Rcut / pLsdft->delta_x);
        // pLsdft->SCnumPoints_y = pLsdft->numPoints_y + 2.0 * pLsdft->delta_y * ceil(pLsdft->Rcut / pLsdft->delta_y);
        // pLsdft->SCnumPoints_z = pLsdft->numPoints_z;

        // if(pLsdft->Nkpts > 1) // allocate memory for storing eigenvalues of all kpoints
        //    {
        //      PetscMalloc(sizeof(PetscScalar*)*(pLsdft->Nkpts_sym),&pLsdft->kptWts);
        //      PetscMalloc(sizeof(PetscScalar*)*(pLsdft->Nkpts_sym),&pLsdft->k1);
        //      PetscMalloc(sizeof(PetscScalar*)*(pLsdft->Nkpts_sym),&pLsdft->k2);
        //      PetscMalloc(sizeof(PetscScalar*)*(pLsdft->Nkpts_sym),&pLsdft->k3);
        //      PetscMalloc(sizeof(PetscScalar*)*(pLsdft->Nkpts_sym*pLsdft->Nstates),&pLsdft->lambdakpt);
        //      /*
        //       * calculate and store the k point weights
        //       */
        //       Calculate_kpoints(pLsdft);
        //    }
    } else if (pLsdft->BC == 4) // 1D periodic boundary condition
    {
        pLsdft->delta_x = 2.0 * pLsdft->range_x / (pLsdft->numPoints_x);
        pLsdft->delta_y = 2.0 * pLsdft->range_y / (pLsdft->numPoints_y - 1);
        pLsdft->delta_z = 2.0 * pLsdft->range_z / (pLsdft->numPoints_z - 1);

        // pLsdft->SCnumPoints_x = pLsdft->numPoints_x + 2.0 * pLsdft->delta_x * ceil(pLsdft->Rcut / pLsdft->delta_x);
        // pLsdft->SCnumPoints_y = pLsdft->numPoints_y;
        // pLsdft->SCnumPoints_z = pLsdft->numPoints_z;

        // if(pLsdft->Nkpts > 1) // allocate memory for storing eigenvalues of all kpoints
        //    {
        //      PetscMalloc(sizeof(PetscScalar*)*(pLsdft->Nkpts_sym),&pLsdft->kptWts);
        //      PetscMalloc(sizeof(PetscScalar*)*(pLsdft->Nkpts_sym),&pLsdft->k1);
        //      PetscMalloc(sizeof(PetscScalar*)*(pLsdft->Nkpts_sym),&pLsdft->k2);
        //      PetscMalloc(sizeof(PetscScalar*)*(pLsdft->Nkpts_sym),&pLsdft->k3);
        //      PetscMalloc(sizeof(PetscScalar*)*(pLsdft->Nkpts_sym*pLsdft->Nstates),&pLsdft->lambdakpt);
        //      /*
        //       * calculate and store the k point weights
        //       */
        //       Calculate_kpoints(pLsdft);
        //    }
    }


     /// / additional number of points for supercell (always)
      //   pLsdft->SCnumPoints_x = pLsdft->numPoints_x + 2 * ceil(pLsdft->Rcut / pLsdft->delta_x);
      //   pLsdft->SCnumPoints_y = pLsdft->numPoints_y + 2 * ceil(pLsdft->Rcut / pLsdft->delta_y);
      //   pLsdft->SCnumPoints_z = pLsdft->numPoints_z + 2 * ceil(pLsdft->Rcut / pLsdft->delta_z);

    pLsdft->nx_loc=2*ceil(pLsdft->Rcut/pLsdft->delta_x)+1;
    pLsdft->ny_loc=2*ceil(pLsdft->Rcut/pLsdft->delta_y)+1;
    pLsdft->nz_loc=2*ceil(pLsdft->Rcut/pLsdft->delta_z)+1;

    pLsdft->delVol = pLsdft->delta_x * pLsdft->delta_y * pLsdft->delta_z;
    pLsdft->HalfdelVol = 0.5 * pLsdft->delVol;
    pLsdft->OrbScalingFactor = 1.0 / (sqrt(pLsdft->delVol));

    /*
     * finite difference coefficients of -(1/2)*laplacian operator
     */
    pLsdft->coeffs_x[0] = 0;
    pLsdft->coeffs_y[0] = 0;
    pLsdft->coeffs_z[0] = 0;
    for (a = 1; a <= pLsdft->order; a++) {
        pLsdft->coeffs_x[0] += ((PetscScalar) 1.0 / (a * a));
    }
    pLsdft->coeffs_y[0] = pLsdft->coeffs_x[0];
    pLsdft->coeffs_z[0] = pLsdft->coeffs_x[0];

    pLsdft->coeffs_x[0] *= ((PetscScalar) 1.0 / (pLsdft->delta_x * pLsdft->delta_x));
    pLsdft->coeffs_y[0] *= ((PetscScalar) 1.0 / (pLsdft->delta_y * pLsdft->delta_y));
    pLsdft->coeffs_z[0] *= ((PetscScalar) 1.0 / (pLsdft->delta_z * pLsdft->delta_z));
    //printf("pLsdft->coeffs_x[0]=%lf    ", pLsdft->coeffs_x[0]);

    for (p = 1; p <= pLsdft->order; p++) {
        pLsdft->coeffs_x[p] =
            (PetscScalar) (-1 * pow(-1, p + 1) * fract(pLsdft->order, p) / (p * p * pLsdft->delta_x * pLsdft->delta_x));
        pLsdft->coeffs_y[p] =
            (PetscScalar) (-1 * pow(-1, p + 1) * fract(pLsdft->order, p) / (p * p * pLsdft->delta_y * pLsdft->delta_y));
        pLsdft->coeffs_z[p] =
            (PetscScalar) (-1 * pow(-1, p + 1) * fract(pLsdft->order, p) / (p * p * pLsdft->delta_z * pLsdft->delta_z));

	//  printf("pLsdft->coeffs_x[%d]=%lf    ", p, pLsdft->coeffs_x[p]);
    }

    /*
     * finite difference coefficients of -(1/4*Pi)*laplacian operator
     */
    for (p = 0; p <= pLsdft->order; p++) {
        pLsdft->stencil_coeffs_x[p] = pLsdft->coeffs_x[p] / (2.0 * M_PI);
        pLsdft->stencil_coeffs_y[p] = pLsdft->coeffs_y[p] / (2.0 * M_PI);
        pLsdft->stencil_coeffs_z[p] = pLsdft->coeffs_z[p] / (2.0 * M_PI);
    }

    /*
     * finite difference coefficients for the gradient operator
     */
    pLsdft->coeffs_grad_x[0] = 0;
    pLsdft->coeffs_grad_y[0] = 0;
    pLsdft->coeffs_grad_z[0] = 0;
    for (p = 1; p <= pLsdft->order; p++) {
        pLsdft->coeffs_grad_x[p] = (PetscScalar) (pow(-1, p + 1) * fract(pLsdft->order, p) / (p * pLsdft->delta_x));
        pLsdft->coeffs_grad_y[p] = (PetscScalar) (pow(-1, p + 1) * fract(pLsdft->order, p) / (p * pLsdft->delta_y));
        pLsdft->coeffs_grad_z[p] = (PetscScalar) (pow(-1, p + 1) * fract(pLsdft->order, p) / (p * pLsdft->delta_z));
    }

    /*
     * print initialization details to output file (only rank 0)
     */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /*
     * assert k-point parallelization is indeed for a k-point calculation
     */
    // if(pLsdft->kptParalFlag==1)
    //  {
    //    if(pLsdft->Nkpts_sym==1)
    //    {
    //      PetscPrintf(PETSC_COMM_WORLD,"Number of symmetricized k points is 1. k-point parallelization is not needed.");
    //    exit(1);
    //    }
    //  }

    // pLsdft->NprocGroup = size/pLsdft->NkptsGroup;  //number of processors in each k-point group
    // pLsdft->Nkpts_symGroup =  pLsdft->Nkpts_sym/pLsdft->NkptsGroup; //number of symmetricized k-points per group

    if (rank == 0) {
        if ((fOutFile = fopen((const char *) pLsdft->OutFilename, "w")) == NULL) {
            cout << "Couldn't create output file...\n" << pLsdft->OutFilename;
            cout << "Exiting...\n";
            exit(1);
        }
        fprintf(fOutFile, "***************************************************************************\n");
        fprintf(fOutFile, "*                            LSDFT-0.2d                                    *\n");
        fprintf(fOutFile, "*                Linear Scaling Density Functional Theory                 *\n");
        fprintf(fOutFile, "*                      Copyright (C) 2018 Caltech                         *\n");
        fprintf(fOutFile, "*      Developed from SPARC-0.2, originally distributed under GNU (GPL) 3 *\n");
        fprintf(fOutFile, "*                 Start time: %0.24s                    *\n", c_time_str);
        fprintf(fOutFile, "***************************************************************************\n");
        fprintf(fOutFile, "                           Input parameters                                \n");
        fprintf(fOutFile, "***************************************************************************\n");
        fprintf(fOutFile, "CELL: %f %f %f \n", pLsdft->range_x, pLsdft->range_y, pLsdft->range_z);
        fprintf(fOutFile, "FD_GRID: %d %d %d\n", pLsdft->numPoints_x, pLsdft->numPoints_y, pLsdft->numPoints_z);
        fprintf(fOutFile, "FD_ORDER: %d\n", 2 * pLsdft->order);
        fprintf(fOutFile, "BOUNDARY_CONDITION: %d\n", pLsdft->BC);
        // if(pLsdft->BC==2 || pLsdft->BC==3 || pLsdft->BC==4)
        //   {
        // fprintf(fOutFile,"KPOINT_GRID: %d %d %d\n",pLsdft->Kx,pLsdft->Ky,pLsdft->Kz);
        /*
         * check k point parallelization parameters
         */
        // if(pLsdft->kptParalFlag==1)
        //   {
        //     fprintf(fOutFile,"KPOINT_PARAL_GROUP: %d\n",pLsdft->NkptsGroup);
        //     if(size%pLsdft->NkptsGroup!=0 || pLsdft->Nkpts_sym%pLsdft->NkptsGroup !=0)
        //    {
        //      fprintf(fOutFile,"KPOINT_PARAL_GROUP is not consistent with number of processors and number of k-points.\n Choose appropriate parameters and rerun.\n");
        //      fprintf(fOutFile,"size=%d, pLsdft->NkptsGroup=%d, pLsdft->Nkpts_sym=%d\n",size,pLsdft->NkptsGroup,pLsdft->Nkpts_sym);
        //      exit(1);
        //    }
        //   }

        // }

        fprintf(fOutFile, "BETA: %lf\n", pLsdft->Beta);
        fprintf(fOutFile,"CHEB_DEGREE: %d\n",pLsdft->NplCC);
	fprintf(fOutFile,"LANCZOS_DEGREE: %d\n",pLsdft->N_qp);
	fprintf(fOutFile, "RADIUS_TRUNC: %lf\n", pLsdft->Rcut);
        fprintf(fOutFile,"NSTATES: %d\n",pLsdft->NstatesBC);
        fprintf(fOutFile, "NTYPE: %d\n", pLsdft->Ntype);
        fprintf(fOutFile, "EXCHANGE_CORRELATION: %s\n", pLsdft->XC);
        fprintf(fOutFile, "ION_RELAX: %d\n", pLsdft->RelaxFlag);
        fprintf(fOutFile, "MAXIT_SCF: %d\n", pLsdft->MAXITSCF);
        fprintf(fOutFile, "TOL_SCF: %E\n", pLsdft->TOLSCF);
        fprintf(fOutFile, "TOL_POISSON: %E\n", pLsdft->KSPTOL);
        fprintf(fOutFile, "TOL_LANCZOS: %E\n", pLsdft->LANCZOSTOL);
        fprintf(fOutFile, "TOL_PSEUDOCHARGE: %E\n", pLsdft->PseudochargeRadiusTOL);
        fprintf(fOutFile, "MIXING_PARAMETER: %f\n", pLsdft->MixingParameter);
        fprintf(fOutFile, "MIXING_HISTORY: %d\n", pLsdft->MixingHistory);
        fprintf(fOutFile, "PULAY_FREQUENCY: %d\n", pLsdft->PulayFrequency);
        fprintf(fOutFile, "PULAY_RESTART: %d\n", pLsdft->PulayRestartFlag);
        fprintf(fOutFile, "REFERENCE_CUTOFF: %f\n", pLsdft->REFERENCE_CUTOFF);
        fprintf(fOutFile, "SCF_RHO_CALC_COUNT: %d\n", pLsdft->SCFNewRhoCalcCtr);
        fprintf(fOutFile, "PRINT_FORCES: %d\n", pLsdft->PrintForceFlag);
        fprintf(fOutFile, "PRINT_ATOMS: %d\n", pLsdft->PrintAtomPosFlag);
        fprintf(fOutFile, "PRINT_EIGEN: %d\n", pLsdft->PrintEigenFlag);
        fprintf(fOutFile, "PRINT_DENSITY: %d\n", pLsdft->PrintElecDensFlag);
        if (pLsdft->DipoleMomentFlag == 1) {
            fprintf(fOutFile, "CALCULATE_DIPOLE: %d\n", pLsdft->DipoleMomentFlag);
        }
        if (pLsdft->ReadElecDensFlag == 1) {
            fprintf(fOutFile, "READ_DEN_FILE: %s\n", pLsdft->ReadDenFilename);
        }
        if (pLsdft->RelaxFlag == 1) {
            fprintf(fOutFile, "TOL_NLCG: %e\n", pLsdft->NLCGTOL);
        }

        fprintf(fOutFile, "PSEUDOPOTENTIAL_FILE:");
        for (i = 0; i < pLsdft->Ntype; i++) {
            fprintf(fOutFile, " %s", pLsdft->psdName[i]);
        }
        fprintf(fOutFile, "\n");

        fprintf(fOutFile, "PSEUDOPOTENTIAL_LOCAL:");
        for (i = 0; i < pLsdft->Ntype; i++) {
            fprintf(fOutFile, " %d", pLsdft->localPsd[i]);
        }
        fprintf(fOutFile, "\n");

        fprintf(fOutFile, "OUTPUT_FILE: %s\n", fname);

        fprintf(fOutFile, "***************************************************************************\n");
        fprintf(fOutFile, "                             Initialization                                \n");
        fprintf(fOutFile, "***************************************************************************\n");
        fprintf(fOutFile, "Number of processors               :  %d\n", size);
        if ((fabs(pLsdft->delta_x - pLsdft->delta_y) <= 1e-15) && (fabs(pLsdft->delta_x - pLsdft->delta_z) <= 1e-15)
            && (fabs(pLsdft->delta_y - pLsdft->delta_z) <= 1e-15)) {
            fprintf(fOutFile, "Mesh spacing                       : % E (Bohr)\n", pLsdft->delta_x);
        } else {
            fprintf(fOutFile, "Mesh spacing in x-direction        : % E (Bohr)\n", pLsdft->delta_x);
            fprintf(fOutFile, "Mesh spacing in y-direction        : % E (Bohr)\n", pLsdft->delta_y);
            fprintf(fOutFile, "Mesh spacing in z direction        : % E (Bohr)\n", pLsdft->delta_z);
        }
        // if(pLsdft->BC==2 || pLsdft->BC==3 || pLsdft->BC==4)
        //  {
        //  fprintf(fOutFile,"Number of symmetry adapted k-points: %d\n",pLsdft->Nkpts_sym);
        //  fprintf(fOutFile,"Number of k-point parallel groups  : %d\n",pLsdft->NkptsGroup);
        //  fprintf(fOutFile,"Number of k-points in each group   : %d\n",pLsdft->Nkpts_symGroup);
        //  fprintf(fOutFile,"Number of processors for each group: %d\n",pLsdft->NprocGroup);
        //  }
        fprintf(fOutFile, "Output printed to                  :  %s\n", pLsdft->OutFilename);
        if (pLsdft->PrintAtomPosFlag == 1) {
            fprintf(fOutFile, "Atom positions printed to          :  %s\n", pLsdft->AtomFilename);
        }
        if (pLsdft->PrintForceFlag == 1) {
            fprintf(fOutFile, "Forces printed to                  :  %s\n", pLsdft->ForceFilename);
        }
        if (pLsdft->PrintEigenFlag == 1) {
            fprintf(fOutFile, "Final eigenvalues  printed to      :  %s\n", pLsdft->EigenFilename);
        }
        
        fclose(fOutFile);

        /*
         * create atom positions and forces output files
         */
        if (pLsdft->PrintAtomPosFlag == 1) {
            if ((fOutFile = fopen((const char *) pLsdft->AtomFilename, "w")) == NULL) {
                cout << "Couldn't create output file for atom positions...\n" << pLsdft->AtomFilename;
                cout << "Exiting...\n";
                exit(1);
            }
            fclose(fOutFile);
        }

        if (pLsdft->PrintForceFlag == 1) {
            if ((fOutFile = fopen((const char *) pLsdft->ForceFilename, "w")) == NULL) {
                cout << "Couldn't create output file for Forces...\n" << pLsdft->ForceFilename;
                cout << "Exiting...\n";
                exit(1);
            }
            fclose(fOutFile);
        }

        if (pLsdft->PrintEigenFlag == 1) {
            if ((fOutFile = fopen((const char *) pLsdft->EigenFilename, "w")) == NULL) {
                cout << "Couldn't create output file for eigenvalues...\n" << pLsdft->EigenFilename;
                cout << "Exiting...\n";
                exit(1);
            }
            fclose(fOutFile);
        }

        


    }
    return;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//                     Read_Ion: Reads the .ion file for ionic positions                     //
///////////////////////////////////////////////////////////////////////////////////////////////
void Read_Ion(LSDFT_OBJ *pLsdft)
{
    FILE *fConfFile;
    char ConfigFile[100] = "./";
    strcat(ConfigFile, pLsdft->file);
    strcat(ConfigFile, ".ion");
    pLsdft->nAtoms = 0;
    PetscInt poscnt, index = 0;
    PetscScalar x0, y0, z0;
    int rank;

    pLsdft->noetot = 0;

    int at = 0;
    if ((fConfFile = fopen((const char *) ConfigFile, "rb")) == NULL) {
        cout << "Couldn't open config file...\n" << ConfigFile;
        cout << "Exiting for this configuration...\n";
        exit(1);
    }

    /*
     * first read total number of atoms
     */
    fscanf(fConfFile, "%d", &pLsdft->nAtoms);

    /*
     * create local vectors for storing forces and atom positions
     */
    VecCreate(PETSC_COMM_SELF, &pLsdft->forces);
    VecSetSizes(pLsdft->forces, PETSC_DECIDE, 3 * pLsdft->nAtoms);
    VecSetFromOptions(pLsdft->forces);
    VecDuplicate(pLsdft->forces, &pLsdft->Atompos);
    VecDuplicate(pLsdft->forces, &pLsdft->nonlocalforces);
    VecZeroEntries(pLsdft->forces);
    VecZeroEntries(pLsdft->nonlocalforces);

    pLsdft->noetot = 0;
    do {
        if (fscanf(fConfFile, "%s", pLsdft->atomType[at]) == EOF)   // atom type (string)
            break;
        if (fscanf(fConfFile, "%d", &pLsdft->noe[at]) == EOF)   // number of electrons
            break;
        if (fscanf(fConfFile, "%d", &pLsdft->noa[at]) == EOF)   // number of atoms
            break;
        pLsdft->startPos[at] = index;

        /*
         * loop over number of atoms of a particular type and read their positions
         */
        for (poscnt = 0; poscnt < pLsdft->noa[at]; poscnt++)
	  {
            pLsdft->noetot += pLsdft->noe[at];
            fscanf(fConfFile, "%lf", &x0);
            fscanf(fConfFile, "%lf", &y0);
            fscanf(fConfFile, "%lf", &z0);

            VecSetValues(pLsdft->Atompos, 1, &index, &x0, INSERT_VALUES);
            index++;
            VecSetValues(pLsdft->Atompos, 1, &index, &y0, INSERT_VALUES);
            index++;
            VecSetValues(pLsdft->Atompos, 1, &index, &z0, INSERT_VALUES);
            index++;
        }
        pLsdft->endPos[at] = index - 1;

        at++;
    } while (!feof(fConfFile));

    fclose(fConfFile);
    /*
     * open output file to write data
     */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        if ((fConfFile = fopen((const char *) pLsdft->OutFilename, "a")) == NULL) {
            cout << "Couldn't open output file...\n" << pLsdft->OutFilename;
            cout << "Exiting for this configuration...\n";
            exit(1);
        }
    }


    if (at != pLsdft->Ntype) {
        cout << "Number of blocks: " << at << " in .ion file not same as number of atom type:" << pLsdft->
            Ntype << "in .inpt file\n";
        exit(1);
    }

    if (rank == 0) {
        fprintf(fConfFile, "Total number of atoms              :  %d\n", pLsdft->nAtoms);
        fprintf(fConfFile, "Total number of electrons          :  %d\n", pLsdft->noetot);
    }

    for (at = 0; at < pLsdft->Ntype; at++) {
        if (rank == 0) {
            fprintf(fConfFile, "Atom type %d                        :  %s\n", at + 1, pLsdft->atomType[at]);
            fprintf(fConfFile, "Number of atoms of type %d          :  %d\n", at + 1,
                (pLsdft->endPos[at] - pLsdft->startPos[at] + 1) / 3);

        }
    }

    if (rank == 0) {
        fclose(fConfFile);
    }
#ifdef _DEBUG
    VecView(pLsdft->Atompos, PETSC_VIEWER_STDOUT_SELF);
    cout << "noetot" << pLsdft->noetot << '\n';
    cout << "natoms" << pLsdft->nAtoms << '\n';
#endif

    return;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//         Read_Relax: reads the .relax file for constraint on movement of atoms             //
///////////////////////////////////////////////////////////////////////////////////////////////
void Read_Relax(LSDFT_OBJ *pLsdft)
{
    FILE *fConfFile;
    char ConfigFile[100] = "./";
    strcat(ConfigFile, pLsdft->file);
    strcat(ConfigFile, ".relax");
    PetscInt poscnt, index = 0;
    PetscScalar xct, yct, zct;
    char name[10];
    int at = 0;
    PetscInt noa;

    VecDuplicate(pLsdft->forces, &pLsdft->mvAtmConstraint);
    VecSet(pLsdft->mvAtmConstraint, 1);

    if (pLsdft->RelaxFlag == 1) {
        if ((fConfFile = fopen((const char *) ConfigFile, "rb")) == NULL) {
            cout << "Couldn't open config file... " << ConfigFile << '\n';
            cout << "Exiting for this configuration...\n";
            exit(1);
        }

        do {
            if (fscanf(fConfFile, "%s", name) == EOF)
                break;
            if (fscanf(fConfFile, "%d", &noa) == EOF)
                break;
            /*
             * check if .relax file is consistent with .ion file
             */
            if (strcmp(name, pLsdft->atomType[at]) != 0) {
                cout << "Incorrect atom type in .relax file. Give same type and same order as .ion file\n";
                exit(1);
            }
            for (poscnt = 0; poscnt < noa; poscnt++) {
                fscanf(fConfFile, "%lf", &xct);
                fscanf(fConfFile, "%lf", &yct);
                fscanf(fConfFile, "%lf", &zct);

                VecSetValues(pLsdft->mvAtmConstraint, 1, &index, &xct, INSERT_VALUES);
                index++;
                VecSetValues(pLsdft->mvAtmConstraint, 1, &index, &yct, INSERT_VALUES);
                index++;
                VecSetValues(pLsdft->mvAtmConstraint, 1, &index, &zct, INSERT_VALUES);
                index++;
            }
            at++;
        } while (!feof(fConfFile));
        fclose(fConfFile);

        if (at != pLsdft->Ntype) {
            cout << "Number of blocks: " << at << " in .relax file not same as number of atom type: " << pLsdft->
                Ntype << " in .inpt file\n";
            exit(1);
        }
    }

    return;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//                   Read_Pseudopotential: reads the pseudopotential file                     //
////////////////////////////////////////////////////////////////////////////////////////////////
void Read_Pseudopotential(LSDFT_OBJ *pLsdft)
{
    FILE *fRvsVFile;
    int at, i, l, lineCtr, count = 0;
    char str[60], RvsVFile[100];
    PetscScalar value, rc;

    /*
     * allocate memory for pseudopotential object in the structure
     */
    pLsdft->psd = (PSD_OBJ *) malloc(sizeof(PSD_OBJ) * pLsdft->Ntype);
    assert(pLsdft->psd != NULL);

    for (at = 0; at < pLsdft->Ntype; at++) {
        strcpy(RvsVFile, "./pseudopotential/");
        strcat(RvsVFile, pLsdft->psdName[at]);
        if ((fRvsVFile = fopen((const char *) RvsVFile, "rb")) == NULL) {
            printf("Couldn't open pseudopotential file... %s\n", RvsVFile);
            exit(1);
        }
        count = 0;
        /*
         * check if pseudopotential file is consistent with input atom type
         */
        fscanf(fRvsVFile, "%s", str);
        if (strcmp(str, pLsdft->atomType[at]) != 0) {
            printf("pseudopotential file does not match with input atom type %d\n", at);
            exit(1);
        }
        /*
         * read the pseudopotential file (must be in Troulier-Martins format)
         */
        do {
            fgets(str, 60, fRvsVFile);
        } while (strcmp(str, " Radial grid follows\n"));
        do {
            fscanf(fRvsVFile, "%s", str);
            count++;
        } while (strcmp(str, "Pseudopotential") != 0);
        fclose(fRvsVFile);
        count = count - 1;
        pLsdft->psd[at].size = count + 1;   // store the size of the pseudopotential list for future

        /*
         * allocate memory for arrays storing radial grid, pseudopotentials and
         * pseudowavefunctions
         */
        pLsdft->psd[at].RadialGrid = (PetscScalar *) malloc(sizeof(PetscScalar) * pLsdft->psd[at].size);
        assert(pLsdft->psd[at].RadialGrid != NULL);

        pLsdft->psd[at].V = (PetscScalar **) malloc(sizeof(PetscScalar *) * 4);
        assert(pLsdft->psd[at].V != NULL);

        pLsdft->psd[at].U = (PetscScalar **) malloc(sizeof(PetscScalar *) * 4);
        assert(pLsdft->psd[at].U != NULL);

        pLsdft->psd[at].rc = (PetscScalar *) malloc(sizeof(PetscScalar) * 4);
        assert(pLsdft->psd[at].rc != NULL);

        pLsdft->psd[at].uu = (PetscScalar *) malloc(sizeof(PetscScalar) * pLsdft->psd[at].size);
        assert(pLsdft->psd[at].uu != NULL);

        pLsdft->psd[at].Vloc = (PetscScalar *) malloc(sizeof(PetscScalar) * pLsdft->psd[at].size);
        assert(pLsdft->psd[at].Vloc != NULL);

        /*
         * allocate memory for the pseudopotentials and pseudowavefunctions for each quantum number
         */
        for (i = 0; i < 4; i++) {
            pLsdft->psd[at].V[i] = (PetscScalar *) malloc(sizeof(PetscScalar) * pLsdft->psd[at].size);
            assert(pLsdft->psd[at].V[i] != NULL);

            pLsdft->psd[at].U[i] = (PetscScalar *) malloc(sizeof(PetscScalar) * pLsdft->psd[at].size);
            assert(pLsdft->psd[at].U[i] != NULL);
        }

        /*
         * open file again and read the pseudopotentials and pseudo wave functions
         */
        if ((fRvsVFile = fopen((const char *) RvsVFile, "rb")) == NULL) {
            printf("Couldn't open pseudopotential file... %s \n", RvsVFile);
            exit(1);
        }

        do {
            fgets(str, 60, fRvsVFile);
        } while (strcmp(str, " Radial grid follows\n"));

        pLsdft->psd[at].RadialGrid[0] = 0.0;
        for (i = 1; i <= count; i++) {
            fscanf(fRvsVFile, "%lf", &value);
            pLsdft->psd[at].RadialGrid[i] = value;
        }
        lineCtr = 0;
        while (strcmp(str, " Pseudopotential follows (l on next line)\n")) {
            fgets(str, 60, fRvsVFile);
            lineCtr++;
        }

        /*
         * read pseudopotential
         */
        while (strcmp(str, " Pseudopotential follows (l on next line)\n") == 0) {
            fscanf(fRvsVFile, "%d", &l);
            if (l == 0) {
                pLsdft->psd[at].lmax = 0;
                for (i = 1; i <= count; i++) {
                    fscanf(fRvsVFile, "%lf", &value);
                    value = 0.5 * value / (pLsdft->psd[at].RadialGrid[i]);
                    pLsdft->psd[at].V[0][i] = value;

                    if (pLsdft->localPsd[at] == 0)
                        pLsdft->psd[at].Vloc[i] = value;
                }
                pLsdft->psd[at].V[0][0] = pLsdft->psd[at].V[0][1];
                if (pLsdft->localPsd[at] == 0)
                    pLsdft->psd[at].Vloc[0] = pLsdft->psd[at].V[0][1];
            }
            if (l == 1) {
                pLsdft->psd[at].lmax = 1;
                for (i = 1; i <= count; i++) {
                    fscanf(fRvsVFile, "%lf", &value);
                    value = 0.5 * value / (pLsdft->psd[at].RadialGrid[i]);
                    pLsdft->psd[at].V[1][i] = value;

                    if (pLsdft->localPsd[at] == 1)
                        pLsdft->psd[at].Vloc[i] = value;
                }
                pLsdft->psd[at].V[1][0] = pLsdft->psd[at].V[1][1];
                if (pLsdft->localPsd[at] == 1)
                    pLsdft->psd[at].Vloc[0] = pLsdft->psd[at].V[1][1];
            }
            if (l == 2) {
                pLsdft->psd[at].lmax = 2;
                for (i = 1; i <= count; i++) {
                    fscanf(fRvsVFile, "%lf", &value);
                    value = 0.5 * value / (pLsdft->psd[at].RadialGrid[i]);
                    pLsdft->psd[at].V[2][i] = value;

                    if (pLsdft->localPsd[at] == 2)
                        pLsdft->psd[at].Vloc[i] = value;
                }
                pLsdft->psd[at].V[2][0] = pLsdft->psd[at].V[2][1];
                if (pLsdft->localPsd[at] == 2)
                    pLsdft->psd[at].Vloc[0] = pLsdft->psd[at].V[2][1];
            }
            if (l == 3) {
                pLsdft->psd[at].lmax = 3;
                for (i = 1; i <= count; i++) {
                    fscanf(fRvsVFile, "%lf", &value);
                    value = 0.5 * value / (pLsdft->psd[at].RadialGrid[i]);
                    pLsdft->psd[at].V[3][i] = value;

                    if (pLsdft->localPsd[at] == 3)
                        pLsdft->psd[at].Vloc[i] = value;
                }
                pLsdft->psd[at].V[3][0] = pLsdft->psd[at].V[3][1];
                if (pLsdft->localPsd[at] == 3)
                    pLsdft->psd[at].Vloc[0] = pLsdft->psd[at].V[3][1];
            }

            for (i = 0; i < lineCtr; i++)
                fgets(str, 60, fRvsVFile);

        }
        /*
         * read until valence charge block is found
         */
        while (strcmp(str, " Valence charge follows\n")) {
            fgets(str, 60, fRvsVFile);
        }

        /*
         * read valence charge
         */
        while (strcmp(str, " Valence charge follows\n") == 0) {
            for (i = 1; i <= count; i++) {
                fscanf(fRvsVFile, "%lf", &value);
                value = value / (4 * M_PI * pLsdft->psd[at].RadialGrid[i] * pLsdft->psd[at].RadialGrid[i]);
                pLsdft->psd[at].uu[i] = value;
            }
            pLsdft->psd[at].uu[0] = 0.0;
            for (i = 0; i < lineCtr; i++)
                fgets(str, 60, fRvsVFile);

        }
        /*
         * read pseudowavefunction
         */
        while (strcmp(str, " Pseudo-wave-function follows (l, zelect, rc)\n") == 0) {
            fscanf(fRvsVFile, "%d", &l);
            fscanf(fRvsVFile, "%lf", &value);
            fscanf(fRvsVFile, "%lf", &rc);
            if (l == 0)         //s orbital
            {
                pLsdft->psd[at].rc[0] = rc;
                for (i = 1; i <= count; i++) {
                    fscanf(fRvsVFile, "%lf", &value);
                    value = value / pLsdft->psd[at].RadialGrid[i];
                    pLsdft->psd[at].U[0][i] = value;
                }
                pLsdft->psd[at].U[0][0] = pLsdft->psd[at].U[0][1];
            }
            if (l == 1)         //p orbital
            {
                pLsdft->psd[at].rc[1] = rc;
                for (i = 1; i <= count; i++) {
                    fscanf(fRvsVFile, "%lf", &value);
                    value = value / pLsdft->psd[at].RadialGrid[i];
                    pLsdft->psd[at].U[1][i] = value;
                }
                pLsdft->psd[at].U[1][0] = pLsdft->psd[at].U[1][1];
            }
            if (l == 2)         //d orbital
            {
                pLsdft->psd[at].rc[2] = rc;
                for (i = 1; i <= count; i++) {
                    fscanf(fRvsVFile, "%lf", &value);
                    value = value / pLsdft->psd[at].RadialGrid[i];
                    pLsdft->psd[at].U[2][i] = value;
                }
                pLsdft->psd[at].U[2][0] = pLsdft->psd[at].U[2][1];
            }
            if (l == 3)         //f orbital
            {
                pLsdft->psd[at].rc[3] = rc;
                for (i = 1; i <= count; i++) {
                    fscanf(fRvsVFile, "%lf", &value);
                    value = value / pLsdft->psd[at].RadialGrid[i];
                    pLsdft->psd[at].U[3][i] = value;
                }
                pLsdft->psd[at].U[3][0] = pLsdft->psd[at].U[3][1];
            }

            for (i = 0; i < lineCtr; i++) {
                if (feof(fRvsVFile))
                    break;
                fgets(str, 60, fRvsVFile);
            }
        }
        fclose(fRvsVFile);
    }
    return;
}

////////////////////////////////////////////
/// print electron density
////////////////////////////////////////////
void Print_ElecDens(LSDFT_OBJ *pLsdft)
{
    FILE *fDenFile;
    char str[100] = "./", DenFile[100];
    int rank, k, j, i;
    PetscInt rcor, thetacor, zcor, lrdim, lthetadim, lzdim;
    PetscScalar ***RhoArray;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    strcat(str, pLsdft->DenFilename);
    strcat(str, "_");
    sprintf(DenFile, "%s%d", str, rank);
    strcat(DenFile, ".den");

    DMDAGetCorners(pLsdft->da, &rcor, &thetacor, &zcor, &lrdim, &lthetadim, &lzdim);

    if ((fDenFile = fopen((const char *) DenFile, "w")) == NULL) {
        printf("Couldn't create density file for printing density...\n");
        printf("Exiting...\n");
        exit(1);
    }

    fprintf(fDenFile, "%d %d %d %d %d %d\n", rcor, thetacor, zcor, lrdim, lthetadim, lzdim);

    DMDAVecGetArray(pLsdft->da, pLsdft->elecDensRho, &RhoArray);
    for (k = zcor; k < zcor + lzdim; k++)
        for (j = thetacor; j < thetacor + lthetadim; j++)
            for (i = rcor; i < rcor + lrdim; i++) {
                fprintf(fDenFile, "%d %d %d %0.16lf\n", k, j, i, RhoArray[k][j][i]);
            }
    DMDAVecRestoreArray(pLsdft->da, pLsdft->elecDensRho, &RhoArray);
    fclose(fDenFile);

}

////////////////////////////////////////////
/// transfer electron density
////////////////////////////////////////////
void TransferFields(LSDFT_OBJ *pLsdft)
{
    
  int rank, k, j, i,LI,K,J,I;
    PetscInt rcor, thetacor, zcor, lrdim, lthetadim, lzdim;
    PetscScalar ***RhoArray, ***PhiArray,*pRho,*pPhi;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 
    DMDAGetCorners(pLsdft->da, &rcor, &thetacor, &zcor, &lrdim, &lthetadim, &lzdim);
    DMDAVecGetArray(pLsdft->da, pLsdft->elecDensRho, &RhoArray);
    DMDAVecGetArray(pLsdft->da, pLsdft->potentialPhi, &PhiArray);
    VecGetArray(pLsdft->BpotentialPhi[0],&pPhi);
    VecGetArray(pLsdft->BelecDensRho[0],&pRho);
    for (k = zcor; k < zcor + lzdim; k++)
        for (j = thetacor; j < thetacor + lthetadim; j++)
            for (i = rcor; i < rcor + lrdim; i++)
	      {
		K=k,J=j;I=i;
		// except the corner point
		if(k==pLsdft->numPoints_z-1)
		  {
		    K=0;	
		  }
		if(j==pLsdft->numPoints_y-1)
		  {
		    J=0;	
		  }
		if(i==pLsdft->numPoints_x-1)
		  {
		    I=0;	
		  }  

		  LI=K*pLsdft->BNx[0]*pLsdft->BNy[0]+J*pLsdft->BNx[0]+I;
		  RhoArray[k][j][i] = pRho[LI];
		  PhiArray[k][j][i]= pPhi[LI];
            }
    DMDAVecRestoreArray(pLsdft->da, pLsdft->elecDensRho,&RhoArray);
    DMDAVecRestoreArray(pLsdft->da, pLsdft->potentialPhi,&PhiArray);
    VecRestoreArray(pLsdft->BpotentialPhi[0],&pPhi);
    VecRestoreArray(pLsdft->BelecDensRho[0],&pRho);
   

}
////////////////////////////////////////////
/// transfer electron density
////////////////////////////////////////////
void DisplayPotential(LSDFT_OBJ *pLsdft)
{
    
  int rank, k, j, i,LI,K,J,I;
    PetscInt rcor, thetacor, zcor, lrdim, lthetadim, lzdim;
    PetscScalar ***RhoArray, ***PhiArray,*pRho,*pPhi,Diff;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 
    DMDAGetCorners(pLsdft->da, &rcor, &thetacor, &zcor, &lrdim, &lthetadim, &lzdim);
    //DMDAVecGetArray(pLsdft->da, pLsdft->elecDensRho, &RhoArray);
    DMDAVecGetArray(pLsdft->da, pLsdft->potentialPhi, &PhiArray);
     VecGetArray(pLsdft->BpotentialPhi[0],&pPhi);
     //VecGetArray(pLsdft->BelecDensRho[0],&pRho);
    for (k = zcor; k < zcor + lzdim; k++)
        for (j = thetacor; j < thetacor + lthetadim; j++)
            for (i = rcor; i < rcor + lrdim; i++)
	      {
		K=k,J=j;I=i;
		// except the corner point
		if(k==pLsdft->numPoints_z-1)
		  {
		    K=0;	
		  }
		if(j==pLsdft->numPoints_y-1)
		  {
		    J=0;	
		  }
		if(i==pLsdft->numPoints_x-1)
		  {
		    I=0;	
		  }  

		  LI=K*pLsdft->BNx[0]*pLsdft->BNy[0]+J*pLsdft->BNx[0]+I;
		  // RhoArray[k][j][i] = pRho[LI];
		  // PhiArray[k][j][i] = pPhi[LI];
		  Diff=fabs(PhiArray[k][j][i]-pPhi[LI]);
		  printf("SG PhiArray[%d][%d][%d]=%0.16lf, pPhi[%d]=%0.16lf  Diff=%0.16lf \n",k,j,i,PhiArray[k][j][i],LI,pPhi[LI],Diff);
		  if(Diff>=1e-5)
		    {
		      printf("BAM!\n");
		    }
            }
    // DMDAVecRestoreArray(pLsdft->da, pLsdft->elecDensRho,&RhoArray);
    DMDAVecRestoreArray(pLsdft->da, pLsdft->potentialPhi,&PhiArray);
    VecRestoreArray(pLsdft->BpotentialPhi[0],&pPhi);
    //VecRestoreArray(pLsdft->BelecDensRho[0],&pRho);
   

}

////////////////////////////////////////////
/// read electron density
////////////////////////////////////////////
void Read_ElecDens(LSDFT_OBJ *pLsdft)
{
    FILE *fDenFile;
    char str[100] = "./", DenFile[100];
    int rank, k, j, i, K, J, I;
    PetscInt rcor, thetacor, zcor, lrdim, lthetadim, lzdim;
    PetscInt Rcor, Thetacor, Zcor, Lrdim, Lthetadim, Lzdim;
    PetscScalar ***RhoArray, rhoval;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    strcat(str, pLsdft->ReadDenFilename);
    strcat(str, "_");
    sprintf(DenFile, "%s%d", str, rank);
    strcat(DenFile, ".den");

    DMDAGetCorners(pLsdft->da, &rcor, &thetacor, &zcor, &lrdim, &lthetadim, &lzdim);

    if ((fDenFile = fopen((const char *) DenFile, "r")) == NULL) {
        printf("Couldn't open density file for reading...\n");
        printf("Exiting...\n");
        exit(1);
    }

    fscanf(fDenFile, "%d %d %d %d %d %d", &Rcor, &Thetacor, &Zcor, &Lrdim, &Lthetadim, &Lzdim);

    if ((rcor != Rcor) || (thetacor != Thetacor) || (zcor != Zcor) || (lrdim != Lrdim) || (lthetadim != Lthetadim)
        || (lzdim != Lzdim)) {
        printf("rank=%d, Different domain decomposition used in density input file\n", rank);
        printf("rank=%d: This simulation: rcor=%d,thetacor=%d,zcor=%d,lrdim=%d,lthetadim=%d,lzdim=%d \n", rank, rcor,
            thetacor, zcor, lrdim, lthetadim, lzdim);
        printf("rank=%d: This simulation: Rcor=%d,Thetacor=%d,Zcor=%d,Lrdim=%d,Lthetadim=%d,Lzdim=%d \n", rank, Rcor,
            Thetacor, Zcor, Lrdim, Lthetadim, Lzdim);
        exit(1);

    }

    DMDAVecGetArray(pLsdft->da, pLsdft->elecDensRho, &RhoArray);
    for (k = zcor; k < zcor + lzdim; k++)
        for (j = thetacor; j < thetacor + lthetadim; j++)
            for (i = rcor; i < rcor + lrdim; i++) {
                fscanf(fDenFile, "%d %d %d %lf", &K, &J, &I, &rhoval);
                RhoArray[K][J][I] = rhoval;
                // printf("%d %d %d %0.16lf\n",K,J,I,RhoArray[K][J][I]);
            }
    DMDAVecRestoreArray(pLsdft->da, pLsdft->elecDensRho, &RhoArray);
    fclose(fDenFile);

}

////////////////////////////////////////////
/// print max min eigenvalue
////////////////////////////////////////////
void Print_MinMaxLambda(LSDFT_OBJ *pLsdft)
{
    FILE *fDenFile;
    char str[100] = "./", DenFile[100];
    int rank, k, j, i, LIpos;
    PetscInt rcor, thetacor, zcor, lrdim, lthetadim, lzdim;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    strcat(str, pLsdft->DenFilename);
    strcat(str, "_");
    sprintf(DenFile, "%s%d", str, rank);
    strcat(DenFile, ".lambdas");

    DMDAGetCorners(pLsdft->da, &rcor, &thetacor, &zcor, &lrdim, &lthetadim, &lzdim);

    if ((fDenFile = fopen((const char *) DenFile, "w")) == NULL) {
        printf("Couldn't create density file for printing density...\n");
        printf("Exiting...\n");
        exit(1);
    }

    fprintf(fDenFile, "%d %d %d %d %d %d\n", rcor, thetacor, zcor, lrdim, lthetadim, lzdim);
    for (k = zcor; k < zcor + lzdim; k++)
        for (j = thetacor; j < thetacor + lthetadim; j++)
            for (i = rcor; i < rcor + lrdim; i++) {
	       LIpos = (k - zcor) * lrdim * lthetadim + (j - thetacor) * lrdim + (i - rcor);
	       fprintf(fDenFile, "%d %0.16lf %0.16lf\n",LIpos,pLsdft->lambda_min[LIpos],pLsdft->lambda_max[LIpos]);
            }
    fclose(fDenFile);

}

////////////////////////////////////////////
/// read electron density
////////////////////////////////////////////
void Read_MinMaxLambda(LSDFT_OBJ *pLsdft)
{
    FILE *fDenFile;
    char str[100] = "./", DenFile[100];
    int rank, k, j, i, K, J, I,LIpos;
    PetscInt rcor, thetacor, zcor, lrdim, lthetadim, lzdim;
    PetscInt Rcor, Thetacor, Zcor, Lrdim, Lthetadim, Lzdim;
    PetscScalar ***RhoArray, minval,maxval;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    strcat(str, pLsdft->ReadDenFilename);
    strcat(str, "_");
    sprintf(DenFile, "%s%d", str, rank);
    strcat(DenFile, ".lambdas");

    DMDAGetCorners(pLsdft->da, &rcor, &thetacor, &zcor, &lrdim, &lthetadim, &lzdim);

    if ((fDenFile = fopen((const char *) DenFile, "r")) == NULL) {
        printf("Couldn't open density file for reading...\n");
        printf("Exiting...\n");
        exit(1);
    }

    fscanf(fDenFile, "%d %d %d %d %d %d", &Rcor, &Thetacor, &Zcor, &Lrdim, &Lthetadim, &Lzdim);

    if ((rcor != Rcor) || (thetacor != Thetacor) || (zcor != Zcor) || (lrdim != Lrdim) || (lthetadim != Lthetadim)
        || (lzdim != Lzdim)) {
        printf("rank=%d, Different domain decomposition used in density input file\n", rank);
        printf("rank=%d: This simulation: rcor=%d,thetacor=%d,zcor=%d,lrdim=%d,lthetadim=%d,lzdim=%d \n", rank, rcor,
            thetacor, zcor, lrdim, lthetadim, lzdim);
        printf("rank=%d: This simulation: Rcor=%d,Thetacor=%d,Zcor=%d,Lrdim=%d,Lthetadim=%d,Lzdim=%d \n", rank, Rcor,
            Thetacor, Zcor, Lrdim, Lthetadim, Lzdim);
        exit(1);

    }

    // DMDAVecGetArray(pLsdft->da, pLsdft->elecDensRho, &RhoArray);
    for (k = zcor; k < zcor + lzdim; k++)
        for (j = thetacor; j < thetacor + lthetadim; j++)
            for (i = rcor; i < rcor + lrdim; i++) {
	      LIpos = (k - zcor) * lrdim * lthetadim + (j - thetacor) * lrdim + (i - rcor);
	      //fprintf(fDenFile, "%d %0.16lf %0.16lf\n",LIpos,,pLsdft->lambda_max[LIpos]);
	       fscanf(fDenFile, "%d %lf %lf", &LIpos, &minval, &maxval);
	       pLsdft->lambda_min[LIpos] = minval;
	       pLsdft->lambda_max[LIpos] = maxval;
          
                // printf("%d %d %d %0.16lf\n",K,J,I,RhoArray[K][J][I]);
            }
    // DMDAVecRestoreArray(pLsdft->da, pLsdft->elecDensRho, &RhoArray);
    fclose(fDenFile);

}

