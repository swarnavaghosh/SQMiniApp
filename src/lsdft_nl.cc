#include "lsdft.h"
#include "petscsys.h"
#include <petsctime.h> 
#include <iostream>
//using namespace std;
///////////////////////////////////////////////////////////////////////////////////////////////
//                 CalculateNonlocalProjectors: Calculate nonlocal projectors                //
///////////////////////////////////////////////////////////////////////////////////////////////
PetscScalar CalculateNonlocalProjectors(LSDFT_OBJ *pLsdft)
{
    PetscInt at, poscnt, l, m, k, j, i, nz, nx, ny;
    PetscScalar x, y, z, r, zz, yy, SpHarmonic, pUlDeltaVl;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("ok1\n");
    pLsdft->Projector_nonlocal = (PROJECTOR_OBJ *) malloc(sizeof(PROJECTOR_OBJ) * pLsdft->Ntype);
    assert(pLsdft->Projector_nonlocal != NULL);
    printf("ok2\n");
    for (at = 0; at < pLsdft->Ntype; at++) {
        pLsdft->Projector_nonlocal[at].Chi =
            (double ******) malloc(sizeof(double *****) * (pLsdft->AtomOverlap_nonlocal[at].Natoms));
        assert(pLsdft->Projector_nonlocal[at].Chi != NULL);
	printf("ok3\n");
        for (poscnt = 0; poscnt < pLsdft->AtomOverlap_nonlocal[at].Natoms; poscnt++) {
            nz = pLsdft->AtomOverlap_nonlocal[at].zend[poscnt] - pLsdft->AtomOverlap_nonlocal[at].zstart[poscnt] + 1;
            ny = pLsdft->AtomOverlap_nonlocal[at].yend[poscnt] - pLsdft->AtomOverlap_nonlocal[at].ystart[poscnt] + 1;
            nx = pLsdft->AtomOverlap_nonlocal[at].xend[poscnt] - pLsdft->AtomOverlap_nonlocal[at].xstart[poscnt] + 1;

            pLsdft->Projector_nonlocal[at].Chi[poscnt] =
                (double *****) malloc(sizeof(double ****) * (pLsdft->psd[at].lmax + 1));
            assert(pLsdft->Projector_nonlocal[at].Chi[poscnt] != NULL);
            for (l = 0; l <= pLsdft->psd[at].lmax; l++) {
                if (l != pLsdft->localPsd[at]) {
                    pLsdft->Projector_nonlocal[at].Chi[poscnt][l] =
                        (double ****) malloc(sizeof(double ***) * (2 * l + 1));
                    assert(pLsdft->Projector_nonlocal[at].Chi[poscnt][l] != NULL);
                    for (m = -l; m <= l; m++) {
                        pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m + l] =
                            (double ***) malloc(sizeof(double **) * nz);
                        assert(pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m + l] != NULL);
                        for (k = 0; k < nz; k++) {
                            z = pLsdft->delta_z * (k + pLsdft->AtomOverlap_nonlocal[at].zstart[poscnt]) -
                                pLsdft->range_z - pLsdft->AtomOverlap_nonlocal[at].Z0[poscnt];
                            zz = z * z;
                            pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m + l][k] =
                                (double **) malloc(sizeof(double *) * ny);
                            assert(pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m + l][k] != NULL);

                            for (j = 0; j < ny; j++) {
                                y = pLsdft->delta_y * (j + pLsdft->AtomOverlap_nonlocal[at].ystart[poscnt]) -
                                    pLsdft->range_y - pLsdft->AtomOverlap_nonlocal[at].Y0[poscnt];
                                yy = y * y;
                                pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m + l][k][j] =
                                    (double *) malloc(sizeof(double) * nx);
                                assert(pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m + l][k][j] != NULL);
                                for (i = 0; i < nx; i++) {
                                    x = pLsdft->delta_x * (i + pLsdft->AtomOverlap_nonlocal[at].xstart[poscnt]) -
                                        pLsdft->range_x - pLsdft->AtomOverlap_nonlocal[at].X0[poscnt];

				    //printf("SS: i=%d, X=%lf, X0=%lf,x=%lf \n",i,pLsdft->delta_x*(i + pLsdft->AtomOverlap_nonlocal[at].xstart[poscnt])-pLsdft->range_x,pLsdft->AtomOverlap_nonlocal[at].X0[poscnt],x);
                                    r = sqrt(x * x + yy + zz);
                                    SplineEvaluate(pLsdft->psd[at].RadialGrid, pLsdft->psd[at].UDeltaV[l],
                                        pLsdft->psd[at].size, &r, &pUlDeltaVl, 1, pLsdft->psd[at].SplineFitUDeltaV[l]);
                                    SpHarmonic = RealSphericalHarmonic(x, y, z, l, m);
                                    pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m + l][k][j][i] =
                                        pUlDeltaVl * SpHarmonic;
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
