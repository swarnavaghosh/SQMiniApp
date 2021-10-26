#include "lsdft.h"
#include "petscsys.h"
#include <petsctime.h> 
#include <iostream>
//using namespace std;
int NonlocalProjectors(LSDFT_OBJ *pLsdft)
{
  int at, poscnt, l,l2, m, k, j, i, nz, nx, ny,LI,Nlmxyz;
  int lmax,lloc,Nlm,zstart,zend,ystart,yend,xstart,xend,K,J,I,lmctr;
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   
  pLsdft->ProjectorNonlocal = (NEWPROJECTOR_OBJ *) malloc(sizeof(NEWPROJECTOR_OBJ) * pLsdft->Ntype);
  assert(pLsdft->ProjectorNonlocal != NULL);
   
  for (at = 0; at < pLsdft->Ntype; at++)
    {
      printf("at=%d \n",at);
       
      // allocate chi
      pLsdft->ProjectorNonlocal[at].Chi=(double **) malloc(sizeof(double *) * (pLsdft->AtomOverlap_nonlocal[at].Natoms));
      assert(pLsdft->ProjectorNonlocal[at].Chi != NULL);
      // allocate l index
      pLsdft->ProjectorNonlocal[at].lId=(int **) malloc(sizeof(int *) * (pLsdft->AtomOverlap_nonlocal[at].Natoms));
      assert(pLsdft->ProjectorNonlocal[at].lId != NULL);
      // allocate m index
      pLsdft->ProjectorNonlocal[at].mId=(int **) malloc(sizeof(int *) * (pLsdft->AtomOverlap_nonlocal[at].Natoms));
      assert(pLsdft->ProjectorNonlocal[at].mId != NULL);
      // allocate z index
      pLsdft->ProjectorNonlocal[at].kId=(int **) malloc(sizeof(int *) * (pLsdft->AtomOverlap_nonlocal[at].Natoms));
      assert(pLsdft->ProjectorNonlocal[at].kId != NULL);
      // allocate y index
       pLsdft->ProjectorNonlocal[at].jId=(int **) malloc(sizeof(int *) * (pLsdft->AtomOverlap_nonlocal[at].Natoms));
      assert(pLsdft->ProjectorNonlocal[at].jId != NULL);
      // allocate x index
       pLsdft->ProjectorNonlocal[at].iId=(int **) malloc(sizeof(int *) * (pLsdft->AtomOverlap_nonlocal[at].Natoms));
      assert(pLsdft->ProjectorNonlocal[at].iId != NULL);

      // allocate Nlmxyz index
       pLsdft->ProjectorNonlocal[at].Ndim=(int *) malloc(sizeof(int) * (pLsdft->AtomOverlap_nonlocal[at].Natoms));
      assert(pLsdft->ProjectorNonlocal[at].Ndim != NULL);
      
      // define sizes
      for (poscnt = 0; poscnt < pLsdft->AtomOverlap_nonlocal[at].Natoms; poscnt++)
	{
	  printf("poscnt=%d \n",poscnt);
	  
	  // number of nx,ny,nz
	  nz = pLsdft->AtomOverlap_nonlocal[at].zend[poscnt] - pLsdft->AtomOverlap_nonlocal[at].zstart[poscnt] + 1;
	  ny = pLsdft->AtomOverlap_nonlocal[at].yend[poscnt] - pLsdft->AtomOverlap_nonlocal[at].ystart[poscnt] + 1;
	  nx = pLsdft->AtomOverlap_nonlocal[at].xend[poscnt] - pLsdft->AtomOverlap_nonlocal[at].xstart[poscnt] + 1;
	  // number of lm index, leaving local pseudopotential out
	  lmax=pLsdft->psd[at].lmax;
	  lloc=pLsdft->localPsd[at];
	  Nlm=lmax*lmax+2*(lmax-lloc);
	  Nlmxyz=Nlm*nz*ny*nx;
	  pLsdft->ProjectorNonlocal[at].Ndim[poscnt] = Nlmxyz;	  
	  pLsdft->ProjectorNonlocal[at].Chi[poscnt]=(double *)malloc(sizeof(double)*(Nlmxyz));
	  assert(pLsdft->ProjectorNonlocal[at].Chi[poscnt] != NULL);
	  pLsdft->ProjectorNonlocal[at].lId[poscnt]=(int *)malloc(sizeof(int)*(Nlmxyz));
	  assert(pLsdft->ProjectorNonlocal[at].lId[poscnt] != NULL);
	  pLsdft->ProjectorNonlocal[at].mId[poscnt]=(int *)malloc(sizeof(int)*(Nlmxyz));
	  assert(pLsdft->ProjectorNonlocal[at].mId[poscnt]!=NULL);
	  pLsdft->ProjectorNonlocal[at].kId[poscnt]=(int *)malloc(sizeof(int)*(Nlmxyz));
	  assert(pLsdft->ProjectorNonlocal[at].kId[poscnt]!=NULL);
	  pLsdft->ProjectorNonlocal[at].jId[poscnt]=(int *)malloc(sizeof(int)*(Nlmxyz));
	  assert(pLsdft->ProjectorNonlocal[at].jId[poscnt]!=NULL);
	  pLsdft->ProjectorNonlocal[at].iId[poscnt]=(int *)malloc(sizeof(int)*(Nlmxyz));
          assert(pLsdft->ProjectorNonlocal[at].iId[poscnt]!=NULL);
	  
	   xstart = pLsdft->AtomOverlap_nonlocal[at].xstart[poscnt];
            ystart = pLsdft->AtomOverlap_nonlocal[at].ystart[poscnt];
            zstart = pLsdft->AtomOverlap_nonlocal[at].zstart[poscnt];
            xend = pLsdft->AtomOverlap_nonlocal[at].xend[poscnt];
            yend = pLsdft->AtomOverlap_nonlocal[at].yend[poscnt];
            zend = pLsdft->AtomOverlap_nonlocal[at].zend[poscnt];
	  
	  // go over the nested for-loops and copy data
	  LI=0;
	  // loop over influencing finite difference nodes
	  for (k = zstart; k <= zend; k++)
	    {
    	      K = k - zstart;
    	      for (j = ystart; j <= yend; j++)
		{
		  J = j - ystart;
		  for (i = xstart; i <= xend; i++)
		    {
		      I = i - xstart;
		      //check if this point is within the region of the truncated hamiltonian 
		      // if((k<=Zend && k>=Zstart) && (j<=Yend && j>=Ystart) && (i<=Xend && i>=Xstart))
		      // 	{
		      lmctr=0;
			  for (l = 0; l <= pLsdft->psd[at].lmax; l++)
			    {
			      if (l != pLsdft->localPsd[at])
				{
				  l2 = 2 * l;
				  for (m = 0; m <= l2; m++)
				    {
				     pLsdft->ProjectorNonlocal[at].Chi[poscnt][LI]=pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m][K][J][I];
				     pLsdft->ProjectorNonlocal[at].lId[poscnt][LI]=l;
				     pLsdft->ProjectorNonlocal[at].mId[poscnt][LI]=m-l;
				     pLsdft->ProjectorNonlocal[at].kId[poscnt][LI]=k;
				     pLsdft->ProjectorNonlocal[at].jId[poscnt][LI]=j;
				     pLsdft->ProjectorNonlocal[at].iId[poscnt][LI]=i;
				     // if(poscnt==4 && LI==0)
				     //   {
				     // 	 printf("HERE! l=%d\n",l);
				     //   }
				     // if(rank==0)
				     //   {
				     // 	 printf("LI=%d, poscnt=%d, l=%d, m=%d, k=%d, j=%d, i=%d \n",LI,poscnt,l,m,k,j,i);
				     //   }
				     LI++; // increment linear index
				     lmctr++;
				    } // for m
				} // if l!
			    }//for l
			  //}// if k<=Zend
		    }// for i
		}// for j
	    }// for k
	  // printf("at=%d, poscnt=%d, LI=%d, Nlmxyz=%d, nx=%d, ny=%d, nz=%d, Nlm=%d, lmctr=%d\n",at,poscnt,LI,Nlmxyz,nx,ny,nz,Nlm,lmctr);
	  
	}// for poscnt
    }// for at

  return 0;
}







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
