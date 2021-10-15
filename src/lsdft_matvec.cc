
/*=============================================================================================
  |
  | file name: matvec.cc
  |
  | Description: Mat Vec routines
  | Authors: Swarnava Ghosh
  |
  |
  |-------------------------------------------------------------------------------------------*/
#include "lsdft.h"
#include "petscsys.h"
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define SIGN(a,b) ((b) >= 0.0 ? PetscAbsScalar(a) : -PetscAbsScalar(a)) 

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
//               Mult_HamiltonianVector: function to multiply Hamiltonian and orbitals       //
///////////////////////////////////////////////////////////////////////////////////////////////
PetscErrorCode Mult_HamiltonianVector(LSDFT_OBJ *pLsdft, Vec *V1, Vec *V2,int Kp, int Jp, int Ip)
{
  PetscScalar ***Varray, ***V2array, *alpha;
  PetscInt i, j, k, l, m, n, p, l2, I, J, K, at, poscnt, xstart, ystart, zstart, xend, yend, zend, shift, rowStart,rowEnd, Nrows, Index, posIdx;
  PetscInt xcor,ycor,zcor,lxdim,lydim,lzdim;
    PetscErrorCode ierr;
    // AO aodmda1;
    int rank;
    PetscInt Xstart,Ystart,Zstart,Xend,Yend,Zend; // starting and ending nodes of the truncated hamiltonian

    VecZeroEntries(pLsdft->VeffVktemp);

    Xstart=Ip-ceil(pLsdft->Rcut/pLsdft->delta_x);
    Xend=Ip+ceil(pLsdft->Rcut/pLsdft->delta_x);

    Ystart=Jp-ceil(pLsdft->Rcut/pLsdft->delta_y);
    Yend=Jp+ceil(pLsdft->Rcut/pLsdft->delta_y);

    Zstart=Kp-ceil(pLsdft->Rcut/pLsdft->delta_z);
    Zend=Kp+ceil(pLsdft->Rcut/pLsdft->delta_z);

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // DMDAGetCorners(pLsdft->daloccell, &xstartcell, &ystartcell, &zstartcell, &lxdimcell, &lydimcell,&lzdimcell);
    DMDAGetCorners(pLsdft->da, &xcor, &ycor, &zcor, &lxdim, &lydim,&lzdim);
   
     alpha = (PetscScalar *) calloc(pLsdft->IP_length, sizeof(PetscScalar));
     
    // // printf("MV2 mult here at rank=%d",rank);

     // (Lap+Veff)*V1=V2                                                                                                                                                              
     MatMult(pLsdft->LapPlusVeffOprloc,*V1,*V2);

     
     DMDAVecGetArray(pLsdft->daloc, *V1, &Varray);
    DMDAVecGetArray(pLsdft->daloc, *V2, &V2array);

    

    // temp comment from here

    for (at = 0; at < pLsdft->Ntype; at++) {
        
        // loop over influencing atoms
        for (poscnt = 0; poscnt < pLsdft->AtomOverlap_nonlocal[at].Natoms; poscnt++) 
    	  {
            xstart = pLsdft->AtomOverlap_nonlocal[at].xstart[poscnt];
            ystart = pLsdft->AtomOverlap_nonlocal[at].ystart[poscnt];
            zstart = pLsdft->AtomOverlap_nonlocal[at].zstart[poscnt];
            xend = pLsdft->AtomOverlap_nonlocal[at].xend[poscnt];
            yend = pLsdft->AtomOverlap_nonlocal[at].yend[poscnt];
            zend = pLsdft->AtomOverlap_nonlocal[at].zend[poscnt];
            p = (int) (pLsdft->AtomOverlap_nonlocal[at].xindex[poscnt] / 3);

            
             // loop over influencing finite difference nodes
            for (k = zstart; k <= zend; k++) {
    	      K = k - zstart;
    	      for (j = ystart; j <= yend; j++) {
    		J = j - ystart;
    		for (i = xstart; i <= xend; i++) {
    		  I = i - xstart;
    		  //check if this point is within the region of the truncated hamiltonian 
    		  if((k<=Zend && k>=Zstart) && (j<=Yend && j>=Ystart) && (i<=Xend && i>=Xstart))
    		    {
    		      posIdx = pLsdft->IP_index[p];
		        // printf("rank=%d,posIdx=%d, Varray[%d][%d][%d]=%0.16lf \n",rank,posIdx,k-Zstart,j-Ystart,i-Xstart,Varray[k-Zstart][j-Ystart][i-Xstart]);
    		      
    		       // loop over angular momentums
    		      for (l = 0; l <= pLsdft->psd[at].lmax; l++) {
    			if (l != pLsdft->localPsd[at]) {
    			  l2 = 2 * l;
    			  for (m = 0; m <= l2; m++) {
    			    alpha[posIdx] +=
    			      Varray[k-Zstart][j-Ystart][i-Xstart]*pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m][K][J][I] *
    			      pLsdft->psd[at].Gamma[l];
    			    posIdx++;
    			  }
    			}
    		      }
    		    }

    		}
    	      }
            }
    	  }
    }

   

    // printf("At here1 \n");
    
     // compute the action of nonlocal operator on the orbitals
    
    for (at = 0; at < pLsdft->Ntype; at++) {
        
         // loop over influencing atoms
        for (poscnt = 0; poscnt < pLsdft->AtomOverlap_nonlocal[at].Natoms; poscnt++) {
            xstart = pLsdft->AtomOverlap_nonlocal[at].xstart[poscnt];
            ystart = pLsdft->AtomOverlap_nonlocal[at].ystart[poscnt];
            zstart = pLsdft->AtomOverlap_nonlocal[at].zstart[poscnt];
            xend = pLsdft->AtomOverlap_nonlocal[at].xend[poscnt];
            yend = pLsdft->AtomOverlap_nonlocal[at].yend[poscnt];
            zend = pLsdft->AtomOverlap_nonlocal[at].zend[poscnt];
            p = (int) (pLsdft->AtomOverlap_nonlocal[at].xindex[poscnt] / 3);
            
             // loop over influencing finite difference nodes
            for (k = zstart; k <= zend; k++) {
                K = k - zstart;
                for (j = ystart; j <= yend; j++) {
                    J = j - ystart;
                    for (i = xstart; i <= xend; i++) {
                        I = i - xstart;

    				//check if this point is within the region of the truncated hamiltonian 
    			if((k<=Zend && k>=Zstart) && (j<=Yend && j>=Ystart) && (i<=Xend && i>=Xstart))
    			  {

                        posIdx = pLsdft->IP_index[p];
                        
                         // loop over angular momentums
                        for (l = 0; l <= pLsdft->psd[at].lmax; l++) {
                            if (l != pLsdft->localPsd[at]) {
                                l2 = 2 * l;
                                for (m = 0; m <= l2; m++) {
                                    V2array[k-Zstart][j-Ystart][i-Xstart] +=
                                        pLsdft->Projector_nonlocal[at].Chi[poscnt][l][m][K][J][I] * alpha[posIdx];
                                    posIdx++;
                                }
                            }
                        }
    			  }
                    }
                }
            }
        }
    }

   
    DMDAVecRestoreArray(pLsdft->daloc, *V1, &Varray);
    DMDAVecRestoreArray(pLsdft->daloc, *V2, &V2array);
 
    free(alpha);
    return ierr;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
