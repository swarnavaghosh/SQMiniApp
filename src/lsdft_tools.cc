/*=============================================================================================
  |
  | file name: tools.cc
  |
  | Authors: Swarnava Ghosh
  |
  |
  |-------------------------------------------------------------------------------------------*/
#include "lsdft.h"
#include "petscsys.h"
#define SIGN(a,b) ((b) >= 0.0 ? PetscAbsScalar(a) : -PetscAbsScalar(a))
////////////////////////////////////////////////////////////////////////////////////////////////
//                SplineEvaluate: cubic spline evaluation from precalculated data             //
////////////////////////////////////////////////////////////////////////////////////////////////
void SplineEvaluate(PetscScalar *X1, PetscScalar *Y1, int len1, PetscScalar *X2, PetscScalar *Y2, int len2,
    PetscScalar *YD)
{
    int i, j;
    PetscScalar A0, A1, A2, A3, x, dx, dy, p1, p2, p3;
    if (X2[0] < X1[0] || X2[len2 - 1] > X1[len1 - 1]) {
        printf
            ("First input X in table=%lf, interpolate at x[first]=%lf, last input X in table=%lf, interpolate at x[last]=%lf\n",
            X1[0], X2[0], X1[len1 - 1], X2[len2 - 1]);
        printf("out of range in spline interpolation\n");
        exit(1);
    }

    /*
     * p1 is left endpoint of the interval
     * p2 is resampling position
     * p3 is right endpoint of interval
     * j is input index of current interval
     */
    p3 = X2[0] - 1;
    for (i = j = 0; i < len2; i++) {
        /*
         * check if in new interval
         */
        p2 = X2[i];
        if (p2 > p3) {
            /*
             * find interval which contains p2
             */
            for (; j < len1 && p2 > X1[j]; j++);
            if (p2 < X1[j])
                j--;
            p1 = X1[j];
            p3 = X1[j + 1];

            /*
             * coefficients
             */
            dx = 1.0 / (X1[j + 1] - X1[j]);
            dy = (Y1[j + 1] - Y1[j]) * dx;
            A0 = Y1[j];
            A1 = YD[j];
            A2 = dx * (3.0 * dy - 2.0 * YD[j] - YD[j + 1]);
            A3 = dx * dx * (-2.0 * dy + YD[j] + YD[j + 1]);
        }
        /*
         * use Horner's rule to calculate cubic polynomial
         */
        x = p2 - p1;
        Y2[i] = ((A3 * x + A2) * x + A1) * x + A0;
    }

    return;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//                  getYD_gen: derivatives required for spline interpolation                  //
//      Reference: Cubic Spline Interpolation: A Review, George Wolberg, Technical Report,    //
//                      Department of Computer Science, Columbia University                   //
////////////////////////////////////////////////////////////////////////////////////////////////
void getYD_gen(PetscScalar *X, PetscScalar *Y, PetscScalar *YD, int len)
{
    int i;
    PetscScalar h0, h1, r0, r1, *A, *B, *C;

    PetscMalloc(sizeof(PetscScalar) * len, &A);
    PetscMalloc(sizeof(PetscScalar) * len, &B);
    PetscMalloc(sizeof(PetscScalar) * len, &C);


    h0 = X[1] - X[0];
    h1 = X[2] - X[1];
    r0 = (Y[1] - Y[0]) / h0;
    r1 = (Y[2] - Y[1]) / h1;
    B[0] = h1 * (h0 + h1);
    C[0] = (h0 + h1) * (h0 + h1);
    YD[0] = r0 * (3 * h0 * h1 + 2 * h1 * h1) + r1 * h0 * h0;

    for (i = 1; i < len - 1; i++) {
        h0 = X[i] - X[i - 1];
        h1 = X[i + 1] - X[i];
        r0 = (Y[i] - Y[i - 1]) / h0;
        r1 = (Y[i + 1] - Y[i]) / h1;
        A[i] = h1;
        B[i] = 2 * (h0 + h1);
        C[i] = h0;
        YD[i] = 3 * (r0 * h1 + r1 * h0);
    }

    A[i] = (h0 + h1) * (h0 + h1);
    B[i] = h0 * (h0 + h1);
    YD[i] = r0 * h1 * h1 + r1 * (3 * h0 * h1 + 2 * h0 * h0);

    tridiag_gen(A, B, C, YD, len);

    PetscFree(A);
    PetscFree(B);
    PetscFree(C);

    return;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//                  tridiag_gen: Solves a tridiagonal system using Gauss Elimination          //
//      Reference: Cubic Spline Interpolation: A Review, George Wolberg, Technical Report,    //
//                      Department of Computer Science, Columbia University                   //
////////////////////////////////////////////////////////////////////////////////////////////////
void tridiag_gen(PetscScalar *A, PetscScalar *B, PetscScalar *C, PetscScalar *D, int len)
{
    int i;
    PetscScalar b, *F;
    PetscMalloc(sizeof(PetscScalar) * len, &F);

    /*
     * Gauss elimination; forward substitution
     */
    b = B[0];
    D[0] = D[0] / b;
    for (i = 1; i < len; i++) {
        F[i] = C[i - 1] / b;
        b = B[i] - A[i] * F[i];
        if (b == 0) {
            printf("Divide by zero in tridiag_gen\n");
            exit(1);
        }
        D[i] = (D[i] - D[i - 1] * A[i]) / b;
    }

    /*
     * backsubstitution
     */
    for (i = len - 2; i >= 0; i--)
        D[i] -= (D[i + 1] * F[i + 1]);

    PetscFree(F);
    return;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//       RealSphericalHarmonic: returns the real spherical harmonic for some given l,m       //
///////////////////////////////////////////////////////////////////////////////////////////////
PetscScalar RealSphericalHarmonic(PetscScalar x, PetscScalar y, PetscScalar z, int l, int m)
{
    /*
     * only l=0,1,2,3,4,5,6 implemented for now
     */

    /*
     * l=0
     */
    PetscScalar C00 = 0.282094791773878;    // 0.5*sqrt(1/pi)
    /*
     * l=1
     */
    PetscScalar C1m1 = 0.488602511902920;   // sqrt(3/(4*pi))
    PetscScalar C10 = 0.488602511902920;    // sqrt(3/(4*pi))
    PetscScalar C1p1 = 0.488602511902920;   // sqrt(3/(4*pi))
    /*
     * l=2
     */
    PetscScalar C2m2 = 1.092548430592079;   // 0.5*sqrt(15/pi)
    PetscScalar C2m1 = 1.092548430592079;   // 0.5*sqrt(15/pi)
    PetscScalar C20 = 0.315391565252520;    // 0.25*sqrt(5/pi)
    PetscScalar C2p1 = 1.092548430592079;   // 0.5*sqrt(15/pi)
    PetscScalar C2p2 = 0.546274215296040;   // 0.25*sqrt(15/pi)
    /*
     * l=3
     */
    PetscScalar C3m3 = 0.590043589926644;   // 0.25*sqrt(35/(2*pi))
    PetscScalar C3m2 = 2.890611442640554;   // 0.5*sqrt(105/(pi))
    PetscScalar C3m1 = 0.457045799464466;   // 0.25*sqrt(21/(2*pi))
    PetscScalar C30 = 0.373176332590115;    // 0.25*sqrt(7/pi)
    PetscScalar C3p1 = 0.457045799464466;   // 0.25*sqrt(21/(2*pi))
    PetscScalar C3p2 = 1.445305721320277;   //  0.25*sqrt(105/(pi))
    PetscScalar C3p3 = 0.590043589926644;   //  0.25*sqrt(35/(2*pi))

    PetscScalar pi = M_PI;
    PetscScalar p;
    PetscScalar r = sqrt(x * x + y * y + z * z);
    PetscScalar SH = 0.0;

    if (l == 0)
        SH = C00;

    if (r >= 1e-10) {
        if (l == 0)
            SH = C00;

        else if (l == 1) {
            if (m == -1)
                SH = C1m1 * (y / r);
            else if (m == 0)
                SH = C10 * (z / r);
            else if (m == 1)
                SH = C1p1 * (x / r);
            else {
                printf("incorrect l: %d,m: %d\n", l, m);
                exit(1);
            }
        } else if (l == 2) {
            if (m == -2)
                SH = C2m2 * (x * y) / (r * r);
            else if (m == -1)
                SH = C2m1 * (y * z) / (r * r);
            else if (m == 0) {
                SH = C20 * (-x * x - y * y + 2.0 * z * z) / (r * r);
            } else if (m == 1)
                SH = C2p1 * (z * x) / (r * r);
            else if (m == 2)
                SH = C2p2 * (x * x - y * y) / (r * r);
            else {
                printf("incorrect l: %d,m: %d\n", l, m);
                exit(1);
            }
        } else if (l == 3) {
            if (m == -3)
                SH = C3m3 * (3 * x * x - y * y) * y / (r * r * r);
            else if (m == -2)
                SH = C3m2 * (x * y * z) / (r * r * r);
            else if (m == -1)
                SH = C3m1 * y * (4 * z * z - x * x - y * y) / (r * r * r);
            else if (m == 0)
                SH = C30 * z * (2 * z * z - 3 * x * x - 3 * y * y) / (r * r * r);
            else if (m == 1)
                SH = C3p1 * x * (4 * z * z - x * x - y * y) / (r * r * r);
            else if (m == 2)
                SH = C3p2 * z * (x * x - y * y) / (r * r * r);
            else if (m == 3)
                SH = C3p3 * x * (x * x - 3 * y * y) / (r * r * r);
            else {
                printf("incorrect l: %d,m: %d\n", l, m);
                exit(1);
            }
        } else if (l == 4) {
            if (m == -4)
                SH = (3.0 / 4.0) * sqrt(35.0 / pi) * (x * y * (x * x - y * y)) / (r * r * r * r);
            else if (m == -3)
                SH = (3.0 / 4.0) * sqrt(35.0 / (2.0 * pi)) * (3.0 * x * x - y * y) * y * z / (r * r * r * r);
            else if (m == -2)
                SH = (3.0 / 4.0) * sqrt(5.0 / pi) * x * y * (7.0 * z * z - r * r) / (r * r * r * r);
            else if (m == -1)
                SH = (3.0 / 4.0) * sqrt(5.0 / (2.0 * pi)) * y * z * (7.0 * z * z - 3.0 * r * r) / (r * r * r * r);
            else if (m == 0)
                SH = (3.0 / 16.0) * sqrt(1.0 / pi) * (35.0 * z * z * z * z - 30.0 * z * z * r * r +
                    3.0 * r * r * r * r) / (r * r * r * r);
            else if (m == 1)
                SH = (3.0 / 4.0) * sqrt(5.0 / (2.0 * pi)) * x * z * (7.0 * z * z - 3.0 * r * r) / (r * r * r * r);
            else if (m == 2)
                SH = (3.0 / 8.0) * sqrt(5.0 / (pi)) * (x * x - y * y) * (7.0 * z * z - r * r) / (r * r * r * r);
            else if (m == 3)
                SH = (3.0 / 4.0) * sqrt(35.0 / (2.0 * pi)) * (x * x - 3.0 * y * y) * x * z / (r * r * r * r);
            else if (m == 4)
                SH = (3.0 / 16.0) * sqrt(35.0 / pi) * (x * x * (x * x - 3.0 * y * y) - y * y * (3.0 * x * x -
                        y * y)) / (r * r * r * r);
            else {
                printf("incorrect l: %d,m: %d\n", l, m);
                exit(1);
            }
        } else if (l == 5) {
            p = sqrt(x * x + y * y);
            if (m == -5)
                SH = (3.0 * sqrt(2.0 * 77.0 / pi) / 32.0) * (8.0 * x * x * x * x * y - 4.0 * x * x * y * y * y +
                    4.0 * pow(y, 5) - 3.0 * y * p * p * p * p) / (r * r * r * r * r);
            else if (m == -4)
                SH = (3.0 / 16.0) * sqrt(385.0 / pi) * (4.0 * x * x * x * y -
                    4.0 * x * y * y * y) * z / (r * r * r * r * r);
            else if (m == -3)
                SH = (sqrt(2.0 * 385.0 / pi) / 32.0) * (3.0 * y * p * p - 4.0 * y * y * y) * (9 * z * z -
                    r * r) / (r * r * r * r * r);
            else if (m == -2)
                SH = (1.0 / 8.0) * sqrt(1155.0 / pi) * 2.0 * x * y * (3.0 * z * z * z -
                    z * r * r) / (r * r * r * r * r);
            else if (m == -1)
                SH = (1.0 / 16.0) * sqrt(165.0 / pi) * y * (21.0 * z * z * z * z - 14.0 * r * r * z * z +
                    r * r * r * r) / (r * r * r * r * r);
            else if (m == 0)
                SH = (1.0 / 16.0) * sqrt(11.0 / pi) * (63.0 * z * z * z * z * z - 70.0 * z * z * z * r * r +
                    15.0 * z * r * r * r * r) / (r * r * r * r * r);
            else if (m == 1)
                SH = (1.0 / 16.0) * sqrt(165.0 / pi) * x * (21.0 * z * z * z * z - 14.0 * r * r * z * z +
                    r * r * r * r) / (r * r * r * r * r);
            else if (m == 2)
                SH = (1.0 / 8.0) * sqrt(1155.0 / pi) * (x * x - y * y) * (3.0 * z * z * z -
                    r * r * z) / (r * r * r * r * r);
            else if (m == 3)
                SH = (sqrt(2.0 * 385.0 / pi) / 32.0) * (4.0 * x * x * x - 3.0 * p * p * x) * (9.0 * z * z -
                    r * r) / (r * r * r * r * r);
            else if (m == 4)
                SH = (3.0 / 16.0) * sqrt(385.0 / pi) * (4.0 * (x * x * x * x + y * y * y * y) -
                    3.0 * p * p * p * p) * z / (r * r * r * r * r);
            else if (m == 5)
                SH = (3.0 * sqrt(2.0) / 32.0) * sqrt(77.0 / pi) * (4.0 * x * x * x * x * x + 8.0 * x * y * y * y * y -
                    4.0 * x * x * x * y * y - 3.0 * x * p * p * p * p) / (r * r * r * r * r);
            else {
                printf("incorrect l: %d,m: %d\n", l, m);
                exit(1);
            }
        } else if (l == 6) {
            p = sqrt(x * x + y * y);
            if (m == -6)
                SH = (sqrt(2.0 * 3003.0 / pi) / 64.0) * (12.0 * pow(x, 5) * y + 12.0 * x * pow(y,
                        5) - 8.0 * x * x * x * y * y * y - 6.0 * x * y * pow(p, 4)) / (r * r * r * r * r * r);
            else if (m == -5)
                SH = (3.0 / 32.0) * sqrt(2.0 * 1001.0 / pi) * (8.0 * pow(x,
                        4) * y - 4.0 * x * x * y * y * y + 4.0 * pow(y, 5) - 3.0 * y * pow(p,
                        4)) * z / (r * r * r * r * r * r);
            else if (m == -4)
                SH = (3.0 / 32.0) * sqrt(91.0 / pi) * (4.0 * x * x * x * y - 4.0 * x * y * y * y) * (11.0 * z * z -
                    r * r) / (r * r * r * r * r * r);
            else if (m == -3)
                SH = (sqrt(2.0 * 1365.0 / pi) / 32.0) * (-4.0 * y * y * y + 3.0 * y * p * p) * (11.0 * z * z * z -
                    3.0 * z * r * r) / (r * r * r * r * r * r);
            else if (m == -2)
                SH = (sqrt(2.0 * 1365 / pi) / 64.0) * (2.0 * x * y) * (33.0 * pow(z, 4) - 18.0 * z * z * r * r + pow(r,
                        4)) / (r * r * r * r * r * r);
            else if (m == -1)
                SH = (sqrt(273.0 / pi) / 16.0) * y * (33.0 * pow(z, 5) - 30.0 * z * z * z * r * r + 5.0 * z * pow(r,
                        4)) / (r * r * r * r * r * r);
            else if (m == 0)
                SH = (sqrt(13.0 / pi) / 32.0) * (231.0 * pow(z, 6) - 315 * pow(z, 4) * r * r + 105.0 * z * z * pow(r,
                        4) - 5.0 * pow(r, 6)) / (r * r * r * r * r * r);
            else if (m == 1)
                SH = (sqrt(273.0 / pi) / 16.0) * x * (33.0 * pow(z, 5) - 30.0 * z * z * z * r * r + 5 * z * pow(r,
                        4)) / (r * r * r * r * r * r);
            else if (m == 2)
                SH = (sqrt(2.0 * 1365 / pi) / 64.0) * (x * x - y * y) * (33.0 * pow(z,
                        4) - 18.0 * z * z * r * r + pow(r, 4)) / (r * r * r * r * r * r);
            else if (m == 3)
                SH = (sqrt(2.0 * 1365.0 / pi) / 32.0) * (4.0 * x * x * x - 3.0 * x * p * p) * (11.0 * z * z * z -
                    3.0 * z * r * r) / (r * r * r * r * r * r);
            else if (m == 4)
                SH = (3.0 / 32.0) * sqrt(91.0 / pi) * (4.0 * pow(x, 4) + 4.0 * pow(y, 4) - 3.0 * pow(p,
                        4)) * (11.0 * z * z - r * r) / (r * r * r * r * r * r);
            else if (m == 5)
                SH = (3.0 / 32.0) * sqrt(2.0 * 1001.0 / pi) * (4.0 * pow(x, 5) + 8.0 * x * pow(y,
                        4) - 4.0 * x * x * x * y * y - 3.0 * x * pow(p, 4)) * z / (r * r * r * r * r * r);
            else if (m == 6)
                SH = (sqrt(2.0 * 3003.0 / pi) / 64.0) * (4.0 * pow(x, 6) - 4.0 * pow(y, 6) + 12.0 * x * x * pow(y,
                        4) - 12.0 * pow(x, 4) * y * y + 3.0 * y * y * pow(p, 4) - 3.0 * x * x * pow(p,
                        4)) / (r * r * r * r * r * r);
            else {
                printf("incorrect l: %d,m: %d\n", l, m);
                exit(1);
            }
        } else {
            printf("Only l=0,1,2,3,4,5,6 supported. Input l:%d\n", l);
            exit(1);
        }

    }

    return SH;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//                ComplexSphericalHarmonic: returns the complex spherical harmonics          //
///////////////////////////////////////////////////////////////////////////////////////////////
PetscScalar ComplexSphericalHarmonic(PetscScalar x, PetscScalar y, PetscScalar z, int l, int m, PetscScalar *SHreal,
    PetscScalar *SHimag)
{
    if (m < 0) {
        *SHreal = (1.0 / sqrt(2.0)) * RealSphericalHarmonic(x, y, z, l, abs(m));
        *SHimag = -(1.0 / sqrt(2.0)) * RealSphericalHarmonic(x, y, z, l, -abs(m));
    } else if (m == 0) {
        *SHreal = RealSphericalHarmonic(x, y, z, l, 0);
        *SHimag = 0.0;
    } else if (m > 0) {
        *SHreal = (pow(-1.0, m) / sqrt(2.0)) * RealSphericalHarmonic(x, y, z, l, abs(m));
        *SHimag = (pow(-1.0, m) / sqrt(2.0)) * RealSphericalHarmonic(x, y, z, l, -abs(m));
    } else {
        printf("incorrect m=%d in complex spherical harmonic \n", m);
    }

    return 0;

}

///////////////////////////////////////////////////////////////////////////////////////////////
//                 Mult_ComplexScalar: returns the product of two complex scalars            //
///////////////////////////////////////////////////////////////////////////////////////////////
PetscErrorCode Mult_ComplexScalar(PetscScalar A1, PetscScalar A2, PetscScalar B1, PetscScalar B2, PetscScalar *C1,
    PetscScalar *C2)
{
    *C1 = A1 * B1 - A2 * B2;
    *C2 = A1 * B2 + A2 * B1;
    return 0;
}

// ///////////////////////////////////////////////////////////////////////////////////////////////
// //                 Mult_ComplexMatMat: returns the product of two complex matrices           //
// ///////////////////////////////////////////////////////////////////////////////////////////////
// PetscErrorCode Mult_ComplexMatMat(LSDFT_OBJ* pLsdft,Mat *A1,Mat *A2,Mat *B1,Mat *B2,Mat *C1,Mat *C2)
// {
//   /*
//    * Real part
//    */
//   MatMatMultNumeric(*A1,*B1,*C1);
//   MatMatMultNumeric(*A2,*B2,pLsdft->ZOrb1);
//   MatAXPY(*C1,-1.0,pLsdft->ZOrb1,DIFFERENT_NONZERO_PATTERN);

//   /*
//    * Imaginary part
//    */
//   MatMatMultNumeric(*A1,*B2,*C2);
//   MatMatMultNumeric(*A2,*B1,pLsdft->ZOrb1);
//   MatAXPY(*C2,1.0,pLsdft->ZOrb1,DIFFERENT_NONZERO_PATTERN);

//   return 0;
// }
///////////////////////////////////////////////////////////////////////////////////////////////
//                 Mult_ComplexMatVec: returns the product of complex matrix vector           //
///////////////////////////////////////////////////////////////////////////////////////////////
PetscErrorCode Mult_ComplexMatVec(LSDFT_OBJ *pLsdft, Mat *A1, Mat *A2, Vec *V1, Vec *V2, Vec *C1, Vec *C2)
{
    /*
     * Real part
     */
    MatMult(*A1, *V1, *C1);
    MatMult(*A2, *V2, pLsdft->tempVec);
    VecAXPY(*C1, -1.0, pLsdft->tempVec);

    /*
     * Imaginary part
     */
    MatMult(*A1, *V2, *C2);
    MatMult(*A2, *V1, pLsdft->tempVec);
    VecAXPY(*C2, 1.0, pLsdft->tempVec);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//                              VecNorm_Complex: norm of a complex vector                     //
////////////////////////////////////////////////////////////////////////////////////////////////
PetscErrorCode VecNorm_Complex(Vec *Vr, Vec *Vi, PetscScalar *norm)
{
    PetscScalar rNorm, iNorm;

    VecNorm(*Vr, NORM_2, &rNorm);
    VecNorm(*Vi, NORM_2, &iNorm);

    *norm = sqrt(rNorm * rNorm + iNorm * iNorm);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//                        VecDot_Complex: dot product of two complex vectors                  //
////////////////////////////////////////////////////////////////////////////////////////////////
PetscErrorCode VecDot_Complex(Vec *V1r, Vec *V1i, Vec *V2r, Vec *V2i, PetscScalar *rvecdot, PetscScalar *ivecdot)
{
    PetscScalar rDot, iDot;

    VecDot(*V1r, *V2r, &rDot);
    VecDot(*V1i, *V2i, &iDot);
    *rvecdot = rDot + iDot;

    VecDot(*V1i, *V2r, &rDot);
    VecDot(*V1r, *V2i, &iDot);
    *ivecdot = -rDot + iDot;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//                        VecScale_Complex: scale a complex vector by a scalar                //
////////////////////////////////////////////////////////////////////////////////////////////////
PetscErrorCode VecScale_Complex(Vec *Vr, Vec *Vi, PetscScalar scale)
{

    VecScale(*Vr, scale);
    VecScale(*Vi, scale);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//                        ScalarNorm_Complex: norm of a complex scalar                        //
////////////////////////////////////////////////////////////////////////////////////////////////
PetscScalar ScalarNorm_Complex(PetscScalar ar, PetscScalar ai)
{
    PetscScalar norm;
    norm = sqrt(ar * ar + ai * ai);
    return norm;

}

///////////////////////////////////////////////////////////////////////////////////////////////
//  Exp_ComplexScalar: returns the exponential of i times dot product of two complex numbers //
//                                 Eulers formula                                            //
///////////////////////////////////////////////////////////////////////////////////////////////
PetscErrorCode Exp_ComplexScalar(PetscScalar k1, PetscScalar k2, PetscScalar k3, PetscScalar v1, PetscScalar v2,
    PetscScalar v3, PetscScalar *C1, PetscScalar *C2)
{
    PetscScalar kv = k1 * v1 + k2 * v2 + k3 * v3;
    *C1 = cos(kv);
    *C2 = sin(kv);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//               fract: calculates (n!)^2/((n-k)!(n+k)!), used for calculating the           //
//                    finite-difference weights of the gradient and lplacian                 //
///////////////////////////////////////////////////////////////////////////////////////////////
PetscScalar fract(PetscInt n, PetscInt k)
{
    int i;
    PetscScalar Nr = 1, Dr = 1, val;

    for (i = n - k + 1; i <= n; i++)
        Nr *= i;
    for (i = n + 1; i <= n + k; i++)
        Dr *= i;
    val = Nr / Dr;

    return (val);
}

//////////////////////////////////////////////////////////////////////////////////////////////
//         Lanczos:  Lanczos algorithm for calculating the minimum and maximum              //
//                                  eigenvalues of Hamiltonian                              //
//        see W.H. Press, Numerical recepies 3rd edition: The art of scientific computing,  //
//                                 Cambridge university press, 2007                         //
/////////////////////////////////////////////////////////////////////////////////////////////
// void Lanczos(LSDFT_OBJ* pLsdft,PetscScalar* EigenMin,PetscScalar* EigenMax,int kpt)
// {

//   PetscScalar tolLanczos=pLsdft->LANCZOSTOL;
//   PetscScalar Vnorm,norm;
//   PetscInt Nx,Ny,Nz;
//   PetscScalar Lk,Mk,Lkp1,Mkp1,deltaL,deltaM;
//   int k;
//   Vec Vk,iVk;
//   Vec Vkm1,iVkm1;
//   Vec Vkp1,iVkp1;
//   PetscScalar *a,*b,*ia;

//   Nx = pLsdft->numPoints_x;
//   Ny = pLsdft->numPoints_y;
//   Nz = pLsdft->numPoints_z;

//   PetscMalloc(sizeof(PetscScalar)*(Nx*Ny*Nz),&a);
//   PetscMalloc(sizeof(PetscScalar)*(Nx*Ny*Nz),&b);

//   VecDuplicate(pLsdft->elecDensRho,&Vk);
//   VecDuplicate(pLsdft->elecDensRho,&Vkm1);
//   VecDuplicate(pLsdft->elecDensRho,&Vkp1);

//   // if(pLsdft->Nkpts==1)
//   //   {
//       VecSet(Vkm1,1.0);
//       VecNorm(Vkm1,NORM_2,&Vnorm);
//       VecScale(Vkm1,1.0/Vnorm);
//       Mult_HamiltonianVector(pLsdft,&Vkm1,&Vk);
//       VecDot(Vkm1,Vk,&a[0]);
//       VecAXPY(Vk,-a[0],Vkm1);
//       VecNorm(Vk,NORM_2,&b[0]);
//       VecScale(Vk,1.0/b[0]);

//     // }else
//     // {
//     //   PetscMalloc(sizeof(PetscScalar)*(Nx*Ny*Nz),&ia);
//     //    VecDuplicate(pLsdft->elecDensRho,&iVk);
//     //    VecDuplicate(pLsdft->elecDensRho,&iVkm1);
//     //    VecDuplicate(pLsdft->elecDensRho,&iVkp1);

//     //    VecSet(Vkm1,1.0);
//     //    VecSet(iVkm1,1.0);
//     //    VecNorm_Complex(&Vkm1,&iVkm1,&Vnorm);
//     //    VecScale_Complex(&Vkm1,&iVkm1,1.0/Vnorm);

//     //    Mult_HamiltonianVector_Complex(pLsdft,pLsdft->k1[kpt],pLsdft->k2[kpt],pLsdft->k3[kpt],&Vkm1,&iVkm1,&Vk,&iVk);
//     //    VecDot_Complex(&Vkm1,&iVkm1,&Vk,&iVk,&a[0],&ia[0]);
//     //    VecAXPY(Vk,-a[0],Vkm1);
//     //    VecAXPY(Vk,ia[0],iVkm1);

//     //    VecAXPY(iVk,-a[0],iVkm1);
//     //    VecAXPY(iVk,-ia[0],Vkm1);

//     //    VecNorm_Complex(&Vk,&iVk,&b[0]);
//     //    VecScale_Complex(&Vk,&iVk,1.0/b[0]);

//     // }

//   k=0;
//   Lk=0.0;
//   Mk=0.0;
//   deltaL=1.0;
//   deltaM=1.0;

//   while((deltaL>tolLanczos || deltaM>tolLanczos))
//     {
//       // if(pLsdft->Nkpts == 1)
//       //     {
//    Mult_HamiltonianVector(pLsdft,&Vk,&Vkp1);
//    VecDot(Vk,Vkp1,&a[k+1]);
//    VecAXPY(Vkp1,-a[k+1],Vk);
//    VecAXPY(Vkp1,-b[k],Vkm1);
//    VecCopy(Vk,Vkm1);
//    VecNorm(Vkp1,NORM_2,&b[k+1]);
//    VecCopy(Vkp1,Vk);
//    VecScale(Vk,1.0/b[k+1]);

//  // }else
//  // {

//  //   Mult_HamiltonianVector_Complex(pLsdft,pLsdft->k1[kpt],pLsdft->k2[kpt],pLsdft->k3[kpt],&Vk,&iVk,&Vkp1,&iVkp1);
//  //   VecDot_Complex(&Vk,&iVk,&Vkp1,&iVkp1,&a[k+1],&ia[k+1]);

//  //   VecAXPY(Vkp1,-a[k+1],Vk);
//  //   VecAXPY(Vkp1,ia[k+1],iVk);

//  //   VecAXPY(iVkp1,-a[k+1],iVk);
//  //   VecAXPY(iVkp1,-ia[k+1],Vk);

//  //   VecAXPY(Vkp1,-b[k],Vkm1);
//  //   VecAXPY(iVkp1,-b[k],iVkm1);

//  //   VecCopy(Vk,Vkm1);
//  //   VecCopy(iVk,iVkm1);

//  //   VecNorm_Complex(&Vkp1,&iVkp1,&b[k+1]);
//  //   VecCopy(Vkp1,Vk);
//  //   VecCopy(iVkp1,iVk);
//  //   VecScale_Complex(&Vk,&iVk,1.0/b[k+1]);
//  // }

//       /*
//        * Call function to find eigenvalue of Tridiagonal matrix here minimum eigenvalue is
//        * Lkp1, maximum eigenvalue is Mkp1
//        */

//       TridiagEigenSolve(a,b,k+2,&Lkp1,&Mkp1);
//       deltaL = PetscAbsScalar(Lkp1-Lk);
//       deltaM = PetscAbsScalar(Mkp1-Mk);
//       Lk=Lkp1;
//       Mk=Mkp1;
//       k++;
//     }
//   *EigenMin = Lkp1;
//   *EigenMax = Mkp1;

//   PetscFree(a);
//   PetscFree(b);
//   VecDestroy(&Vk);
//   VecDestroy(&Vkm1);
//   VecDestroy(&Vkp1);
//   // if(pLsdft->Nkpts >1)
//   //   {
//   //    PetscFree(ia);
//   //    VecDestroy(&iVk);
//   //    VecDestroy(&iVkm1);
//   //    VecDestroy(&iVkp1);
//   //   }

//   return;
// }
// //////////////////////////////////////////////////////////////////////////////////////////////
// //   kPointLanczos:  Lanczos algorithm for calculating the minimum and maximum              //
// //                                eigenvalues of k-point Hamiltonian                        //
// //        see W.H. Press, Numerical recepies 3rd edition: The art of scientific computing,  //
// //                                 Cambridge university press, 2007                         //
// /////////////////////////////////////////////////////////////////////////////////////////////
// void kPointLanczos(LSDFT_OBJ* pLsdft,PetscScalar* EigenMin,PetscScalar* EigenMax,int kpt)
// {

//   PetscScalar tolLanczos=pLsdft->LANCZOSTOL;
//   PetscScalar rVnorm,iVnorm;
//   PetscInt Nx,Ny,Nz;
//   PetscScalar Lk,Mk,Lkp1,Mkp1,deltaL,deltaM;
//   PetscScalar temp;
//   int k;

//   Nx = pLsdft->numPoints_x;
//   Ny = pLsdft->numPoints_y;
//   Nz = pLsdft->numPoints_z;

//   PetscScalar *a,*b;
//   PetscMalloc(sizeof(PetscScalar)*(Nx*Ny*Nz),&a);
//   PetscMalloc(sizeof(PetscScalar)*(Nx*Ny*Nz),&b);

//   PetscScalar c1,c2;
//   Mult_ComplexScalar(1,2,2,-3,&c1,&c2);
//   Exp_ComplexScalar(1,2,5,3,4,8,&c1,&c2);

//   /*
//    * real
//    */
//   Vec rVk;
//   Vec rVkm1;
//   Vec rVkp1;

//   /*
//    * imaginary
//    */
//   Vec iVk;
//   Vec iVkm1;
//   Vec iVkp1;
//   Vec tempVec;

//   VecDuplicate(pLsdft->elecDensRho,&rVk);
//   VecDuplicate(pLsdft->elecDensRho,&rVkm1);
//   VecDuplicate(pLsdft->elecDensRho,&rVkp1);

//   VecDuplicate(pLsdft->elecDensRho,&iVk);
//   VecDuplicate(pLsdft->elecDensRho,&iVkm1);
//   VecDuplicate(pLsdft->elecDensRho,&iVkp1);

//   VecDuplicate(pLsdft->elecDensRho,&tempVec);

//   VecSet(rVkm1,1.0);
//   VecSet(iVkm1,1.0);

//   VecNorm(rVkm1,NORM_2,&rVnorm);
//   VecNorm(iVkm1,NORM_2,&iVnorm);
//   VecScale(rVkm1,1.0/(sqrt(rVnorm*rVnorm + iVnorm*iVnorm)));
//   VecScale(iVkm1,1.0/(sqrt(rVnorm*rVnorm + iVnorm*iVnorm)));

//    Mult_HamiltonianVector_Complex(pLsdft,pLsdft->k1[kpt],pLsdft->k2[kpt],pLsdft->k3[kpt],&rVkm1,&iVkm1,&rVk,&iVk);

//   VecDot(rVkm1,rVk,&a[0]);
//   VecDot(iVkm1,iVk,&temp);
//   a[0]=a[0]+temp;

//   VecAXPY(rVk,-a[0],rVkm1);
//   VecAXPY(iVk,-a[0],iVkm1);

//   VecNorm(rVk,NORM_2,&b[0]);
//   VecNorm(iVk,NORM_2,&temp);
//   b[0] = sqrt(b[0]*b[0]+ temp*temp);

//   VecScale(rVk,1.0/b[0]);
//   VecScale(iVk,1.0/b[0]);

//   k=0;
//   Lk=0.0;
//   Mk=0.0;
//   deltaL=1.0;
//   deltaM=1.0;

//   while((deltaL>tolLanczos || deltaM>tolLanczos))
//     {

//       Mult_HamiltonianVector_Complex(pLsdft,pLsdft->k1[kpt],pLsdft->k2[kpt],pLsdft->k3[kpt],&rVk,&iVk,&rVkp1,&iVkp1);

//       VecDot(rVk,rVkp1,&a[k+1]);
//       VecDot(iVk,iVkp1,&temp);
//       a[k+1]=a[k+1]+temp;

//       VecAXPY(rVkp1,-a[k+1],rVk); // Vkp1 = Vkp1 - ak+1 Vk
//       VecAXPY(iVkp1,-a[k+1],iVk); // Vkp1 = Vkp1 - ak+1 Vk

//       VecAXPY(rVkp1,-b[k],rVkm1); // Vkp1 = Vkp1 -bk Vkm1
//       VecAXPY(iVkp1,-b[k],iVkm1); // Vkp1 = Vkp1 -bk Vkm1

//       VecCopy(rVk,rVkm1); // Vkm1=Vk
//       VecCopy(iVk,iVkm1); // Vkm1=Vk

//       VecNorm(rVkp1,NORM_2,&b[k+1]);
//       VecNorm(iVkp1,NORM_2,&temp);
//       b[k+1] = sqrt(b[k+1]*b[k+1]+ temp*temp);

//       VecCopy(rVkp1,rVk);
//       VecCopy(iVkp1,iVk);

//       VecScale(rVk,1.0/b[k+1]); // Vk=Vkp1/b[k+1]
//       VecScale(iVk,1.0/b[k+1]); // Vk=Vkp1/b[k+1]

//       /*
//        * Call function to find eigenvalue of Tridiagonal matrix here minimum eigenvalue is
//        * Lkp1, maximum eigenvalue is Mkp1
//        */

//       TridiagEigenSolve(a,b,k+2,&Lkp1,&Mkp1);

//       deltaL = PetscAbsScalar(Lkp1-Lk);
//       deltaM = PetscAbsScalar(Mkp1-Mk);
//       Lk=Lkp1;
//       Mk=Mkp1;
//       k++;

//     }
//   *EigenMin = Lkp1;
//   *EigenMax = Mkp1;

//   PetscFree(a);
//   PetscFree(b);
//   VecDestroy(&rVk);
//   VecDestroy(&rVkm1);
//   VecDestroy(&rVkp1);

//   VecDestroy(&iVk);
//   VecDestroy(&iVkm1);
//   VecDestroy(&iVkp1);
//   VecDestroy(&tempVec);

//   return;
// }
///////////////////////////////////////////////////////////////////////////////////////////////
//       TridiagEigenSolve: Tridiagonal eigen solver for calculating the eigenvalues of      //
//                                         tridiagonal matrix.                               //
//       see W.H. Press, Numerical recepies 3rd edition: The art of scientific computing,    //
//                                  Cambridge university press, 2007                         //
///////////////////////////////////////////////////////////////////////////////////////////////
void TridiagEigenSolve(PetscScalar diag[], PetscScalar subdiag[], int n, PetscScalar *EigenMin, PetscScalar *EigenMax)
{


    int m, l, iter, i;
    PetscScalar s, r, p, g, f, dd, c, b;

    PetscScalar *d, *e;
    PetscMalloc(sizeof(PetscScalar) * n, &d);
    PetscMalloc(sizeof(PetscScalar) * n, &e);
    /*
     * create copy of diag and subdiag in d and e
     */
    for (i = 0; i < n; i++) {
        d[i] = diag[i];
        e[i] = subdiag[i];
    }

    /*
     * e has the subdiagonal elements, ignore last element(n-1) of e by making it zero
     */
    e[n - 1] = 0.0;
    for (l = 0; l <= n - 1; l++) {
        iter = 0;
        do {
            for (m = l; m <= n - 2; m++) {
                dd = PetscAbsScalar(d[m]) + PetscAbsScalar(d[m + 1]);
                if ((PetscScalar) (PetscAbsScalar(e[m]) + dd) == dd)
                    break;
            }
            if (m != l) {
                if (iter++ == 50) {
                    PetscPrintf(PETSC_COMM_SELF, "Too many iterations in Tridiagonal solver\n");
                    exit(1);
                }
                g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                r = sqrt(g * g + 1.0);
                g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
                s = c = 1.0;
                p = 0.0;

                for (i = m - 1; i >= l; i--) {
                    f = s * e[i];
                    b = c * e[i];
                    e[i + 1] = (r = sqrt(g * g + f * f));
                    if (r == 0.0) {
                        d[i + 1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    d[i + 1] = g + (p = s * r);
                    g = c * r - b;

                }
                if (r == 0.0 && i >= l)
                    continue;
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m != l);
    }

    /*
     * go over the array d to find the smallest and largest eigenvalue
     */
    *EigenMin = d[0];
    *EigenMax = d[0];

    for (i = 1; i < n; i++) {
        if (d[i] > *EigenMax) {
            *EigenMax = d[i];
        } else if (d[i] < *EigenMin) {
            *EigenMin = d[i];
        }
    }

    PetscFree(d);
    PetscFree(e);

    return;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//       TridiagEigenVecSolve: Tridiagonal eigen solver for calculating the eigenvalues of      //
//                                         tridiagonal matrix.                               //
//       see W.H. Press, Numerical recepies 3rd edition: The art of scientific computing,    //
//                                  Cambridge university press, 2007                         //
///////////////////////////////////////////////////////////////////////////////////////////////
void TridiagEigenVecSolve_NodesAndWeights(LSDFT_OBJ *pLsdft, PetscScalar diag[], PetscScalar subdiag[], int n, int LIp)
{
    int m, l, iter, i, k;
    PetscScalar s, r, p, g, f, dd, c, b, EigenMin, EigenMax;

    PetscScalar *d, *e, **z;
    PetscMalloc(sizeof(PetscScalar) * n, &d);
    PetscMalloc(sizeof(PetscScalar) * n, &e);
    // allocate memory for z
    PetscMalloc(sizeof(PetscScalar *) * n, &z);
    /*
     * create copy of diag and subdiag in d and e
     */
    for (i = 0; i < n; i++) {
        PetscMalloc(sizeof(PetscScalar) * n, &z[i]);
        d[i] = diag[i];
        e[i] = subdiag[i];
        // z as identity
        for (k = 0; k < n; k++) {
            z[i][k] = 0.0;
        }
        z[i][i] = 1.0;
        //  printf("d[%d]=%lf, e[%d]=%lf \n",i,d[i],i,e[i]);
    }

    /*
     * e has the subdiagonal elements, ignore last element(n-1) of e by making it zero
     */
    e[n - 1] = 0.0;
    for (l = 0; l <= n - 1; l++) {
        iter = 0;
        do {
            for (m = l; m <= n - 2; m++) {
                dd = PetscAbsScalar(d[m]) + PetscAbsScalar(d[m + 1]);
                if ((PetscScalar) (PetscAbsScalar(e[m]) + dd) == dd)
                    break;
            }
            if (m != l) {
                if (iter++ == 50) {
                    PetscPrintf(PETSC_COMM_SELF, "Too many iterations in Tridiagonal solver\n");
                    exit(1);
                }
                g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                r = sqrt(g * g + 1.0);
                g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
                s = c = 1.0;
                p = 0.0;

                for (i = m - 1; i >= l; i--) {
                    f = s * e[i];
                    b = c * e[i];
                    e[i + 1] = (r = sqrt(g * g + f * f));
                    if (r == 0.0) {
                        d[i + 1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    d[i + 1] = g + (p = s * r);
                    g = c * r - b;
                    /*calculation of eigenvectors */
                    for (k = 0; k < n; k++) {
                        f = z[k][i + 1];
                        z[k][i + 1] = s * z[k][i] + c * f;
                        z[k][i] = c * z[k][i] - s * f;
                    }
                }
                if (r == 0.0 && i >= l)
                    continue;
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m != l);
    }

     /*
   * go over the array d to find the smallest and largest eigenvalue
   */
  EigenMin = d[0];
  EigenMax = d[0];
    /*
     * go over the array to find nodes and weights and min and max eigenvalues
     */
    for (i = 0; i < n; i++) 
      {
        pLsdft->lssgq_lambda[LIp][i] = d[i];
        pLsdft->lssgq_weights[LIp][i] = z[0][i] * z[0][i];
	if(d[i]>EigenMax)
	{
	  EigenMax =d[i];
	}
      else if(d[i]<EigenMin)
	{
	  EigenMin =d[i];
	}          

      }
    
    // store eigenmin and eigenmax
    pLsdft->lambda_max[LIp] = EigenMax;
    pLsdft->lambda_min[LIp] = EigenMin;
    
    PetscFree(d);
    PetscFree(e);
    // need to free z
    for (i = 0; i < n; i++) {
        PetscFree(z[i]);
    }
    PetscFree(z);
    return;
}
///////////////////////////////////////////////////////////////////////////////////////

PetscScalar TrilinearInterpolation(PetscScalar x0,PetscScalar y0,PetscScalar z0,PetscScalar x1,PetscScalar y1,PetscScalar z1,PetscScalar c000,PetscScalar c001,PetscScalar c010,PetscScalar c011,PetscScalar c100,PetscScalar c101,PetscScalar c110,PetscScalar c111,PetscScalar x,PetscScalar y,PetscScalar z)
{

  PetscScalar a0,a1,a2,a3,a4,a5,a6,a7,val;
  PetscScalar D,c0,c1,c00,c01,c10,c11,A00,A10,A01,A11;


   if (fabs(x0-x1) <=1e-4 && fabs(y0-y1) <=1e-4 && fabs(z0-z1)>1e-3)
    {
      // 1d interpolation in z
      c0=c000;
      c1=c001;
      val = c0 + (z-z0)*(c1-c0)/(z1-z0);
  
    }
   else if (fabs(x0-x1) <=1e-4 && fabs(z0-z1) <=1e-4 && fabs(y0-y1)>1e-3)
    {
      // 1d interpolation in y
      c0=c000;
      c1=c010;
      val = c0 + (y-y0)*(c1-c0)/(y1-y0);

    } 
   else if (fabs(y0-y1) <=1e-4 && fabs(z0-z1) <=1e-4 && fabs(x0-x1)>1e-3)
    {
      // 1d interpolation in x
      c0=c000;
      c1=c100;
      val = c0 + (x-x0)*(c1-c0)/(x1-x0);
  
    }else if (fabs(x0-x1) <=1e-4 && fabs(y0-y1) >1e-3 && fabs(z0-z1)>1e-3)
    {
      // Bilinear interpolation in yz
      D=1.0/((y1-y0)*(z1-z0));
      
      c00=c000;
      c01=c001;
      c10=c010;
      c11=c011;
      
      A00=(y1-y)*(z1-z);
      A10=(y-y0)*(z1-z);
      A01=(y1-y)*(z-z0);
      A11=(y-y0)*(z-z0);
      
     val= D*(c00*A00 + c10*A10 + c01*A01 + c11*A11);
      
    } else if (fabs(z0-z1) <=1e-4 && fabs(y0-y1) >1e-3 && fabs(x0-x1)>1e-3)
    {
      // Bilinear interpolation in xy
       D=1.0/((x1-x0)*(y1-y0));
      
      c00=c000;
      c01=c010;
      c10=c100;
      c11=c110;
      
      A00=(x1-x)*(y1-y);
      A10=(x-x0)*(y1-y);
      A01=(x1-x)*(y-y0);
      A11=(x-x0)*(y-y0);
      
     val= D*(c00*A00 + c10*A10 + c01*A01 + c11*A11);

    } else if (fabs(y0-y1) <=1e-4 && fabs(z0-z1) >1e-3 && fabs(x0-x1)>1e-3)
    {
      // Bilinear interpolation in xz
       D=1.0/((x1-x0)*(z1-z0));
      
      c00=c000;
      c01=c001;
      c10=c100;
      c11=c101;
      
      A00=(x1-x)*(z1-z);
      A10=(x-x0)*(z1-z);
      A01=(x1-x)*(z-z0);
      A11=(x-x0)*(z-z0);
      
     val= D*(c00*A00 + c10*A10 + c01*A01 + c11*A11);
  
    } else if(fabs(x0-x1) >1e-3 && fabs(y0-y1) >1e-3 && fabs(z0-z1) >1e-3)
    {
      // trilinear interpolation in xyz
      D=1.0/((x0-x1)*(y0-y1)*(z0-z1));
  a0 = (-c000*x1*y1*z1 + c001*x1*y1*z0 + c010*x1*y0*z1 - c011*x1*y0*z0 + c100*x0*y1*z1 - c101*x0*y1*z0 - c110*x0*y0*z1 + c111*x0*y0*z0)*D;
 
  a1 = (+c000*y1*z1 - c001*y1*z0 - c010*y0*z1 + c011*y0*z0 - c100*y1*z1 + c101*y1*z0 + c110*y0*z1 - c111*y0*z0)*D; 
  a2 = (+c000*x1*z1 - c001*x1*z0 - c010*x1*z1 + c011*x1*z0 - c100*x0*z1 + c101*x0*z0 + c110*x0*z1 - c111*x0*z0)*D;
  a3 = (+c000*x1*y1 - c001*x1*y1 - c010*x1*y0 + c011*x1*y0 - c100*x0*y1 + c101*x0*y1 + c110*x0*y0 - c111*x0*y0)*D;

  a4 = (-c000*z1 + c001*z0 + c010*z1 - c011*z0 + c100*z1 - c101*z0 - c110*z1 + c111*z0)*D;
  a5 = (-c000*y1 + c001*y1 + c010*y0 - c011*y0 + c100*y1 - c101*y1 - c110*y0 + c111*y0)*D;
  a6 = (-c000*x1 + c001*x1 + c010*x1 - c011*x1 + c100*x0 - c101*x0 - c110*x0 + c111*x0)*D;

  a7 = (+c000 - c001 - c010 + c011 - c100 + c101 + c110 - c111)*D;

  val = a0 + a1*x + a2*y + a3*z + a4*x*y + a5*x*z + a6*y*z + a7*x*y*z;

    }else if (fabs(x0-x1) <=1e-4 && fabs(y0-y1) <=1e-4 && fabs(z0-z1) <=1e-4)
     {

       val=c000;
     }

  return val;

}
// PetscScalar NearestNeighborInterpolation(PetscScalar x0,PetscScalar y0,PetscScalar z0,PetscScalar x1,PetscScalar y1,PetscScalar z1,PetscScalar c000,PetscScalar c001,PetscScalar c010,PetscScalar c011,PetscScalar c100,PetscScalar c101,PetscScalar c110,PetscScalar c111,PetscScalar x,PetscScalar y,PetscScalar z)
// {




//   return val;
// }
///////////////////////////////////////////////////////////////////////////////////////
void FindMeshCorners(PetscScalar x, PetscScalar y,PetscScalar z,PetscScalar Rx,PetscScalar Ry,PetscScalar Rz, PetscScalar rx,PetscScalar ry,PetscScalar rz, int nx,int ny,int nz,int *i0,int *j0,int *k0,int *i1,int *j1,int *k1,PetscScalar *X, PetscScalar *Y,PetscScalar *Z)
{

  // Rx,Ry,Rz is the simulation domain
  // rx,ry,rz is the boundary domain
  // nx,ny,nz is number of nodes in boundary domain

  int ii,jj,kk,LI,k,j,i;
  PetscScalar dx,dy,dz,Lx,Ly,Lz;
  PetscScalar hx,hy,hz;

  Lx=2.0*rx;
  Ly=2.0*ry;
  Lz=2.0*rz;

  hx=Lx/nx;
  hy=Ly/ny;
  hz=Lz/nz;

  if(x<0)
    dx=x+Rx;
  else if(x>0)
    dx=x-Rx;
		  
  if(y<0)
    dy=y+Ry;
  else if(y>0)
    dy=y-Ry;

  if(z<0)
    dz=z+Rz;
  else if(z>0)
    dz=z-Rz;

  *X=fmod(dx,Lx);
  *Y=fmod(dy,Ly);
  *Z=fmod(dz,Lz);

   if(*X<0)
    *X=*X+Lx;
  if(*Y<0)
    *Y=*Y+Ly;
  if(*Z<0)
    *Z=*Z+Lz;

  *i0=floorf(*X/hx);
  *j0=floorf(*Y/hy);
  *k0=floorf(*Z/hz);
		  
  if(*i0>=nx)
    *i0=*i0-nx;
  if(*j0>=ny)
    *j0=*j0-ny;
  if(*k0>=nz)
    *k0=*k0-nz;

  *i1=ceilf(*X/hx);
  *j1=ceilf(*Y/hy);
  *k1=ceilf(*Z/hz);
		  
  if(*i1>=nx)
    *i1=*i1-nx;
  if(*j1>=ny)
    *j1=*j1-ny;
  if(*k1>=nz)
    *k1=*k1-nz;
  
}

///////////////////////////////////////////////////////////////////////////////////////
void EvaluateCorners(PetscScalar *arrVals,int nx,int ny,int nz,int i0,int j0,int k0,int i1,int j1,int k1,PetscScalar *c000,PetscScalar *c001,PetscScalar *c010,PetscScalar *c011,PetscScalar *c100,PetscScalar *c101,PetscScalar *c110,PetscScalar *c111)
{

  int LI;
  // c000
  LI= k0*nx*ny + j0*nx + i0;
  *c000=arrVals[LI];

   // c001
  LI= k1*nx*ny + j0*nx + i0;
  *c001=arrVals[LI];

   // c010
  LI= k0*nx*ny + j1*nx + i0;
  *c010=arrVals[LI];

 // c011
  LI= k1*nx*ny + j1*nx + i0;
  *c011=arrVals[LI];

 // c100
  LI= k0*nx*ny + j0*nx + i1;
  *c100=arrVals[LI];

 // c101
  LI= k1*nx*ny + j0*nx + i1;
  *c101=arrVals[LI];

 // c110
  LI= k0*nx*ny + j1*nx + i1;
  *c110=arrVals[LI];

 // c111
  LI= k1*nx*ny + j1*nx + i1;
  *c111=arrVals[LI];

}
///////////////////////////////////////////////////////////////////////////////////////
