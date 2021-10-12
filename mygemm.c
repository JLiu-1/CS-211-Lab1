#include "mygemm.h"

/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/

//Register Reuse part 1
void dgemm0(const double* A, const double* B, double* C, const int n)
{
    int i,j,k;
    for (i=0; i<n; i++) 
        for (j=0; j<n; j++) 
            for (k=0; k<n; k++) 
               C[i*n+j] += A[i*n+k] * B[k*n+j];
}

void dgemm1(const double *A, const double *B, double *C, const int n) 
{
    int i,j,k;
    for (i=0; i<n; i++) 
        for (j=0; j<n; j++) { 
            register double r = C[i*n+j] ; 
        for (k=0; k<n; k++) 
            r += A[i*n+k] * B[k*n+j]; 
        C[i*n+j] = r; 
    }
}
//Register Reuse part 1 End

//Register Reuse part 2
void dgemm2(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i = 0; i < n; i+=2) {
        for (j = 0; j < n; j+=2) {
            for (k = 0; k < n; k+=2) {
                C[i*n + j] = A[i*n + k]*B[k*n + j] + A[i*n + k+1]*B[(k+1)*n + j] + C[i*n + j];
                C[(i+1)*n + j] = A[(i+1)*n + k]*B[k*n + j] + A[(i+1)*n + k+1]*B[(k+1)*n + j] + C[(i+1)*n + j];
                C[i*n + (j+1)] = A[i*n + k]*B[k*n + (j+1)] + A[i*n + k+1]*B[(k+1)*n + (j+1)] + C[i*n + (j+1)];
                C[(i+1)*n + (j+1)] = A[(i+1)*n + k]*B[k*n + (j+1)] + A[(i+1)*n + k+1]*B[(k+1)*n + (j+1)] + C[(i+1)*n + (j+1)];
            }
        } 
    } 
}
//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i = 0; i < n; i += 3) {
        for (j = 0; j < n; j += 3) {
            register int t = i*n+j; // 00
            register int tt = t+n; // 01
            register int ttt = tt+n; // 02
            register double rc00 = C[t];
            register double rc01 = C[t+1];
            register double rc02 = C[t+2];
            register double rc10 = C[tt];
            register double rc11 = C[tt+1];
            register double rc12 = C[tt+2];
            register double rc20 = C[ttt];
            register double rc21 = C[ttt+1];
            register double rc22 = C[ttt+2];

            for (k = 0; k < n; k += 3) {
                register int ta = i*n+k;
                register int tta = ta+n;
                register int ttta = tta+n;
                register int tb = k*n+j;
                register int ttb = tb+n;
                register int tttb = ttb+n;

                register double R1 = A[ta]; // ra00
                register double R2 = A[tta]; // ra10
                register double R3 = A[ttta]; // ra20
                register double R4 = B[tb]; // rb00
                register double R5 = B[tb+1]; // rb01
                register double R6 = B[tb+2]; // rb02
                
                rc00 += R1 * R4;
                rc01 += R1 * R5;
                rc02 += R1 * R6;
                rc10 += R2 * R4;
                rc11 += R2 * R5;
                rc12 += R2 * R6;
                rc20 += R3 * R4;
                rc21 += R3 * R5;
                rc22 += R3 * R6;

                R1 = A[ta+1];
                R2 = A[tta+1];
                R3 = A[ttta+1];
                R4 = B[ttb];
                R5 = B[ttb+1];
                R6 = B[ttb+2];
                rc00 += R1 * R4;
                rc01 += R1 * R5;
                rc02 += R1 * R6;
                rc10 += R2 * R4;
                rc11 += R2 * R5;
                rc12 += R2 * R6;
                rc20 += R3 * R4;
                rc21 += R3 * R5;
                rc22 += R3 * R6;

                R1 = A[ta+2];
                R2 = A[tta+2];
                R3 = A[ttta+2];
                R4 = B[tttb];
                R5 = B[tttb+1];
                R6 = B[tttb+2];
                rc00 += R1 * R4;
                rc01 += R1 * R5;
                rc02 += R1 * R6;
                rc10 += R2 * R4;
                rc11 += R2 * R5;
                rc12 += R2 * R6;
                rc20 += R3 * R4;
                rc21 += R3 * R5;
                rc22 += R3 * R6;
            }
            C[t] = rc00;
            C[t+1] = rc01;
            C[t+2] = rc02;
            C[tt] = rc10;
            C[tt+1] = rc11;
            C[tt+2] = rc12;
            C[ttt] = rc20;
            C[ttt+1] = rc21;
            C[ttt+2] = rc22;
        }
    }

}
//Register Reuse part 3 End

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n) 
{
    /* ijk – simple triple loop algorithm with simple single register reuse*/ 
    int i, j, k;
    for (i=0; i<n; i++) 
        for (j=0; j<n; j++) { 
            register double sum = 0.0; 
            for (k=0; k<n; k++) 
                sum += A[i*n+k] * B[i*k+j]; 
            C[i*n+j] = sum;
        }
}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i, j, k, i1, j1, k1;
    for (i = 0; i < n; i+=b)
        for (j = 0; j < n; j+=b)
            for (k = 0; k < n; k+=b)
                /* B x B mini matrix multiplications */
                for (i1 = i; i1 < i+b; i++)
                    for (j1 = j; j1 < j+b; j++) {
                        register double r = C[i1*n+j1];
                        for (k1 = k; k1 < k+b; k++)
                            r += A[i1*n + k1] * B[k1*n + j1];
                        C[i1*n+j1] = r;
                    }
}

void jik(const double *A, const double *B, double *C, const int n) 
{
    /* jik */
    int i, j, k;
    for (j=0; j<n; j++) {
        for (i=0; i<n; i++) {
            register double sum = 0.0; 
            for (k=0; k<n; k++)
                sum += A[i*n+k] * B[i*k+j];
            C[i*n+j] = sum;
        }
    }
}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i, j, k, i1, j1, k1;
    for (j = 0; j < n; j+=b)
        for (i = 0; i < n; i+=b)
            for (k = 0; k < n; k+=b)
                /* B x B mini matrix multiplications */
                for (j1 = j; j1 < j+b; j++)
                    for (i1 = i; i1 < i+b; i++) {
                        register double r = C[i1*n+j1];
                        for (k1 = k; k1 < k+b; k++)
                            r += A[i1*n + k1] * B[k1*n + j1];
                        C[i1*n+j1] = r;
                    }
}

void kij(const double *A, const double *B, double *C, const int n) 
{
    /* kij */
    int i, j, k;
    for (k=0; k<n; k++) {
        for (i=0; i<n; i++) {
            register double r = A[i*n+k];
            for (j=0; j<n; j++)
                C[i*n+j] += r * B[i*k+j];   
        }
    }
}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i, j, k, i1, j1, k1;
    for (k = 0; k < n; k+=b)
        for (i = 0; i < n; i+=b)
            for (j = 0; j < n; j+=b)
                /* B x B mini matrix multiplications */
                for (k1 = k; k1 < k+b; k++)
                    for (i1 = i; i1 < i+b; i++) {
                        register double r = A[i1*n + k1];
                        for (j1 = j; j1 < j+b; j++)
                            C[i1*n+j1] += r * B[k1*n + j1];
                    }
}


void ikj(const double *A, const double *B, double *C, const int n) 
{
    /* ikj */
    int i, j, k;
    for (i=0; i<n; i++) {
        for (k=0; k<n; k++) {
            register double r = A[i*n+k];
            for (j=0; j<n; j++)
                C[i*n+j] += r * B[i*k+j];
        }
    }
}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i, j, k, i1, j1, k1;
    for (i = 0; i < n; i+=b)
        for (k = 0; k < n; k+=b)
            for (j = 0; j < n; j+=b)
                /* B x B mini matrix multiplications */
                for (i1 = i; i1 < i+b; i++)
                    for (k1 = k; k1 < k+b; k++) {
                        register double r = A[i1*n + k1];
                        for (j1 = j; j1 < j+b; j++)
                            C[i1*n+j1] += r * B[k1*n + j1];
                    }
}

void jki(const double *A, const double *B, double *C, const int n) 
{
    /* jki */
    int i, j, k;
    for (j=0; j<n; j++) {
        for (k=0; k<n; k++) {
            register double r = B[i*k+j];
            for (i=0; i<n; i++)
                C[i*n+j] += A[i*n+k] * r;
        }
    }  
}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i, j, k, i1, j1, k1;
    for (j = 0; j < n; j+=b)
        for (k = 0; k < n; k+=b)
            for (i = 0; i < n; i+=b)
                /* B x B mini matrix multiplications */
                for (j1 = j; j1 < j+b; j++)
                    for (k1 = k; k1 < k+b; k++) {
                        register double r = B[k1*n + j1];
                        for (i1 = i; i1 < i+b; i++)
                            C[i1*n+j1] += A[i1*n + k1] * r;
                    }
}

void kji(const double *A, const double *B, double *C, const int n) 
{
    /* kji */
    int i, j, k;
    for (k=0; k<n; k++) {
        for (j=0; j<n; j++) {
            register double r = B[i*k+j];
            for (i=0; i<n; i++)
                C[i*n+j] += A[i*n+k] * r;
        }
    }
}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i, j, k, i1, j1, k1;
    for (k = 0; k < n; k+=b)
        for (j = 0; j < n; j+=b)
            for (i = 0; i < n; i+=b)
                /* B x B mini matrix multiplications */
                for (k1 = k; k1 < k+b; k++)
                    for (j1 = j; j1 < j+b; j++) {
                        register double r = B[k1*n + j1];
                        for (i1 = i; i1 < i+b; i++)
                            C[i1*n+j1] += A[i1*n + k1] * r;
                    }
}
//Cache Reuse part 3 End 

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{
    int i, j, k, i1, j1, k1;
    for (i = 0; i < n; i+=b)
        for (j = 0; j < n; j+=b)
            for (k = 0; k < n; k+=b)
                /* B x B mini matrix multiplications */
                for (i1 = i; i1 < i+b; i += 2)
                    for (j1 = j; j1 < j+b; j += 2) {
                        register int t = i1*n+j1; register int tt = t+n; 
                        register double c00 = C[t]; register double c01 = C[t+1];  register double c10 = C[tt]; register double c11 = C[tt+1];
                        for (k1 = k; k1 < k+b; k += 2) {
                            /* 2 by 2 mini matrix multiplication using registers*/
                            register int ta = i1*n+k1; register int tta = ta+n; register int tb = k1*n+j1; register int ttb = tb+n;
                            register double a00 = A[ta]; register double a10 = A[tta]; register double b00 = B[tb]; register double b01 = B[tb+1]; 

                            c00 += a00*b00 ; c01 += a00*b01 ; c10 += a10*b00 ; c11 += a10*b01 ;
                            a00 = A[ta+1]; a10 = A[tta+1]; b00 = B[ttb]; b01 = B[ttb+1];
                            c00 += a00*b00 ; c01 += a00*b01 ; c10 += a10*b00 ; c11 += a10*b01 ;
                        }
                        C[t] = c00;
                        C[t+1] = c01;
                        C[tt] = c10;
                        C[tt+1] = c11;
                    }
}
