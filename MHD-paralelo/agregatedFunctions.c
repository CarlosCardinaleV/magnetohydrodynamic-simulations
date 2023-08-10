#include "agregatedFunctions.h"
// valores globales
#define GAMMA 2.0
#define NODOS 900

double* linspace(double start, double end, int num) {
    double* linspaced = NULL;
    double delta = 0.0;

    // If there is only one element, then return the start
    if (num == 1) {
        linspaced = (double*) calloc(1, sizeof(double));
        linspaced[0] = start;
        return linspaced;
    }

    delta = (end - start) / (num - 1);
    linspaced = (double*) calloc(num, sizeof(double));
    if(linspaced == NULL){
        fprintf(stderr, "Error: memory not allocated!\n");
        exit(EXIT_FAILURE);
    }
    for(int i=0; i < num; i++){
        linspaced[i] = start + delta * i;
    }
    return linspaced;
}

double* addAndScale(double* U, double factor, double* LU) {
    int size = 8 * NODOS;
    double* result = (double*) calloc(size, sizeof(double));

    #pragma omp parallel for
    for (int i = 0; i < size; i++) {
        result[i] = U[i] + factor * LU[i];
    }

    return result;
}

double* addAndScale4(double* U, double factor, double* LU, \
double LU_factor, double* LU1, double LU1_factor, double* LU2, \
double LU2_factor, double* LU3, double LU3_factor) {
    int size = 8 * NODOS;
    double* result = (double*) calloc(size, sizeof(double));

    #pragma omp parallel for
    for (int i = 0; i < size; i++) {
        result[i] = U[i] + factor * (LU_factor * LU[i] + LU1_factor * \
        LU1[i] + LU2_factor * LU2[i] + LU3_factor * LU3[i]);
    }

    return result;
}


void printMatrix(double* matrix, int nrows, int ncols) {
    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            printf("%f ", matrix[i*ncols + j]);
        }
        printf("\n");
    }
}

void printVector(double* vector, int size) {
    printf("[ ");
    for (int i = 0; i < size; i++) {
        printf("%f ", vector[i]);
    }
    printf("]\n");
}

void printTresU(TresU res) {
    int size = 8 * NODOS;

    printf("U1:\n");
    for (int i = 0; i < size; i++) {
        printf("%f ", res.U1[i]);
    }
    printf("\n");

    printf("U2:\n");
    for (int i = 0; i < size; i++) {
        printf("%f ", res.U2[i]);
    }
    printf("\n");

    printf("U3:\n");
    for (int i = 0; i < size; i++) {
        printf("%f ", res.U3[i]);
    }
    printf("\n");
}

double* sumAndScale(double* U, double dt, double scale1, \
double* LU, double scale2, double* LUnMenos1, double scale3, \
double* LUnMenos2, double scale4, double* LUnMenos3) {
    int nNodos = NODOS;
    double* result = (double*) malloc(8 * nNodos * sizeof(double));
    for (int i = 0; i < 8 * nNodos; i++) {
        result[i] = U[i] + dt * (scale1 * LU[i] + scale2 * \
        LUnMenos1[i] + scale3 * LUnMenos2[i] + scale4 * LUnMenos3[i]);
    }
    return result;
}

double maxAbsDiff(double* A, double* B) {
    int nNodos = NODOS;
    double maxDiff = 0;
    for (int i = 0; i < 8 * nNodos; i++) {
        double diff = fabs(A[i] - B[i]);
        if (diff > maxDiff) {
            maxDiff = diff;
        }
    }
    return maxDiff;
}
