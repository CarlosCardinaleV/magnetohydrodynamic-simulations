#include <float.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "agregatedFunctions.h"


// valores globales
#define GAMMA 2.0
#define NODOS 900


int mhd(void);

/**
 * Calcula una serie de parámetros de fluido a partir de
 * las entradas proporcionadas y las almacena en un solo arreglo.
 *
 * @param Bx Un puntero al inicio de un arreglo que contiene
 * los componentes x del campo magnético.
 * @param By Un puntero al inicio de un arreglo que contiene
 * los componentes y del campo magnético.
 * @param Bz Un puntero al inicio de un arreglo que contiene
 * los componentes z del campo magnético.
 * @param rho Un puntero al inicio de un arreglo que contiene
 * la densidad del fluido en cada nodo.
 * @param p Un puntero al inicio de un arreglo que contiene
 * la presión del fluido en cada nodo.
 * @param Vx Un puntero al inicio de un arreglo que contiene
 * los componentes x de la velocidad del fluido.
 * @param Vy Un puntero al inicio de un arreglo que contiene
 * los componentes y de la velocidad del fluido.
 * @param Vz Un puntero al inicio de un arreglo que contiene
 * los componentes z de la velocidad del fluido.
 *
 * @return Un puntero al inicio de un arreglo de tamaño 8*nNodos
 * que contiene los siguientes parámetros calculados para cada nodo:
 *         - Los componentes x, y, z del campo magnético.
 *         - La densidad del fluido.
 *         - Los componentes x, y, z del momento del fluido
 *           (rho*Vx, rho*Vy, rho*Vz).
 *         - La energía total del fluido, calculada como
 *           p/(gamma-1) + 0.5*rho*(V*V) + 0.5*(B*B).
 *         En caso de error en la asignación de memoria, la función
 *         imprimirá un mensaje de error y terminará el programa.
 */
double* getU(double* Bx, double* By, double* Bz,
double* rho, double* p, double* Vx, double* Vy, double* Vz);


/**
 * @brief Calcula la matriz F basada en la matriz U proporcionada.
 * 
 * Esta función calcula varios parámetros utilizando la matriz U dada, 
 * y guarda los resultados en una nueva matriz F. La función asume 
 * que U es un array 2D aplanado con un tamaño de 8*nNodos.
 * 
 * @param U Matriz de entrada, un array 2D aplanado. 
 *    U[0..nNodos-1] = Bx,
 *    U[nNodos..2*nNodos-1] = By,
 *    U[2*nNodos..3*nNodos-1] = Bz,
 *    U[3*nNodos..4*nNodos-1] = rho,
 *    U[4*nNodos..5*nNodos-1] = Vx,
 *    U[5*nNodos..6*nNodos-1] = Vy,
 *    U[6*nNodos..7*nNodos-1] = Vz,
 *    U[7*nNodos..8*nNodos-1] = E.
 * 
 * @return F Matriz de salida, un nuevo array 2D aplanado.
 * El llamante es responsable de liberar esta memoria.
 *    F[0..nNodos-1] siempre será 0,
 *    F[nNodos..2*nNodos-1] corresponde a By*Vx - Bx*Vy,
 *    F[2*nNodos..3*nNodos-1] corresponde a Bz*Vx - Bx*Vz,
 *    F[3*nNodos..4*nNodos-1] corresponde a rho*Vx,
 *    F[4*nNodos..5*nNodos-1] corresponde a rho*Vx^2 + p + B^2/2 - Bx^2,
 *    F[5*nNodos..6*nNodos-1] corresponde a rho*Vx*Vy - Bx*By,
 *    F[6*nNodos..7*nNodos-1] corresponde a rho*Vx*Vz - Bx*Bz,
 *    F[7*nNodos..8*nNodos-1] corresponde a
 *      (0.5*rho*V^2 + gamma*p/(gamma-1) + B^2)*Vx - Bx*(Bx*Vx +
 *      By*Vy + Bz*Vz).
 */
double* getF(double* U);

/**
 * @brief Genera una matriz D basada en la matriz de entrada U.
 *
 * @param U Puntero a una matriz dinámicamente asignada,
 * representando la matriz de entrada. 
 *        Las ocho filas de esta matriz corresponden a
 *        las variables Bx, By, Bz, rho, Vx, Vy, Vz y E,
 *        cada una de tamaño 'nNodos'.
 * 
 * @return Un puntero a una matriz dinámicamente asignada que
 *         representa la matriz D.
 *         Esta matriz consta de ocho filas, cada una de tamaño
 *         'nNodos'.
 *         La quinta fila corresponde a la velocidad en el eje
 *         X 'Vx' y la última fila es la presión 'p' dividida por
 *         la densidad 'rho'.
 *         El resto de las filas están llenas de ceros.
 *         Es responsabilidad del llamante liberar esta memoria
 *         una vez que ya no se necesite.
 *         Devuelve NULL si la asignación de memoria falla.
 */
double* getD(double* U);

/**
 * @brief Calcula y devuelve el arreglo LU.
 *
 * Esta función toma como entrada un puntero a un
 * arreglo de double (`U`) que representa un vector
 * aplanado. El vector original es una matriz 2D de
 * tamaño 8*nNodos. La función realiza operaciones
 * matemáticas específicas para calcular un nuevo vector
 * (`LU`), que luego se devuelve.
 * 
 * `U` se interpreta como una matriz aplanada y se accede
 * utilizando cálculos de desplazamiento para determinar
 * los índices equivalentes en la matriz 2D.
 * 
 * Las funciones `getF` y `getD` deben estar diseñadas
 * para trabajar con `U` en su forma aplanada. `F` y `D`
 * son vectores intermedios calculados dentro de la función.
 *
 * @param U Un puntero a un arreglo de double que representa
 *          la matriz `U` aplanada.
 *          Se supone que `U` tiene un tamaño de 8*nNodos.
 * @return Un puntero al arreglo `LU` recién calculado,
 *         que también tiene un tamaño de 8*nNodos.
 */
double* getLU(double* U);

/**
 * @brief Calcula los primeros tres pasos del método numérico
 * Runge-Kutta de cuarto orden.
 *
 * Esta función toma una matriz U y un intervalo de tiempo dt,
 * y calcula los primeros tres pasos del método Runge-Kutta
 * de cuarto orden, que se utilizan como entrada para el método
 * Adams-Bashforth-Moulton.
 *
 * @param U Puntero a la matriz de entrada.
 * @param dt El paso de tiempo en la simulación.
 *
 * @return Una estructura TresU que contiene los primeros tres
 * pasos (U1, U2, U3) de la simulación.
 *
 * @note Esta función libera la memoria asociada a las variables
 * temporales internamente, pero la responsabilidad de liberar
 * la memoria de las matrices resultantes U1, U2 y U3 recae sobre
 * el usuario de la función.
 */
TresU getFirst3U(double *U, double dt);



//  main function
int main(int argc, char* argv[]) {
    
    // para tomar el tiempo de ejecucion
    clock_t start_time = clock();

    // Inicializa MPI
    //MPI_Init(&argc, &argv);

    int nNodos = NODOS;
    double W = 0.0085;
    double* x = linspace(-1, 1, nNodos);
    double dt = 0.0005;
    double tMax = 0.20;

    double rho1 = 1.0;
    double Vx1 = 0.0;
    double Vy1 = 0.0;
    double Vz1 = 0.0;
    double p1 = 1.0;
    double Bx1 = 0.75;
    double By1 = 1.0;
    double Bz1 = 0.0;
    double rho2 = 0.125;
    double Vx2 = 0;
    double Vy2 = 0;
    double Vz2 = 0;
    double p2 = 0.1;
    double Bx2 = 0.75;
    double By2 = 0.75;
    double Bz2 = 0.75;
    double* rho = (double*) calloc(nNodos, sizeof(double));
    double* Vx = (double*) calloc(nNodos, sizeof(double));
    double* Vy = (double*) calloc(nNodos, sizeof(double));
    double* Vz = (double*) calloc(nNodos, sizeof(double));
    double* p = (double*) calloc(nNodos, sizeof(double));
    double* Bx = (double*) calloc(nNodos, sizeof(double));
    double* By = (double*) calloc(nNodos, sizeof(double));
    double* Bz = (double*) calloc(nNodos, sizeof(double));
    for (int i = 0; i < nNodos; ++i) {
        rho[i] = (rho2 + rho1) / 2 + ((rho2 - rho1) / 2) * tanh(x[i] / W);
        Vx[i] = (Vx2 + Vx1) / 2 + ((Vx2 - Vx1) / 2) * tanh(x[i] / W);
        Vy[i] = (Vy2 + Vy1) / 2 + ((Vy2 - Vy1) / 2) * tanh(x[i] / W);
        Vz[i] = (Vz2 + Vz1) / 2 + ((Vz2 - Vz1) / 2) * tanh(x[i] / W);
        p[i] = (p2 + p1) / 2 + ((p2 - p1) / 2) * tanh(x[i] / W);
        Bx[i] = (Bx2 + Bx1) / 2 + ((Bx2 - Bx1) / 2) * tanh(x[i] / W);
        By[i] = (By2 + By1) / 2 + ((By2 - By1) / 2) * tanh(x[i] / W);
        Bz[i] = (Bz2 + Bz1) / 2 + ((Bz2 - Bz1) / 2) * tanh(x[i] / W);
    }

    //double* U = (double*) malloc(8 * nNodos * sizeof(double));
    //U = getU(Bx, By, Bz, rho, p, Vx, Vy, Vz);
    double* tempU = getU(Bx, By, Bz, rho, p, Vx, Vy, Vz);
    double* U = calloc(8 * nNodos, sizeof(double));
    memcpy(U, tempU, 8 * nNodos * sizeof(double));
    free(tempU);
    
    TresU resultado = getFirst3U(U, dt);
    double* UNMenos3 = U;
    double* UNMenos2 = resultado.U1;
    double* UNMenos1 = resultado.U2;
    double* Un = resultado.U3;

    double* LUnMenos3 = getLU(UNMenos3);
    double* LUnMenos2 = getLU(UNMenos2);
    double* LUnMenos1 = getLU(UNMenos1);


    FILE *fp;
    fp = fopen("datos.txt", "w");

    double tiempo = 3 * dt;
    int m = 3;
    while (tMax > tiempo) {
        tiempo += dt;
        double* LUn = getLU(Un);
        double* UIM = sumAndScale(Un, dt, 55.0/24, LUn, -59.0/24, \
        LUnMenos1, 37.0/24, LUnMenos2, -9.0/24, LUnMenos3);
        double* LUIM = getLU(UIM);

        double* UNmasUno = sumAndScale(Un, dt, 9.0/24, LUIM, 19.0/24, \
        LUn, -5.0/24, LUnMenos1, 1.0/24, LUnMenos2);
        double error = 1.0;
        int n = 1;

        while (error > 20*DBL_EPSILON) {
            double* LUNmasUno = getLU(UNmasUno);
            double* UNmasUnoP = sumAndScale(Un, dt, 9.0/24, LUNmasUno, \
            19.0/24, LUn, -5.0/24, LUnMenos1, 1.0/24, LUnMenos2);
            error = maxAbsDiff(UNmasUno, UNmasUnoP);
            free(UNmasUno);
            UNmasUno = UNmasUnoP;
            free(LUNmasUno);
            n++;
        }

        free(UNMenos3);
        UNMenos3 = UNMenos2;
        UNMenos2 = UNMenos1;
        UNMenos1 = Un;
        Un = UNmasUno;

        free(LUnMenos3);
        LUnMenos3 = LUnMenos2;
        LUnMenos2 = LUnMenos1;
        LUnMenos1 = LUn;

        // Liberar las variables que no se liberaron en el ciclo
        free(UIM);
        free(LUIM);
    
        m++;
    }

    // Escribimos los datos en el archivo
    for(int i = 0; i < nNodos; i++){
        fprintf(fp, "%lf %lf\n", x[i], Un[i + 3*nNodos]);
    }
    // Separador para múltiples conjuntos de datos en gnuplot
    fprintf(fp, "\n\n");
    fclose(fp); // Cierra el archivo al final de la simulación

    // se libera memoria dinamica
    free(Un); free(UNMenos1); free(UNMenos2); free(UNMenos3);
    free(LUnMenos1); free(LUnMenos2); free(LUnMenos3);
    free(x); free(rho); free(Vx); free(Vy); free(Vz);
    free(p); free(Bx); free(By); free(Bz);

    // Finaliza MPI
    //MPI_Finalize();

    // obtine el tiempo final de ejecucion
    clock_t end_time = clock();
    // Calcula el tiempo de ejecucion en segundos
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // imprime el tiempo de ejecucion en segundos
    printf("Tiempo de ejecución: %.2f segundos\n", execution_time);

    return EXIT_SUCCESS;
}


double* getU(double* Bx, double* By, double* Bz, double* rho,
double* p, double* Vx, double* Vy, double* Vz) {
    int nNodos = NODOS;
    double gamma = GAMMA;

    // creamos los vectores necesarios para U
    double* Mx = (double*) calloc(nNodos, sizeof(double));
    double* My = (double*) calloc(nNodos, sizeof(double));
    double* Mz = (double*) calloc(nNodos, sizeof(double));
    double* V = (double*) calloc(nNodos, sizeof(double));
    double* B = (double*) calloc(nNodos, sizeof(double));
    double* E = (double*) calloc(nNodos, sizeof(double));

    // vector U es un array 1D que contiene todos
    // los vectores de entrada uno después del otro, es decir,
    // primero todos los elementos de Bx, luego todos los elementos
    // de By, y así sucesivamente (en realidad es matriz aplanada).
    double* U = (double*) calloc(8 * nNodos, sizeof(double));

    // revisamos si fueron creados correctamente todos los vectores
    if (Mx == NULL || My == NULL || Mz == NULL || V == NULL ||
    B == NULL || E == NULL || U == NULL) {
        fprintf(stderr, "Eror: failed to allocate memory!\n");
        exit(EXIT_FAILURE);
    }

    // se le asigna el valor par cada uno de los elementos de los
    // vectores
    for (int i = 0; i < nNodos; ++i) {
        Mx[i] = rho[i] * Vx[i];
        My[i] = rho[i] * Vy[i];
        Mz[i] = rho[i] * Vz[i];
        V[i] = sqrt(Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i]);
        B[i] = sqrt(Bx[i]*Bx[i] + By[i]*By[i] + Bz[i]*Bz[i]);
        E[i] = p[i]/(gamma-1) + 0.5*rho[i]*(V[i]*V[i]) + 0.5*(B[i]*B[i]);
    }

    // colocamos por medio de puntero el valor de cada vector para que
    // se guarden en U de manera seguida (recordando que U es una matriz
    // aplanada). 
    memcpy(U, Bx, nNodos*sizeof(double));
    memcpy(U + nNodos, By, nNodos*sizeof(double));
    memcpy(U + 2*nNodos, Bz, nNodos*sizeof(double));
    memcpy(U + 3*nNodos, rho, nNodos*sizeof(double));
    memcpy(U + 4*nNodos, Mx, nNodos*sizeof(double));
    memcpy(U + 5*nNodos, My, nNodos*sizeof(double));
    memcpy(U + 6*nNodos, Mz, nNodos*sizeof(double));
    memcpy(U + 7*nNodos, E, nNodos*sizeof(double));

    // liberamos la memoria de los vectores
    free(Mx); free(My); free(Mz); free(V); free(B); free(E);

    return U;
}


double* getF(double* U) {
    int nNodos = NODOS;
    double gamma = GAMMA;

    // creamos los vectores necesarios para F
    double* Bx = (double*) calloc(nNodos, sizeof(double));
    double* By = (double*) calloc(nNodos, sizeof(double));
    double* Bz = (double*) calloc(nNodos, sizeof(double));
    double* rho = (double*) calloc(nNodos, sizeof(double));
    double* Vx = (double*) calloc(nNodos, sizeof(double));
    double* Vy = (double*) calloc(nNodos, sizeof(double));
    double* Vz = (double*) calloc(nNodos, sizeof(double));
    double* E = (double*) calloc(nNodos, sizeof(double));

    // vector F es un array 1D que contiene todos
    // los vectores de entrada uno después del otro,
    // despues de hacerle los calculos respectivos
    // (en realidad es matriz aplanada).
    double* F = (double*) calloc(8 * nNodos, sizeof(double));

    if (Bx == NULL || By == NULL || Bz == NULL || rho == NULL || Vx == NULL \
    || Vy == NULL || Vz == NULL || E == NULL || F == NULL) {
        fprintf(stderr, "Error: failed to allocate memory!\n");
        exit(EXIT_FAILURE);
    }

    // Assigning input values
    memcpy(Bx, U, nNodos*sizeof(double));
    memcpy(By, U + nNodos, nNodos*sizeof(double));
    memcpy(Bz, U + 2*nNodos, nNodos*sizeof(double));
    memcpy(rho, U + 3*nNodos, nNodos*sizeof(double));
    memcpy(Vx, U + 4*nNodos, nNodos*sizeof(double));
    memcpy(Vy, U + 5*nNodos, nNodos*sizeof(double));
    memcpy(Vz, U + 6*nNodos, nNodos*sizeof(double));
    memcpy(E, U + 7*nNodos, nNodos*sizeof(double));

    // Calculating the values
    for (int i = 0; i < nNodos; ++i) {
        Vx[i] /= rho[i];
        Vy[i] /= rho[i];
        Vz[i] /= rho[i];
        
        double B = sqrt(Bx[i]*Bx[i] + By[i]*By[i] + Bz[i]*Bz[i]);
        double V = sqrt(Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i]);
        double p = (E[i] - 0.5*rho[i]*V*V - 0.5*B*B) * (gamma-1);
        
        F[i] = 0;
        F[i + nNodos] = By[i]*Vx[i] - Bx[i]*Vy[i];
        F[i + 2*nNodos] = Bz[i]*Vx[i] - Bx[i]*Vz[i];
        F[i + 3*nNodos] = rho[i]*Vx[i];
        F[i + 4*nNodos] = rho[i]*Vx[i]*Vx[i] + p + 0.5*B*B - Bx[i]*Bx[i];
        F[i + 5*nNodos] = rho[i]*Vx[i]*Vy[i] - Bx[i]*By[i];
        F[i + 6*nNodos] = rho[i]*Vx[i]*Vz[i] - Bx[i]*Bz[i];
        
        double F8a = (0.5*rho[i]*V*V + gamma*p/(gamma-1) + B*B)*Vx[i];
        double F8b = Bx[i] * (Bx[i]*Vx[i] + By[i]*Vy[i] + Bz[i]*Vz[i]);
        F[i + 7*nNodos] = F8a - F8b;
    }

    // limpiamos memoria de los valores intermediarios
    free(Bx);
    free(By);
    free(Bz);
    free(rho);
    free(Vx);
    free(Vy);
    free(Vz);
    free(E);

    return F;
}


double* getD(double* U) {
    int nNodos = NODOS;
    double gamma = GAMMA;

    // creamos los vectores necesarios para D
    double* Vx = (double*) calloc(nNodos, sizeof(double));
    double* Bx = (double*) calloc(nNodos, sizeof(double));
    double* By = (double*) calloc(nNodos, sizeof(double));
    double* Bz = (double*) calloc(nNodos, sizeof(double));
    double* rho = (double*) calloc(nNodos, sizeof(double));
    double* Vy = (double*) calloc(nNodos, sizeof(double));
    double* Vz = (double*) calloc(nNodos, sizeof(double));
    double* E = (double*) calloc(nNodos, sizeof(double));

    // vector D es un array 1D que contiene todos
    // los vectores de entrada uno después del otro,
    // despues de hacerle los calculos respectivos
    // (en realidad es matriz aplanada).
    double* D = (double*) calloc(8 * nNodos, sizeof(double));

    // revisamos que se hayan creado correctamente
    if (Vx == NULL || Bx == NULL || By == NULL || Bz == NULL || \
    rho == NULL || Vy == NULL || Vz == NULL || E == NULL || D == NULL) {
        fprintf(stderr, "Error: failed to allocate memory!\n");
        exit(EXIT_FAILURE);
    }

    memcpy(Bx, U, nNodos*sizeof(double));
    memcpy(By, U + nNodos, nNodos*sizeof(double));
    memcpy(Bz, U + 2*nNodos, nNodos*sizeof(double));
    memcpy(rho, U + 3*nNodos, nNodos*sizeof(double));
    memcpy(Vx, U + 4*nNodos, nNodos*sizeof(double));
    memcpy(Vy, U + 5*nNodos, nNodos*sizeof(double));
    memcpy(Vz, U + 6*nNodos, nNodos*sizeof(double));
    memcpy(E, U + 7*nNodos, nNodos*sizeof(double));

    double* B = (double*) calloc(nNodos, sizeof(double));
    double* V = (double*) calloc(nNodos, sizeof(double));
    double* p = (double*) calloc(nNodos, sizeof(double));
    
    // revisamos que se hayan creado correctamente
    if (B == NULL || V == NULL || p == NULL) {
        fprintf(stderr, "Error: failed to allocate memory!\n");
        exit(EXIT_FAILURE);
    }

    // calculamos y asignamos los valores a la matriz aplanada D
    for (int i = 0; i < nNodos; ++i) {
        Vx[i] /= rho[i];
        Vy[i] /= rho[i];
        Vz[i] /= rho[i];
        B[i] = sqrt(Bx[i]*Bx[i] + By[i]*By[i] + Bz[i]*Bz[i]);
        V[i] = sqrt(Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i]);
        p[i] = (E[i] - 0.5*rho[i]*V[i]*V[i] - 0.5*B[i]*B[i])*(gamma-1);
        D[4*nNodos + i] = Vx[i];
        D[7*nNodos + i] = p[i]/rho[i];
    }

    // para liberar la memoria dinamica
    free(Vx);
    free(Bx);
    free(By);
    free(Bz);
    free(rho);
    free(Vy);
    free(Vz);
    free(E);
    free(B);
    free(V);
    free(p);

    return D;
}


double* getLU(double* U) {
    int nNodos = NODOS;
    double gamma = GAMMA;

    double dx = 2.0 / (nNodos - 1);
    double etaV = 0.001;
    double etaT = 0.001;
    double a1 = 1.0 / 60, a2 = -3.0 / 20, a3 = 3.0 / 4, \
    a4 = -3.0 / 4, a5 = 3.0 / 20, a6 = -1.0 / 60;
    double b1 = 1.0 / 90, b2 = -3.0 / 20, b3 = 3.0 / 2, \
    b4 = -49.0 / 18, b5 = 3.0 / 2, b6 = -3.0 / 20, b7 = 1.0 / 90;

    double* F = getF(U);
    double* D = getD(U);
    
    double* dummy1 = (double*) calloc(8 * nNodos, sizeof(double));
    double* dummy2 = (double*) calloc(8 * nNodos, sizeof(double));
    
    // calcula dummy1
    #pragma omp parallel for
    for (int i = 3; i < nNodos - 3; ++i) {
        for (int j = 0; j < 8; ++j) {
            dummy1[j*nNodos + i] = a1*F[(j*nNodos) + i+3] + a2*F[(j*nNodos) + i+2] + \
            a3*F[(j*nNodos) + i+1] + a4*F[(j*nNodos) + i-1] + a5*F[(j*nNodos) + i-2] + \
            a6*F[(j*nNodos) + i-3];
        }
    }

    // calcula dummy2
    #pragma omp parallel for
    for (int i = 3; i < nNodos - 3; ++i) {
        for (int j = 0; j < 8; ++j) {
            dummy2[j*nNodos + i] = b1*D[(j*nNodos) + i+3] + b2*D[(j*nNodos) + i+2] + \
            b3*D[(j*nNodos) + i+1] + b4*D[(j*nNodos) + i] + b5*D[(j*nNodos) + i-1] + \
            b6*D[(j*nNodos) + i-2] + b7*D[(j*nNodos) + i-3];
        }
    }

    for (int i = 0; i < nNodos; ++i) {
        dummy2[4*nNodos + i] *= etaV * U[3*nNodos + i];
        dummy2[7*nNodos + i] *= etaT * U[3*nNodos + i];
    }

    // calcula LU
    double* LU = (double*) calloc(8 * nNodos, sizeof(double));
    #pragma omp parallel for
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < nNodos; ++j) {
            LU[i*nNodos + j] = dummy1[i*nNodos + j] / dx + dummy2[i*nNodos + j] / (dx * dx);
        }
    }

    free(F);
    free(D);
    free(dummy1);
    free(dummy2);

    return LU;
}


TresU getFirst3U(double *U, double dt) {
    double *LU, *LU1, *LU2, *LU3;
    double *U1fake, *U2fake, *U3fake;
    double *U1, *U2, *U3;
    TresU resultado;

    LU = getLU(U);
    U1fake = addAndScale(U, dt/2, LU);
    LU1 = getLU(U1fake);
    U2fake = addAndScale(U, dt/2, LU1);
    LU2 = getLU(U2fake);
    U3fake = addAndScale(U, dt, LU2);
    LU3 = getLU(U3fake);
    U1 = addAndScale4(U, dt, LU, 1/6, LU1, 1/3, LU2, 1/3, LU3, 1/6);

    free(LU); free(U1fake); free(U2fake); free(U3fake);
    free(LU1); free(LU2); free(LU3);

    LU = getLU(U1);
    U1fake = addAndScale(U1, dt/2, LU);
    LU1 = getLU(U1fake);
    U2fake = addAndScale(U1, dt/2, LU1);
    LU2 = getLU(U2fake);
    U3fake = addAndScale(U1, dt, LU2);
    LU3 = getLU(U3fake);
    U2 = addAndScale4(U1, dt, LU, 1/6, LU1, 1/3, LU2, 1/3, LU3, 1/6);

    free(LU); free(U1fake); free(U2fake); free(U3fake);
    free(LU1); free(LU2); free(LU3);

    LU = getLU(U2);
    U1fake = addAndScale(U2, dt/2, LU);
    LU1 = getLU(U1fake);
    U2fake = addAndScale(U2, dt/2, LU1);
    LU2 = getLU(U2fake);
    U3fake = addAndScale(U2, dt, LU2);
    LU3 = getLU(U3fake);
    U3 = addAndScale4(U2, dt, LU, 1/6, LU1, 1/3, LU2, 1/3, LU3, 1/6);

    free(LU); free(U1fake); free(U2fake); free(U3fake);
    free(LU1); free(LU2); free(LU3);

    resultado.U1 = U1;
    resultado.U2 = U2;
    resultado.U3 = U3;

    // Devuelve la estructura con los resultados
    return resultado;
}

