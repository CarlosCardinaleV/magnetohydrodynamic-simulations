#ifndef AGREGATEDFUNCTIONS_H
#define AGREGATEDFUNCTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/**
 * Struct TresU
 * 
 * Este struct se utiliza para almacenar tres vectores
 * de números de punto flotante. Los vectores son aplanados
 * (1D) y representan matrices 2D aplanadas. Los campos U1,
 * U2, y U3 son punteros a los vectores, y cada vector tiene 
 * un tamaño de 8 * NODOS.
 * 
 * Los vectores U1, U2, y U3 son generados por la función
 * getFirst3U y representan los primeros tres términos de
 * una secuencia utilizada en un método numérico para resolver
 * ecuaciones diferenciales.
 */
typedef struct {
    double* U1;
    double* U2;
    double* U3;
} TresU;


/**
 * @brief Genera una secuencia de números uniformemente espaciados
 *        entre un valor inicial y un valor final.
 * @param start El valor inicial de la secuencia.
 * @param end El valor final de la secuencia.
 * @param num El número total de valores en la secuencia.
 * @return un puntero a un array dinámicamente asignado 
 *         que contiene la secuencia generada.
 *         La longitud del array es igual al valor de entrada 'num'.
 *         Es responsabilidad del llamante liberar esta memoria
 *         una vez que ya no sea necesaria.
 *         Devuelve NULL si la asignación de memoria falla.
 */
double* linspace(double start, double end, int num);

/**
 * @brief Libera la memoria asignada a una matriz 2D de doubles.
 *
 * @param array2D El puntero a la matriz 2D a liberar.
 * @param rows La cantidad de filas de la matriz.
 * 
 * Esta función libera la memoria previamente asignada a una matriz 2D.
 * Cada fila de la matriz se libera individualmente, y luego el puntero
 * principal a la matriz.
 */
void free2DArray(double** array2D, int rows);

/**
 * Función addAndScale
 * 
 * Esta función toma un vector (array unidimensional) U, un 
 * factor y un vector LU. 
 * Devuelve un nuevo vector que es el resultado de sumar el
 * vector original U y el vector LU 
 * escalado por el factor proporcionado. Esto es equivalente
 * a U + factor * LU en algebra matricial.
 * 
 * @param U Un puntero a la matriz aplanada original
 * @param factor Un número para escalar el vector LU
 * @param LU Un puntero a la matriz aplanada LU
 * @return Un puntero al vector resultante
 */
// nota: este funcion fue recomendada por chatgpt para compactar
// mejor el codigo, pero fue modificada para que funcione en mi caso
double* addAndScale(double* U, double factor, double* LU);

/**
 * Función addAndScale4
 * 
 * Esta función toma un vector (array unidimensional) U, un factor, cuatro
 * vectores LU, LU1, LU2, LU3, y cuatro factores de escala correspondientes.
 * Devuelve un nuevo vector que es el resultado de sumar el vector original
 * U y los cuatro vectores LU, LU1, LU2, LU3 escalados por sus factores respectivos. 
 * Esto es equivalente a U + factor * (LU_factor*LU + LU1_factor*LU1 +
 * LU2_factor*LU2 + LU3_factor*LU3) en algebra matricial.
 * 
 * @param U Un puntero a la matriz aplanada original
 * @param factor Un número para escalar el resultado de la suma de vectores escalados
 * @param LU Un puntero a la matriz aplanada LU
 * @param LU_factor Un número para escalar el vector LU
 * @param LU1 Un puntero a la matriz aplanada LU1
 * @param LU1_factor Un número para escalar el vector LU1
 * @param LU2 Un puntero a la matriz aplanada LU2
 * @param LU2_factor Un número para escalar el vector LU2
 * @param LU3 Un puntero a la matriz aplanada LU3
 * @param LU3_factor Un número para escalar el vector LU3
 * @return Un puntero al vector resultante
 */
// nota: este funcion fue recomendada por chatgpt para compactar
// mejor el codigo, pero fue modificada para que funcione en mi caso
double* addAndScale4(double* U, double factor, double* LU, double LU_factor,
double* LU1, double LU1_factor, double* LU2, double LU2_factor, double* LU3,
double LU3_factor);


/**
 * printMatrix
 *
 * Imprime una matriz 2D de doubles almacenada en memoria como un solo bloque contiguo.
 * La matriz se asume en formato de filas mayores (row-major), que es el estándar en C.
 * En este formato, los elementos en una fila están contiguos en memoria.
 *
 * @param matrix Un puntero al primer elemento de la matriz. 
 *               Esto debe apuntar a un bloque de memoria que contenga al
 *               menos nrows*ncols elementos.
 * @param nrows El número de filas en la matriz.
 * @param ncols El número de columnas en la matriz.
 *
 * La función no devuelve ningún valor. En su lugar, imprime los elementos
 * de la matriz en la salida estándar (por lo general, la terminal),
 * con cada fila en una línea separada y los elementos de cada fila
 * separados por espacios.
 */
void printMatrix(double* matrix, int nrows, int ncols);

/**
 * @brief Imprime los elementos de un vector.
 *
 * Esta función recorre cada elemento del vector e imprime su valor. 
 * Cada valor se imprime con un espacio después y todos los valores
 * están encerrados entre corchetes. El vector se imprime en una única línea.
 *
 * @param vector Un puntero al primer elemento del vector de doubles que se va a imprimir.
 * @param size El número de elementos en el vector. 
 */
void printVector(double* vector, int size);


/**
 * Función printTresU
 * 
 * Esta función imprime los valores contenidos en el struct TresU.
 * 
 * @param res Struct TresU que contiene los punteros a los vectores U1, U2, U3.
 */
void printTresU(TresU res);

/**
 * Realiza la suma y escala de los arreglos U, LU, LUnMenos1, LUnMenos2 y LUnMenos3.
 *
 * @param U            Arreglo de entrada U.
 * @param dt           Valor del paso de tiempo.
 * @param scale1       Escala para U.
 * @param LU           Arreglo de entrada LU.
 * @param scale2       Escala para LU.
 * @param LUnMenos1    Arreglo de entrada LUnMenos1.
 * @param scale3       Escala para LUnMenos1.
 * @param LUnMenos2    Arreglo de entrada LUnMenos2.
 * @param scale4       Escala para LUnMenos2.
 * @param LUnMenos3    Arreglo de entrada LUnMenos3.
 * @return             Arreglo resultado de la suma y escala.
 */
double* sumAndScale(double* U, double dt, double scale1, \
double* LU, double scale2, double* LUnMenos1, double scale3, \
double* LUnMenos2, double scale4, double* LUnMenos3);

/**
 * Calcula la diferencia absoluta máxima entre los arreglos A y B.
 *
 * @param A    Arreglo de entrada A.
 * @param B    Arreglo de entrada B.
 * @return     Diferencia absoluta máxima entre A y B.
 */
double maxAbsDiff(double* A, double* B);

#endif
