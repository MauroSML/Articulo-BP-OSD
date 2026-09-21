#ifndef OSD_H
#define OSD_H

#include <ap_fixed.h>
#include <ap_int.h>
#define AP_INT_MAX_W 1024

int const N = 429;
int const M = 882;

typedef ap_uint<1> bit_t;
typedef double value_t;
typedef ap_uint<N> col_vec_t;

// Declaraciones de funciones OSD
void decode(bit_t synd_in[N], value_t prob_ini_in[M], bit_t sol[M]);
int sort_columns(value_t prob[M], int sol[M]);
/*int make_H(bit_t H_T_cod[M][N], int sorted_cols[M], bit_t H_sol[N][N], int used_cols[N]);
void eliminacion_gaussiana(bit_t col1[N], bit_t col2[N], int elem_estud);
int busca_posicion_columna(bit_t columna_colocar[N], bit_t H_lu[N][N], int numDetec);
void eliminacion_salida(bit_t columna_eliminar[N], bit_t H_lu[N][N], int index_col, int numDetec);
int solve_ecuationSys(bit_t synd[N], bit_t H_s[N][N], int used_cols[N], bit_t sol[M]);*/

int make_H_vec(const col_vec_t H_T_cod[M], const int sorted_cols[M], col_vec_t H_sol[N], int used_cols[N]);
int busca_posicion_columna_vec(col_vec_t &columna_colocar, const col_vec_t H_lu[N], int numDetec);
void eliminacion_salida_vec(col_vec_t &columna_eliminar, col_vec_t H_lu[N], int index_col, int numDetec);
int solve_ecuationSys_vec(col_vec_t synd, const col_vec_t H_s[N], const int used_cols[N], bit_t sol[M]);

// Álgebra lineal
/*void multiplicarMatrizVector(const bit_t synd[N], const bit_t H_s[N][N], bit_t result[N]);
void transponer_matriz_rect(const bit_t entrada[N][M], bit_t salida[M][N]);
void transponer_matriz(bit_t entrada[N][N], bit_t salida[N][N]);*/

#endif // OSD_H