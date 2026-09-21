#include "OSD.h"

/*void read_synd(bit_t synd_in[N], bit_t synd[N]){
    for(int i = 0; i < N; i++){
        synd[i] = synd_in[i];
    }
}*/

inline col_vec_t read_synd_vec(const bit_t synd_in[N]){
    #pragma HLS INLINE

    col_vec_t s = 0;

    READ_SYND: for(int i = 0; i < N;i++){
        #pragma HLS PIPELINE II=1
        s[i] = synd_in[i];
    }

    return s;
}

void read_prob_ini(value_t prob_ini_in[M], value_t prob_ini[M]){
    for(int i = 0; i < M; i++){
        #pragma HLS PIPELINE II=1
        prob_ini[i] = prob_ini_in[i];
    }
}

void decode(bit_t synd_in[N], value_t prob_ini_in[M], bit_t sol[M]) {

#pragma HLS INTERFACE mode=m_axi bundle=gmem0 port=synd_in     depth=N
#pragma HLS INTERFACE mode=m_axi bundle=gmem1 port=prob_ini_in depth=M
#pragma HLS INTERFACE mode=m_axi bundle=gmem2 port=sol         depth=M

    value_t prob_ini[M];
    col_vec_t H_sol[N];
    static const bit_t H[N][M] = {
        #include "H1.h"
    };
    col_vec_t H_T_cod[M];

    int sorted_cols[M];
    int used_cols[N];

#pragma HLS BIND_STORAGE variable=H_T_cod type=RAM_1P impl=BRAM
// Para forzar a que se guarden en BRAM y no gaste todas las LUTS

#pragma HLS ARRAY_PARTITION variable=H           complete dim=1
#pragma HLS ARRAY_PARTITION variable=H_sol       complete dim=1

#pragma HLS ARRAY_PARTITION variable=prob_ini    cyclic factor=16
#pragma HLS ARRAY_PARTITION variable=sorted_cols cyclic factor=16
#pragma HLS ARRAY_PARTITION variable=used_cols   complete dim=1

    TRANSP_H: for (int j = 0; j < M; j++) {
        #pragma HLS PIPELINE II=1
        col_vec_t col_words = 0;
        for (int i = 0; i < N; i++) {
            #pragma HLS UNROLL
            col_words[i] = H[i][j];
        }
        H_T_cod[j] = col_words;
    }

    read_prob_ini(prob_ini_in, prob_ini);
    col_vec_t synd = read_synd_vec(synd_in);

    sort_columns(prob_ini, sorted_cols);
    make_H_vec(H_T_cod, sorted_cols, H_sol, used_cols);
    solve_ecuationSys_vec(synd, H_sol, used_cols, sol);
}

int sort_columns(value_t prob[M], int sol[M]) {

    // Indices
    INIT_INDICES: for (int i = 0; i < M; ++i) {
        #pragma HLS UNROLL
        sol[i] = i;
    }

    int temp_indices[M];
    #pragma HLS ARRAY_PARTITION variable=temp_indices cyclic factor=16

    SORT_OUTER_LOOP: for (int width = 1; width < M; width *= 2) {
        #pragma HLS LOOP_FLATTEN off
        SORT_INNER_LOOP: for (int i = 0; i < M; i += 2 * width) {
            int left = i;
            int right = (i + width < M) ? (i + width) : M;
            int end = (i + 2 * width < M) ? (i + 2 * width) : M;

            int l = left;
            int r = right;
            int t = left;

            SORT_WHILE_1: while (l < right && r < end) {
                #pragma HLS PIPELINE II=1
                if (prob[sol[l]] <= prob[sol[r]]) {
                    temp_indices[t++] = sol[l++];
                }
                else {
                    temp_indices[t++] = sol[r++];
                }
            }

            SORT_WHILE_2: while (l < right) {
                #pragma HLS PIPELINE II=1
                temp_indices[t++] = sol[l++];
            }

            SORT_WHILE_3: while (r < end) {
                #pragma HLS PIPELINE II=1
                temp_indices[t++] = sol[r++];
            }

            UPDATE_SOL: for (int k = left; k < end; k++) {
                #pragma HLS PIPELINE II=1
                sol[k] = temp_indices[k];
            }
        }
    }

    return 0;
}


/*int make_H(bit_t H_T_cod[M][N], int sorted_cols[M], bit_t H_sol[N][N], int used_cols[N]) {

    bit_t H_T_lu[N][N];

    INIT_H_OUTER: for (int i = 0; i < N; i++) {
        used_cols[i] = -1;
        INIT_H_INNER: for (int j = 0; j < N; j++) {
            H_T_lu[i][j] = 0;
            H_sol[i][j] = 0;
        }
    }


    int asigned_cols = 0;
    int i = 0;
    FIND_COLS_WHILE: while (i < M && asigned_cols < N) {

        // Para no modificar H_T_cod
        int col_temp[N];
        for(int k = 0; k < N;k++){
            col_temp[k] = H_T_cod[sorted_cols[i]][k];
        }

        int pos_insert = busca_posicion_columna(H_T_cod[sorted_cols[i]], H_T_lu, N);

        if (pos_insert < N) {
            // En caso de que sea una posicion valida, se coloca la columna
            used_cols[pos_insert] = sorted_cols[i];
            asigned_cols += 1;
            COPY_COL_LOOP: for (int j = 0; j < N; j++) {
                H_T_lu[pos_insert][j] = H_T_cod[sorted_cols[i]][j];
            }
        }

        // Devolvemos al estado original H_T_cod
        for(int k = 0; k < N;k++){
            H_T_cod[sorted_cols[i]][k] = col_temp[k];
        }
        i++;
    }


    ELIM_OUTER_LOOP: for (int col = 0; col < N; col++) {
        eliminacion_salida(H_T_lu[col], H_T_lu, col, N);

        COPY_H_SOL_LOOP: for (int j = 0; j < N; j++) {
            H_sol[j][col] = H_T_lu[col][j];
        }
    }

    return 0;
}*/

int make_H_vec(const col_vec_t H_T_cod[M], const int sorted_cols[M], col_vec_t H_sol[N], int used_cols[N]) {

    col_vec_t H_T_lu[N];
    #pragma HLS BIND_STORAGE variable=H_T_lu type=RAM_1P impl=BRAM

    INIT_H: for(int i = 0; i < N;i++){
        #pragma HLS UNROLL
        used_cols[i] = -1;
        H_T_lu[i] = 0;
        H_sol[i] = 0;
    }

    int asigned_cols = 0;
    int i = 0;

    FIND_COLS: while(i < M && asigned_cols < N) {

        col_vec_t col_curr = H_T_cod[sorted_cols[i]];
        int pos_insert = busca_posicion_columna_vec(col_curr, H_T_lu, N);

        if (pos_insert < N){
            used_cols[pos_insert] = sorted_cols[i];
            asigned_cols++;
            H_T_lu[pos_insert] = col_curr;
        }
        i++;
    }


    ELIM_OUTER_LOOP: for (int col = 0; col < N; col++) {

        col_vec_t col_elim = H_T_lu[col];
        eliminacion_salida_vec(col_elim, H_T_lu, col, N);

        COPY_H_SOL_LOOP: for (int j = 0; j < N; j++) {
            #pragma HLS UNROLL
            H_sol[j][col] = H_T_lu[col][j];
        }
    }

    return 0;
}

/*void eliminacion_gaussiana(bit_t col1[N], bit_t col2[N], int elem_estud) {

    #pragma HLS INLINE

    GAUSS_ELIM_LOOP: for (int i = 0; i < N; i++) {
        #pragma HLS UNROLL 
        col2[i] = col1[i] ^ col2[i];
    }

    col2[elem_estud] = 1;
}*/

inline void eliminacion_gaussiana_vec(col_vec_t col1, col_vec_t &col2, int elem_estud){
    #pragma HLS INLINE
    col2 ^= col1;
    col2[elem_estud] = 1;
}

/*int busca_posicion_columna(bit_t columna_colocar[N], bit_t H_lu[N][N], int numDetec) {
    int j = 0;

    SEARCH_POS_WHILE: while (j < numDetec && (columna_colocar[j] == 0 || H_lu[j][j] == 1)) {
        if (columna_colocar[j] == 1 && H_lu[j][j] == 1) {
            eliminacion_gaussiana(H_lu[j], columna_colocar, j);
        }
        j += 1;
    }
    return j;
}*/

int busca_posicion_columna_vec(col_vec_t &columna_colocar, const col_vec_t H_lu[N], int numDetec){

    int j = 0;
    SEARCH_POS_WHILE: while(j < numDetec && (columna_colocar[j] == 0 || H_lu[j][j] == 1)){
        #pragma HLS PIPELINE II=1
        if(columna_colocar[j] == 1 && H_lu[j][j] == 1){
            eliminacion_gaussiana_vec(H_lu[j], columna_colocar, j);
        }
        j++;
    }
    return j;
}

/*void eliminacion_salida(bit_t columna_eliminar[N], bit_t H_lu[N][N], int index_col, int numDetec) {

    OUTPUT_ELIM_LOOP: for (int j = index_col + 1; j < numDetec; j++) {

        if (columna_eliminar[j] == 1 && H_lu[j][j] == 1) {
            eliminacion_gaussiana(H_lu[j], columna_eliminar, j);
        }
    }

    for(int i = 0;i <numDetec; i++){
        H_lu[index_col][i] = columna_eliminar[i];
    }
}*/

void eliminacion_salida_vec(col_vec_t &columna_eliminar, col_vec_t H_lu[N], int index_col, int numDetec){
    OUTPUT_ELIM: for(int j = 0; j < N; j++){
        #pragma HLS PIPELINE II=1
        if (j >= index_col + 1 && j < numDetec){
            if (columna_eliminar[j] == 1 && H_lu[j][j] == 1){
                eliminacion_gaussiana_vec(H_lu[j], columna_eliminar, j);
            }
        }
    }
    H_lu[index_col] = columna_eliminar;
}

/*int solve_ecuationSys(bit_t synd[N], bit_t H_s[N][N], int used_cols[N], bit_t sol[M]) {
    bit_t result[N];

    multiplicarMatrizVector(synd, H_s, result);

    INIT_SOL_LOOP: for (int i = 0; i < M; i++) {
        sol[i] = 0;
    }

    MAP_SOL_LOOP: for (int i = 0; i < N; i++) {
        if (used_cols[i] >= 0 && used_cols[i] < M) {
            sol[used_cols[i]] = result[i];
        }
    }

    return 0;
}*/

inline void  multiplicarMatrizVector_vec(const col_vec_t synd, const col_vec_t H_s[N], bit_t result[N]){
    #pragma HLS INLINE

    MULT_MAT_VEC: for(int i = 0; i < N;i++){
        #pragma HLS PIPELINE II=1
        col_vec_t combined = H_s[i] & synd;

        bit_t parity = 0;
        PARITY_LOOP: for(int k = 0; k<N;k++){
            #pragma HLS UNROLL
            parity^=combined[k];
        }

        result[i] = parity;
    }
}

int solve_ecuationSys_vec(const col_vec_t synd, const col_vec_t H_s[N], const int used_cols[N], bit_t sol[M]) {
    bit_t result[N];

    #pragma HLS ARRAY_PARTITION variable=result complete dim=1

    multiplicarMatrizVector_vec(synd, H_s, result);

    INIT_SOL: for (int i = 0; i < M; i++) {
        #pragma HLS PIPELINE II=1
        sol[i] = 0;
    }

    MAP_SOL: for (int i = 0; i < N; i++) {
        #pragma HLS PIPELINE II=1
        int idx = used_cols[i]; // Para no acceder dos veces
        if (idx >= 0 && idx < M) {
            sol[idx] = result[i];
        }
    }

    return 0;
}

// ÁLGEBRA LINEAL **************************************************************************************************/
/*void multiplicarMatrizVector(const bit_t synd[N], const bit_t H_s[N][N], bit_t result[N]) {
    MULT_MAT_VEC_OUTER: for (int i = 0; i < N; i++) {
        int sum = 0;
        MULT_MAT_VEC_INNER: for (int j = 0; j < N; j++) {
            sum ^= H_s[i][j] & synd[j];
        }
        result[i] = sum;
    }
}


void transponer_matriz_rect(const bit_t entrada[N][M], bit_t salida[M][N]) {
    TRANS_RECT_ROW: for (int i = 0; i < N; i++) {
        TRANS_RECT_COL: for (int j = 0; j < M; j++) {
            salida[j][i] = entrada[i][j];
        }
    }
}

void transponer_matriz(bit_t entrada[N][N], bit_t salida[N][N]) {
    TRANS_MAT_ROW: for (int i = 0; i < N; i++) {
        TRANS_MAT_COL: for (int j = 0; j < N; j++) {
            salida[j][i] = entrada[i][j];
        }
    }
}*/