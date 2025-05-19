#include <stdio.h>

#include "matrix.h"

void print_mat(matrix_t mat){
    printf("+- Matrix ");
    int max = mat.c * 9 - 9;
    for(int i = 0; i < max; ++i) putc('-', stdout);
    putc('\n', stdout);
    for(int i = 0; i < mat.r; ++i){
        for(int j = 0; j < mat.c; ++j){
            printf("%*.3Lf ", 8, MAT_AT(mat, i, j));
        }
        printf("\n");
    }
}

matrix_t mat(int r, int c, long double * data){
    return create_matrix(r, c);
}

matrix_t create_matrix(int r, int c){
    matrix_t m = {
        .c = c,
        .r = r,
    };
    m.data = ALLOC(sizeof(*m.data)*c*r);
    return m;
}

void free_matrix(matrix_t * matrix){
    if(!matrix) return;
    matrix->c = 0;
    matrix->r = 0;
    FREE(matrix->data);
}

matrix_t gaussian_kernel(int size, double sigma){
    assert((size % 2) != 0);
    matrix_t mat = create_matrix(size, size);

    int center = size/2;
    double dx = 0, dy = 0;
    double sig2 = sigma*sigma;
    double sig2x = 2.0 * sig2;
    double norm = 1.0/(2.0 * sig2 * M_PI);

    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            dx = j - center;
            dy = i - center;
            MAT_AT(mat, i, j) = norm * exp(-(((dx*dx)+(dy*dy))/sig2x));
        }
    }

    return mat;
}

matrix_t mean_kernel(int size){
    assert((size % 2) != 0);
    matrix_t mat = create_matrix(size, size);
    double value = 1.0 / (double) (size*size);
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            MAT_AT(mat, i, j) = value;
        }
    }

    return mat;
}

matrix_t sobel_x_kernel_x5(){
    matrix_t mat = create_matrix(5, 5);
    long double sobel_values[] = {
        -1, -2, 0, 2, 1,   
        -1, -2, 0, 2, 1,  
        -2, -4, 0, 4, 2,  
        -1, -2, 0, 2, 1,  
        -1, -2, 0, 2, 1,      
    };
    for (int i = 0; i < 25; ++i) {
        mat.data[i] = sobel_values[i];
    }
    return mat;
}

matrix_t sobel_y_kernel_x5() {
    matrix_t mat = create_matrix(5, 5);
    long double sobel_values[] = {
        -1, -1, -2, -1, -1,
        -2, -2, -4, -2, -2, 
         2,  2,  4,  2,  2, 
         1,  1,  2,  1,  1,
         0,  0,  0,  0,  0,
    };
    for (int i = 0; i < 25; ++i) {
        mat.data[i] = sobel_values[i];
    }
    return mat;
}

matrix_t sobel_y_kernel() {
    matrix_t mat = create_matrix(3, 3);
    long double sobel_values[] = {
        -1,  0,  1,
        -2,  0,  2,
        -1,  0,  1
    };
    for (int i = 0; i < 9; ++i) {
        mat.data[i] = sobel_values[i];
    }
    return mat;
}

matrix_t sobel_x_kernel() {
    matrix_t mat = create_matrix(3, 3);
    long double sobel_values[] = {
        -1, -2, -1,
         0,  0,  0,
         1,  2,  1
    };
    for (int i = 0; i < 9; ++i) {
        mat.data[i] = sobel_values[i];
    }
    return mat;
}

matrix_t sobel_45_kernel_x5() {
    matrix_t mat = create_matrix(5, 5);
    long double sobel_values[] = {
        0,   -1,  -2,  -1,  0,
        1,    0,  -1,  -2,  -1,
        2,    1,   0,  -1,  -2,
        1,    2,   1,   0,  -1,
        0,    1,   2,   1,   0
    };
    for (int i = 0; i < 25; ++i) {
        mat.data[i] = sobel_values[i];
    }
    return mat;
}

matrix_t sobel_45_kernel() {
    matrix_t mat = create_matrix(3, 3);
    long double sobel_values[] = {
         0, -1, -2,
         1,  0, -1,
         2,  1,  0
    };
    for (int i = 0; i < 9; ++i) {
        mat.data[i] = sobel_values[i];
    }
    return mat;
}

matrix_t sobel_135_kernel_x5() {
    matrix_t mat = create_matrix(5, 5);
    long double sobel_values[] = {
         0,    1,   2,   1,   0,
        -1,    0,   1,   2,   1,
        -2,   -1,   0,   1,   2,
        -1,   -2,  -1,   0,   1,
         0,   -1,  -2,  -1,   0
    };
    for (int i = 0; i < 25; ++i) {
        mat.data[i] = sobel_values[i];
    }
    return mat;
}

matrix_t sobel_135_kernel() {
    matrix_t mat = create_matrix(3, 3);
    long double sobel_values[] = {
         2,  1,  0,
         1,  0, -1,
         0, -1, -2
    };
    for (int i = 0; i < 9; ++i) {
        mat.data[i] = sobel_values[i];
    }
    return mat;
}

matrix_t laplacian_4c_kernel() {
    matrix_t mat = create_matrix(3, 3);
    long double laplacian_values[] = {
         0, -1,  0,
        -1,  4, -1,
         0, -1,  0
    };
    for (int i = 0; i < 9; ++i) {
        mat.data[i] = laplacian_values[i];
    }
    return mat;
}

matrix_t laplacian_kernel() {
    matrix_t mat = create_matrix(3, 3);
    long double laplacian_values[] = {
        -1, -1,  -1,
        -1,  8,  -1,
        -1, -1,  -1
    };
    for (int i = 0; i < 9; ++i) {
        mat.data[i] = laplacian_values[i];
    }
    return mat;
}