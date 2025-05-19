#ifndef __MATRIX_IMAGEC_H__
#define __MATRIX_IMAGEC_H__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "alloc.h"

#define MAT_AT(mtx, x, y) ((mtx).data[((mtx).r * (x) + (y))])
#define M_PI 3.14159265358979323846

typedef struct matrix matrix_t;
typedef struct matrix {
    int r;
    int c;
    long double * data;
} matrix;

void print_mat(matrix_t mat);
void free_matrix(matrix_t * matrix);
matrix_t mat(int r, int c, long double * data);
matrix_t create_matrix(int r, int c);
matrix_t gaussian_kernel(int size, double sigma);
matrix_t mean_kernel(int size);
matrix_t sobel_x_kernel();
matrix_t sobel_y_kernel();
matrix_t sobel_45_kernel();
matrix_t sobel_135_kernel();
matrix_t sobel_x_kernel_x5();
matrix_t sobel_y_kernel_x5();
matrix_t sobel_45_kernel_x5();
matrix_t sobel_135_kernel_x5();
matrix_t laplacian_kernel();
matrix_t laplacian_4c_kernel();
#endif /* __MATRIX_IMAGEC_H__ */