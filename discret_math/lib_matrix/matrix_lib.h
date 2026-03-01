#ifndef MATRIX_LIB_H
#define MATRIX_LIB_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>

enum INPUT_CHOISE
{
    RANDOM      = 1,
    HADNS_INPUT = 0   
};

enum ERROS_TYPE
{
    OK = 0,
    INPUT_OUTPUT_ERROR = 1,
    MEMORY_ALLOCATION_ERROR = 2,
    NULL_POINTER_ERROR = 3,
    SIZE_MISMATCH_ERROR = 4,
    SINGULAR_MATRIX_ERROR = 5
};

typedef struct
{
    size_t rows;
    size_t colums;
    double **matrix;
    enum ERROS_TYPE errors;
} matrix;

//-------------------- База --------------------//
void constructor(matrix *pointer_matrix);
void destructor(matrix *pointer_matrix);

void print_matrix(const matrix *pointer_matrix);
void free_matrix(matrix *pointer_matrix);
enum ERROS_TYPE copy_matrix(const matrix *src, matrix *dst);

//-------------------- Ввод/заполнение --------------------//
enum INPUT_CHOISE input_rows_colums(matrix *pointer_matrix);
enum ERROS_TYPE type_matrix(matrix *pointer_matrix);

enum ERROS_TYPE fill_random(matrix *pointer_matrix);
enum ERROS_TYPE fill_hands(matrix *pointer_matrix);

//-------------------- Операции --------------------//
enum ERROS_TYPE matrix_add(const matrix *a, const matrix *b, matrix *out);
enum ERROS_TYPE matrix_sub(const matrix *a, const matrix *b, matrix *out);
enum ERROS_TYPE matrix_mul(const matrix *a, const matrix *b, matrix *out);
enum ERROS_TYPE matrix_scale(const matrix *a, double k, matrix *out);

//-------------------- Линейная алгебра --------------------//
enum ERROS_TYPE matrix_det_gauss(const matrix *a, double *det_out);
enum ERROS_TYPE matrix_rank_gauss(const matrix *a, size_t *rank_out);

enum ERROS_TYPE matrix_inv_gauss_jordan(const matrix *a, matrix *out);
enum ERROS_TYPE matrix_inv_adjugate(const matrix *a, matrix *out);

//-------------------- Ошибки --------------------//
void matrix_fatal(enum ERROS_TYPE err, const char *where);
enum ERROS_TYPE validate_matrix_ro(const matrix *m);
void free_matrix_data(matrix *m);
enum ERROS_TYPE alloc_matrix(matrix *m, size_t rows, size_t colums);
void swap_rows(double **m, size_t r1, size_t r2);
enum ERROS_TYPE matrix_identity(matrix *out, size_t n);

#endif