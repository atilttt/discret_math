#include "matrix_lib.h"
#include <math.h>
#include <string.h>

//==================== Внутренние хелперы ====================//

static enum ERROS_TYPE validate_matrix_ro(const matrix *m)
{
    if (!m) return NULL_POINTER_ERROR;
    if (m->rows == 0 || m->colums == 0) return INPUT_OUTPUT_ERROR;
    if (!m->matrix) return NULL_POINTER_ERROR;

    for (size_t i = 0; i < m->rows; i++)
        if (!m->matrix[i]) return NULL_POINTER_ERROR;

    return OK;
}

// освобождает ТОЛЬКО данные, НЕ трогает rows/colums
static void free_matrix_data(matrix *m)
{
    if (!m || !m->matrix) return;

    for (size_t i = 0; i < m->rows; i++)
        free(m->matrix[i]);

    free(m->matrix);
    m->matrix = NULL;
}

static enum ERROS_TYPE alloc_matrix(matrix *m, size_t rows, size_t colums)
{
    if (!m) return NULL_POINTER_ERROR;
    if (rows == 0 || colums == 0) return INPUT_OUTPUT_ERROR;

    free_matrix_data(m);

    m->rows = rows;
    m->colums = colums;

    m->matrix = (double**)calloc(rows, sizeof(double*));
    if (!m->matrix) return MEMORY_ALLOCATION_ERROR;

    for (size_t i = 0; i < rows; i++)
    {
        m->matrix[i] = (double*)calloc(colums, sizeof(double));
        if (!m->matrix[i])
        {
            for (size_t k = 0; k < i; k++) free(m->matrix[k]);
            free(m->matrix);
            m->matrix = NULL;
            return MEMORY_ALLOCATION_ERROR;
        }
    }

    return OK;
}

static void swap_rows(double **m, size_t r1, size_t r2)
{
    double *tmp = m[r1];
    m[r1] = m[r2];
    m[r2] = tmp;
}

static enum ERROS_TYPE matrix_identity(matrix *out, size_t n)
{
    if (!out) return NULL_POINTER_ERROR;

    enum ERROS_TYPE e = alloc_matrix(out, n, n);
    if (e != OK) { out->errors = e; return e; }

    for (size_t i = 0; i < n; i++)
        out->matrix[i][i] = 1.0;

    out->errors = OK;
    return OK;
}

//==================== База ====================//

void constructor(matrix *pointer_matrix)
{
    assert(pointer_matrix != NULL && "Ошибка, пришел нулевой указатель\n");

    pointer_matrix->colums = 0;
    pointer_matrix->rows = 0;
    pointer_matrix->matrix = NULL;
    pointer_matrix->errors = OK;
}

void free_matrix(matrix *pointer_matrix)
{
    if (!pointer_matrix) return;

    free_matrix_data(pointer_matrix);

    pointer_matrix->rows = 0;
    pointer_matrix->colums = 0;
}

void destructor(matrix *pointer_matrix)
{
    assert(pointer_matrix != NULL && "Ошибка, пришел нулевой указатель\n");
    free_matrix(pointer_matrix);
    pointer_matrix->errors = OK;
}

void print_matrix(const matrix *pointer_matrix)
{
    enum ERROS_TYPE e = validate_matrix_ro(pointer_matrix);
    if (e != OK)
    {
        fprintf(stderr, "print_matrix: матрица невалидна, ошибка=%d\n", e);
        return;
    }

    for (size_t i = 0; i < pointer_matrix->rows; i++)
    {
        for (size_t j = 0; j < pointer_matrix->colums; j++)
            printf("%10.3lf ", pointer_matrix->matrix[i][j]);
        printf("\n");
    }
}

enum ERROS_TYPE copy_matrix(const matrix *src, matrix *dst)
{
    if (!dst) return NULL_POINTER_ERROR;

    enum ERROS_TYPE e = validate_matrix_ro(src);
    if (e != OK)
    {
        dst->errors = e;
        return e;
    }

    e = alloc_matrix(dst, src->rows, src->colums);
    if (e != OK)
    {
        dst->errors = e;
        return e;
    }

    for (size_t i = 0; i < src->rows; i++)
        for (size_t j = 0; j < src->colums; j++)
            dst->matrix[i][j] = src->matrix[i][j];

    dst->errors = OK;
    return OK;
}

//==================== Ввод/создание матрицы ====================//

enum INPUT_CHOISE input_rows_colums(matrix *pointer_matrix)
{
    assert(pointer_matrix != NULL && "Ошибк, пришел нулевой указатель\n");

    printf("--------------------------------------------------Заполнение кол-ва строк и столбцов--------------------------------------------------\n");
    printf("Вам нужно ввести количество строк и количество столбцов в вашей матрице\n");

    printf("Введите количество строк: ");
    if (scanf("%zu", &pointer_matrix->rows) != 1 || pointer_matrix->rows == 0)
    {
        fprintf(stderr, "Вы ввели не число (нужно целое >= 1)\n");
        pointer_matrix->errors = INPUT_OUTPUT_ERROR;
        exit(INPUT_OUTPUT_ERROR);
    }

    printf("Введите количество столбцов: ");
    if (scanf("%zu", &pointer_matrix->colums) != 1 || pointer_matrix->colums == 0)
    {
        fprintf(stderr, "Вы ввели не число (нужно целое >= 1)\n");
        pointer_matrix->errors = INPUT_OUTPUT_ERROR;
        exit(INPUT_OUTPUT_ERROR);
    }

    printf("Вы ввели [%zu] строк и [%zu] столбцов\n", pointer_matrix->rows, pointer_matrix->colums);

    printf("Вам нужно выбрать, каким способом вы собираетесь заполнять матрицу\n");
    printf("Ваш выбор %d - Рандомное заполнение или %d - Ручной ввод\n", RANDOM, HADNS_INPUT);

    int choise = 0;
    if (scanf("%d", &choise) != 1 || (choise != RANDOM && choise != HADNS_INPUT))
    {
        fprintf(stderr, "Вы ввели не число или число не соответствует ни одному из вариантов\n");
        pointer_matrix->errors = INPUT_OUTPUT_ERROR;
        exit(INPUT_OUTPUT_ERROR);
    }

    printf("Вы выбрали способ заполнения матрицы: %s\n",
           (choise == RANDOM) ? "Рандомное заполнение" : "Ручной ввод");

    printf("--------------------------------------------------Переходим к заполнению матрицы--------------------------------------------------\n");
    pointer_matrix->errors = OK;

    return (enum INPUT_CHOISE)choise;
}

enum ERROS_TYPE type_matrix(matrix *pointer_matrix)
{
    if (!pointer_matrix) return NULL_POINTER_ERROR;

    // ВАЖНО: не сбрасываем rows/colums, только память чистим
    free_matrix_data(pointer_matrix);

    if (pointer_matrix->rows == 0 || pointer_matrix->colums == 0)
    {
        pointer_matrix->errors = INPUT_OUTPUT_ERROR;
        return INPUT_OUTPUT_ERROR;
    }

    enum ERROS_TYPE e = alloc_matrix(pointer_matrix, pointer_matrix->rows, pointer_matrix->colums);
    pointer_matrix->errors = e;
    return e;
}

enum ERROS_TYPE fill_random(matrix *pointer_matrix)
{
    if (!pointer_matrix) return NULL_POINTER_ERROR;

    enum ERROS_TYPE e = validate_matrix_ro(pointer_matrix);
    if (e != OK)
    {
        pointer_matrix->errors = e;
        return e;
    }

    printf("Рандомное заполнение (rand)\n");
    printf("Введите min и max (min < max): ");

    int min = 0, max = 0;
    if (scanf("%d", &min) != 1 || scanf("%d", &max) != 1 || min >= max)
    {
        fprintf(stderr, "Некорректный ввод (min должно быть меньше max)\n");
        pointer_matrix->errors = INPUT_OUTPUT_ERROR;
        return INPUT_OUTPUT_ERROR;
    }

    const int range = max - min + 1;

    for (size_t i = 0; i < pointer_matrix->rows; i++)
        for (size_t j = 0; j < pointer_matrix->colums; j++)
            pointer_matrix->matrix[i][j] = (double)(min + rand() % range);

    pointer_matrix->errors = OK;
    return OK;
}

enum ERROS_TYPE fill_hands(matrix *pointer_matrix)
{
    if (!pointer_matrix) return NULL_POINTER_ERROR;

    enum ERROS_TYPE e = validate_matrix_ro(pointer_matrix);
    if (e != OK)
    {
        pointer_matrix->errors = e;
        return e;
    }

    printf("Ручной ввод матрицы (%zu x %zu)\n", pointer_matrix->rows, pointer_matrix->colums);

    for (size_t i = 0; i < pointer_matrix->rows; i++)
    {
        for (size_t j = 0; j < pointer_matrix->colums; j++)
        {
            printf("a[%zu][%zu] = ", i, j);
            if (scanf("%lf", &pointer_matrix->matrix[i][j]) != 1)
            {
                fprintf(stderr, "Некорректный ввод элемента матрицы\n");
                pointer_matrix->errors = INPUT_OUTPUT_ERROR;
                return INPUT_OUTPUT_ERROR;
            }
        }
    }

    pointer_matrix->errors = OK;
    return OK;
}

//==================== Операции ====================//

enum ERROS_TYPE matrix_add(const matrix *a, const matrix *b, matrix *out)
{
    if (!out) return NULL_POINTER_ERROR;

    enum ERROS_TYPE ea = validate_matrix_ro(a);
    enum ERROS_TYPE eb = validate_matrix_ro(b);
    if (ea != OK) { out->errors = ea; return ea; }
    if (eb != OK) { out->errors = eb; return eb; }

    if (a->rows != b->rows || a->colums != b->colums)
    {
        out->errors = SIZE_MISMATCH_ERROR;
        return SIZE_MISMATCH_ERROR;
    }

    enum ERROS_TYPE e = alloc_matrix(out, a->rows, a->colums);
    if (e != OK) { out->errors = e; return e; }

    for (size_t i = 0; i < a->rows; i++)
        for (size_t j = 0; j < a->colums; j++)
            out->matrix[i][j] = a->matrix[i][j] + b->matrix[i][j];

    out->errors = OK;
    return OK;
}

enum ERROS_TYPE matrix_sub(const matrix *a, const matrix *b, matrix *out)
{
    if (!out) return NULL_POINTER_ERROR;

    enum ERROS_TYPE ea = validate_matrix_ro(a);
    enum ERROS_TYPE eb = validate_matrix_ro(b);
    if (ea != OK) { out->errors = ea; return ea; }
    if (eb != OK) { out->errors = eb; return eb; }

    if (a->rows != b->rows || a->colums != b->colums)
    {
        out->errors = SIZE_MISMATCH_ERROR;
        return SIZE_MISMATCH_ERROR;
    }

    enum ERROS_TYPE e = alloc_matrix(out, a->rows, a->colums);
    if (e != OK) { out->errors = e; return e; }

    for (size_t i = 0; i < a->rows; i++)
        for (size_t j = 0; j < a->colums; j++)
            out->matrix[i][j] = a->matrix[i][j] - b->matrix[i][j];

    out->errors = OK;
    return OK;
}

enum ERROS_TYPE matrix_mul(const matrix *a, const matrix *b, matrix *out)
{
    if (!out) return NULL_POINTER_ERROR;

    enum ERROS_TYPE ea = validate_matrix_ro(a);
    enum ERROS_TYPE eb = validate_matrix_ro(b);
    if (ea != OK) { out->errors = ea; return ea; }
    if (eb != OK) { out->errors = eb; return eb; }

    if (a->colums != b->rows)
    {
        out->errors = SIZE_MISMATCH_ERROR;
        return SIZE_MISMATCH_ERROR;
    }

    enum ERROS_TYPE e = alloc_matrix(out, a->rows, b->colums);
    if (e != OK) { out->errors = e; return e; }

    for (size_t i = 0; i < a->rows; i++)
    {
        for (size_t j = 0; j < b->colums; j++)
        {
            double sum = 0.0;
            for (size_t k = 0; k < a->colums; k++)
                sum += a->matrix[i][k] * b->matrix[k][j];

            out->matrix[i][j] = sum;
        }
    }

    out->errors = OK;
    return OK;
}

enum ERROS_TYPE matrix_scale(const matrix *a, double k, matrix *out)
{
    if (!out) return NULL_POINTER_ERROR;

    enum ERROS_TYPE ea = validate_matrix_ro(a);
    if (ea != OK) { out->errors = ea; return ea; }

    enum ERROS_TYPE e = alloc_matrix(out, a->rows, a->colums);
    if (e != OK) { out->errors = e; return e; }

    for (size_t i = 0; i < a->rows; i++)
        for (size_t j = 0; j < a->colums; j++)
            out->matrix[i][j] = a->matrix[i][j] * k;

    out->errors = OK;
    return OK;
}

//==================== det (Гаусс) ====================//

enum ERROS_TYPE matrix_det_gauss(const matrix *a, double *det_out)
{
    if (!det_out) return NULL_POINTER_ERROR;

    enum ERROS_TYPE ea = validate_matrix_ro(a);
    if (ea != OK) return ea;

    if (a->rows != a->colums) return SIZE_MISMATCH_ERROR;

    const double EPS = 1e-12;
    const size_t n = a->rows;

    matrix tmp;
    constructor(&tmp);
    enum ERROS_TYPE e = copy_matrix(a, &tmp);
    if (e != OK) { destructor(&tmp); return e; }

    int sign = 1;

    for (size_t col = 0; col < n; col++)
    {
        size_t pivot = col;
        double best = fabs(tmp.matrix[col][col]);
        for (size_t r = col + 1; r < n; r++)
        {
            double v = fabs(tmp.matrix[r][col]);
            if (v > best) { best = v; pivot = r; }
        }

        if (best < EPS)
        {
            *det_out = 0.0;
            destructor(&tmp);
            return OK;
        }

        if (pivot != col)
        {
            swap_rows(tmp.matrix, pivot, col);
            sign = -sign;
        }

        double piv = tmp.matrix[col][col];
        for (size_t r = col + 1; r < n; r++)
        {
            double factor = tmp.matrix[r][col] / piv;
            tmp.matrix[r][col] = 0.0;
            for (size_t c = col + 1; c < n; c++)
                tmp.matrix[r][c] -= factor * tmp.matrix[col][c];
        }
    }

    double det = (double)sign;
    for (size_t i = 0; i < n; i++)
        det *= tmp.matrix[i][i];

    *det_out = det;
    destructor(&tmp);
    return OK;
}

//==================== rank (Гаусс) ====================//

enum ERROS_TYPE matrix_rank_gauss(const matrix *a, size_t *rank_out)
{
    if (!rank_out) return NULL_POINTER_ERROR;

    enum ERROS_TYPE ea = validate_matrix_ro(a);
    if (ea != OK) return ea;

    const double EPS = 1e-12;

    matrix tmp;
    constructor(&tmp);
    enum ERROS_TYPE e = copy_matrix(a, &tmp);
    if (e != OK) { destructor(&tmp); return e; }

    const size_t m = tmp.rows;
    const size_t n = tmp.colums;

    size_t rank = 0;
    size_t row = 0;

    for (size_t col = 0; col < n && row < m; col++)
    {
        size_t pivot = row;
        double best = fabs(tmp.matrix[row][col]);
        for (size_t r = row + 1; r < m; r++)
        {
            double v = fabs(tmp.matrix[r][col]);
            if (v > best) { best = v; pivot = r; }
        }

        if (best < EPS)
            continue;

        if (pivot != row)
            swap_rows(tmp.matrix, pivot, row);

        double piv = tmp.matrix[row][col];
        for (size_t r = row + 1; r < m; r++)
        {
            double factor = tmp.matrix[r][col] / piv;
            tmp.matrix[r][col] = 0.0;
            for (size_t c = col + 1; c < n; c++)
                tmp.matrix[r][c] -= factor * tmp.matrix[row][c];
        }

        rank++;
        row++;
    }

    *rank_out = rank;
    destructor(&tmp);
    return OK;
}

//==================== inverse 1: Gauss–Jordan ====================//

enum ERROS_TYPE matrix_inv_gauss_jordan(const matrix *a, matrix *out)
{
    if (!out) return NULL_POINTER_ERROR;

    enum ERROS_TYPE ea = validate_matrix_ro(a);
    if (ea != OK) { out->errors = ea; return ea; }

    if (a->rows != a->colums)
    {
        out->errors = SIZE_MISMATCH_ERROR;
        return SIZE_MISMATCH_ERROR;
    }

    const double EPS = 1e-12;
    const size_t n = a->rows;

    matrix tmp, inv;
    constructor(&tmp);
    constructor(&inv);

    enum ERROS_TYPE e = copy_matrix(a, &tmp);
    if (e != OK) { destructor(&tmp); out->errors = e; return e; }

    e = matrix_identity(&inv, n);
    if (e != OK) { destructor(&tmp); destructor(&inv); out->errors = e; return e; }

    for (size_t col = 0; col < n; col++)
    {
        size_t pivot = col;
        double best = fabs(tmp.matrix[col][col]);
        for (size_t r = col + 1; r < n; r++)
        {
            double v = fabs(tmp.matrix[r][col]);
            if (v > best) { best = v; pivot = r; }
        }

        if (best < EPS)
        {
            destructor(&tmp);
            destructor(&inv);
            out->errors = SINGULAR_MATRIX_ERROR;
            return SINGULAR_MATRIX_ERROR;
        }

        if (pivot != col)
        {
            swap_rows(tmp.matrix, pivot, col);
            swap_rows(inv.matrix, pivot, col);
        }

        double piv = tmp.matrix[col][col];
        for (size_t c = 0; c < n; c++)
        {
            tmp.matrix[col][c] /= piv;
            inv.matrix[col][c] /= piv;
        }

        for (size_t r = 0; r < n; r++)
        {
            if (r == col) continue;
            double factor = tmp.matrix[r][col];
            if (fabs(factor) < EPS) { tmp.matrix[r][col] = 0.0; continue; }

            tmp.matrix[r][col] = 0.0;
            for (size_t c = 0; c < n; c++)
            {
                tmp.matrix[r][c] -= factor * tmp.matrix[col][c];
                inv.matrix[r][c] -= factor * inv.matrix[col][c];
            }
        }
    }

    e = copy_matrix(&inv, out);
    destructor(&tmp);
    destructor(&inv);

    out->errors = e;
    return e;
}

//==================== inverse 2: adj(A)/det(A) ====================//

static enum ERROS_TYPE minor_det(const matrix *a, size_t skip_r, size_t skip_c, double *det_out)
{
    const size_t n = a->rows;

    matrix m;
    constructor(&m);

    enum ERROS_TYPE e = alloc_matrix(&m, n - 1, n - 1);
    if (e != OK) { destructor(&m); return e; }

    size_t rr = 0;
    for (size_t i = 0; i < n; i++)
    {
        if (i == skip_r) continue;
        size_t cc = 0;
        for (size_t j = 0; j < n; j++)
        {
            if (j == skip_c) continue;
            m.matrix[rr][cc] = a->matrix[i][j];
            cc++;
        }
        rr++;
    }

    e = matrix_det_gauss(&m, det_out);
    destructor(&m);
    return e;
}

enum ERROS_TYPE matrix_inv_adjugate(const matrix *a, matrix *out)
{
    if (!out) return NULL_POINTER_ERROR;

    enum ERROS_TYPE ea = validate_matrix_ro(a);
    if (ea != OK) { out->errors = ea; return ea; }

    if (a->rows != a->colums)
    {
        out->errors = SIZE_MISMATCH_ERROR;
        return SIZE_MISMATCH_ERROR;
    }

    const double EPS = 1e-12;
    const size_t n = a->rows;

    double det = 0.0;
    enum ERROS_TYPE e = matrix_det_gauss(a, &det);
    if (e != OK) { out->errors = e; return e; }

    if (fabs(det) < EPS)
    {
        out->errors = SINGULAR_MATRIX_ERROR;
        return SINGULAR_MATRIX_ERROR;
    }

    matrix C;
    constructor(&C);

    e = alloc_matrix(&C, n, n);
    if (e != OK) { destructor(&C); out->errors = e; return e; }

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            double md = 0.0;
            e = minor_det(a, i, j, &md);
            if (e != OK)
            {
                destructor(&C);
                out->errors = e;
                return e;
            }

            C.matrix[i][j] = (((i + j) % 2) ? -md : md);
        }
    }

    e = alloc_matrix(out, n, n);
    if (e != OK)
    {
        destructor(&C);
        out->errors = e;
        return e;
    }

    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            out->matrix[i][j] = C.matrix[j][i] / det;

    destructor(&C);
    out->errors = OK;
    return OK;
}

//==================== Ошибки: печать + exit ====================//

void matrix_fatal(enum ERROS_TYPE err, const char *where)
{
    if (err == OK) return;
    if (!where) where = "matrix_lib";

    fprintf(stderr, "\n[FATAL] %s: ", where);

    switch (err)
    {
        case INPUT_OUTPUT_ERROR:
            fprintf(stderr, "ошибка ввода/вывода (некорректные данные)\n");
            break;

        case MEMORY_ALLOCATION_ERROR:
            fprintf(stderr, "ошибка выделения памяти\n");
            fprintf(stderr, "        errno=%d (%s)\n", errno, strerror(errno));
            break;

        case NULL_POINTER_ERROR:
            fprintf(stderr, "нулевой указатель\n");
            break;

        case SIZE_MISMATCH_ERROR:
            fprintf(stderr, "несовпадение размеров матриц\n");
            break;

        case SINGULAR_MATRIX_ERROR:
            fprintf(stderr, "матрица вырождена (det = 0), обратной нет\n");
            break;

        default:
            fprintf(stderr, "неизвестная ошибка (%d)\n", (int)err);
            break;
    }

    exit((int)err);
}
