// slau.c
#include "SLAU.h"
#include <stdio.h>
#include <math.h>

// ==================== Локальные хелперы ==================== //

static double slau_get_eps(const slau *s)
{
    if (!s) return 1e-12;
    return (s->eps > 0.0) ? s->eps : 1e-12;
}

// out = [A | b], m x (n+1)
static enum ERROS_TYPE slau_build_augmented(const matrix *A, const matrix *b, matrix *out)
{
    if (!out) return NULL_POINTER_ERROR;

    enum ERROS_TYPE ea = validate_matrix_ro(A);
    enum ERROS_TYPE eb = validate_matrix_ro(b);
    if (ea != OK) return ea;
    if (eb != OK) return eb;

    if (b->colums != 1) return SIZE_MISMATCH_ERROR;
    if (A->rows != b->rows) return SIZE_MISMATCH_ERROR;

    enum ERROS_TYPE e = alloc_matrix(out, A->rows, A->colums + 1);
    if (e != OK) return e;

    for (size_t i = 0; i < A->rows; i++)
    {
        for (size_t j = 0; j < A->colums; j++)
            out->matrix[i][j] = A->matrix[i][j];

        out->matrix[i][A->colums] = b->matrix[i][0];
    }

    out->errors = OK;
    return OK;
}

// Проверка: строка "нулевая" по коэффициентам?
static int row_all_zero_coeffs(const matrix *Ab, size_t row, size_t n, double eps)
{
    for (size_t j = 0; j < n; j++)
        if (fabs(Ab->matrix[row][j]) > eps)
            return 0;
    return 1;
}

// ==================== Инфо ==================== //

void reminder(void)
{
    printf("\n==================== Напоминалка линала ====================\n");
    printf("Кронекер–Капелли:\n");
    printf("  rank(A) != rank([A|b]) -> несовместна\n");
    printf("  rank(A) == rank([A|b]) < n -> беск. много решений\n");
    printf("  rank(A) == rank([A|b]) == n -> единственное решение\n");
    printf("============================================================\n\n");
}

void instruction(void)
{
    printf("\n==================== Инструкция ====================\n");
    printf("Решаем методом Гаусса:\n");
    printf("  1) строим [A|b]\n");
    printf("  2) приводим к ступенчатому виду\n");
    printf("  3) проверяем несовместность\n");
    printf("  4) обратный ход -> x (для беск. мн-ва вернем частное: свободные=0)\n");
    printf("====================================================\n\n");
}

// ==================== Конструктор/деструктор ==================== //

void slau_constructor(slau *s)
{
    if (!s) return;

    s->choice = METHOD_NONE;
    s->eps = 1e-12;
    s->err = SLAU_OK;

    constructor(&s->A);
    constructor(&s->b);
    constructor(&s->x);
}

void slau_destructor(slau *s)
{
    if (!s) return;

    destructor(&s->A);
    destructor(&s->b);
    destructor(&s->x);

    s->choice = METHOD_NONE;
    s->eps = 1e-12;
    s->err = SLAU_OK;
}

OSHIBKA_NAFIK slau_validate_basic(const slau *s)
{
    if (!s) return SLAU_BAD_INPUT;

    if (validate_matrix_ro(&s->A) != OK) return SLAU_BAD_INPUT;
    if (validate_matrix_ro(&s->b) != OK) return SLAU_BAD_INPUT;

    if (s->b.colums != 1) return SLAU_DIM_MISMATCH;
    if (s->A.rows != s->b.rows) return SLAU_DIM_MISMATCH;

    return SLAU_OK;
}

OSHIBKA_NAFIK slau_check_kronecker_capelli(const slau *s)
{
    OSHIBKA_NAFIK vb = slau_validate_basic(s);
    if (vb != SLAU_OK) return vb;

    matrix Ab;
    constructor(&Ab);

    enum ERROS_TYPE e = slau_build_augmented(&s->A, &s->b, &Ab);
    if (e != OK) { destructor(&Ab); return SLAU_BAD_INPUT; }

    size_t rA = 0, rAb = 0;

    e = matrix_rank_gauss(&s->A, &rA);
    if (e != OK) { destructor(&Ab); return SLAU_BAD_INPUT; }

    e = matrix_rank_gauss(&Ab, &rAb);
    destructor(&Ab);
    if (e != OK) return SLAU_BAD_INPUT;

    const size_t n = s->A.colums;

    if (rA != rAb) return SLAU_INCOMPATIBLE;
    if (rA < n)    return SLAU_INFINITE_SOL;

    return SLAU_OK;
}

OSHIBKA_NAFIK slau_solve_gauss_jordan(slau *s)
{
    if (!s) return SLAU_BAD_INPUT;

    OSHIBKA_NAFIK vb = slau_validate_basic(s);
    if (vb != SLAU_OK) return (s->err = vb);

    const double EPS = slau_get_eps(s);
    const size_t m = s->A.rows;
    const size_t n = s->A.colums;

    matrix Ab;
    constructor(&Ab);
    enum ERROS_TYPE e = slau_build_augmented(&s->A, &s->b, &Ab);
    if (e != OK)
    {
        destructor(&Ab);
        return (s->err = SLAU_BAD_INPUT);
    }

    size_t *pivot_cols = (size_t*)calloc((m < n ? m : n), sizeof(size_t));
    if (!pivot_cols)
    {
        destructor(&Ab);
        return (s->err = SLAU_BAD_INPUT);
    }

    size_t rank = 0;     
    size_t row = 0;

    for (size_t col = 0; col < n && row < m; col++)
    {
        size_t pivot = row;
        double best = fabs(Ab.matrix[row][col]);
        for (size_t r = row + 1; r < m; r++)
        {
            double v = fabs(Ab.matrix[r][col]);
            if (v > best) { best = v; pivot = r; }
        }

        if (best <= EPS) continue;

        if (pivot != row)
            swap_rows(Ab.matrix, pivot, row);

        pivot_cols[rank] = col;

        const double piv = Ab.matrix[row][col];

        for (size_t r = row + 1; r < m; r++)
        {
            double a_rc = Ab.matrix[r][col];
            if (fabs(a_rc) <= EPS) { Ab.matrix[r][col] = 0.0; continue; }

            double factor = a_rc / piv;
            Ab.matrix[r][col] = 0.0;

            for (size_t c = col + 1; c <= n; c++) // <= n, т.к. последний столбец — b
                Ab.matrix[r][c] -= factor * Ab.matrix[row][c];
        }

        row++;
        rank++;
    }

    for (size_t r = 0; r < m; r++)
    {
        if (row_all_zero_coeffs(&Ab, r, n, EPS) && fabs(Ab.matrix[r][n]) > EPS)
        {
            free(pivot_cols);
            destructor(&Ab);
            return (s->err = SLAU_INCOMPATIBLE);
        }
    }

    destructor(&s->x);
    constructor(&s->x);
    e = alloc_matrix(&s->x, n, 1);
    if (e != OK)
    {
        free(pivot_cols);
        destructor(&Ab);
        return (s->err = SLAU_BAD_INPUT);
    }

    for (size_t i = 0; i < n; i++)
        s->x.matrix[i][0] = 0.0;

    for (size_t k = rank; k-- > 0; )
    {
        size_t r = k;             
        size_t col = pivot_cols[k]; 

        double rhs = Ab.matrix[r][n];

        // вычитаем известные x[j]
        for (size_t j = col + 1; j < n; j++)
            rhs -= Ab.matrix[r][j] * s->x.matrix[j][0];

        double coeff = Ab.matrix[r][col];
        if (fabs(coeff) <= EPS)
        {
            free(pivot_cols);
            destructor(&Ab);
            return (s->err = SLAU_BAD_INPUT);
        }

        s->x.matrix[col][0] = rhs / coeff;
    }

    free(pivot_cols);
    destructor(&Ab);

    if (rank < n)
        return (s->err = SLAU_INFINITE_SOL);  
    else
        return (s->err = SLAU_OK);            
}

OSHIBKA_NAFIK slau_solve_kramer(slau *s)
{
    if (!s) return SLAU_BAD_INPUT;

    OSHIBKA_NAFIK vb = slau_validate_basic(s);
    if (vb != SLAU_OK) return (s->err = vb);

    if (s->A.rows != s->A.colums)
        return (s->err = SLAU_METHOD_UNSUPPORTED);

    const double EPS = slau_get_eps(s);

    double detA = 0.0;
    enum ERROS_TYPE e = matrix_det_gauss(&s->A, &detA);
    if (e != OK)
        return (s->err = SLAU_BAD_INPUT);

    if (fabs(detA) < EPS)
        return (s->err = SLAU_SINGULAR);

    const size_t n = s->A.rows;

    destructor(&s->x);
    constructor(&s->x);
    e = alloc_matrix(&s->x, n, 1);
    if (e != OK)
        return (s->err = SLAU_BAD_INPUT);

    matrix Ai;
    constructor(&Ai);

    for (size_t col = 0; col < n; col++)
    {
        e = copy_matrix(&s->A, &Ai);
        if (e != OK)
        {
            destructor(&Ai);
            return (s->err = SLAU_BAD_INPUT);
        }

        for (size_t r = 0; r < n; r++)
            Ai.matrix[r][col] = s->b.matrix[r][0];

        double detAi = 0.0;
        e = matrix_det_gauss(&Ai, &detAi);

        destructor(&Ai);
        constructor(&Ai);

        if (e != OK)
        {
            destructor(&Ai);
            return (s->err = SLAU_BAD_INPUT);
        }

        s->x.matrix[col][0] = detAi / detA;
    }

    destructor(&Ai);
    return (s->err = SLAU_OK);
}

OSHIBKA_NAFIK slau_solve(slau *s)
{
    if (!s) return SLAU_BAD_INPUT;

    switch (s->choice)
    {
        case METHOD_GAUSS_JORDAN: // теперь это ГАУСС
            return slau_solve_gauss_jordan(s);

        case METHOD_KRAMER:
            return slau_solve_kramer(s);

        default:
            return (s->err = SLAU_BAD_INPUT);
    }
}