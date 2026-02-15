#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include "matrix_lib.h"

// Надёжный ввод пункта меню через fgets/strtol
static int read_int_line(const char *prompt, int *out)
{
    char buf[128];
    for (;;)
    {
        if (prompt) printf("%s", prompt);
        if (!fgets(buf, sizeof(buf), stdin)) return 0;

        char *end = NULL;
        errno = 0;
        long v = strtol(buf, &end, 10);
        if (end == buf || errno != 0)
        {
            fprintf(stderr, "Некорректный ввод. Введите целое число.\n");
            continue;
        }
        *out = (int)v;
        return 1;
    }
}

static int read_double_line(const char *prompt, double *out)
{
    char buf[128];
    for (;;)
    {
        if (prompt) printf("%s", prompt);
        if (!fgets(buf, sizeof(buf), stdin)) return 0;

        char *end = NULL;
        errno = 0;
        double v = strtod(buf, &end);
        if (end == buf || errno != 0)
        {
            fprintf(stderr, "Некорректный ввод. Введите число.\n");
            continue;
        }
        *out = v;
        return 1;
    }
}

// съедаем хвост после scanf из твоих функций ввода матриц
static void eat_stdin_line(void)
{
    int c;
    while ((c = getchar()) != '\n' && c != EOF) {}
}

// критичные ошибки -> exit, а SIZE_MISMATCH/SINGULAR -> просто warn
static void handle_op_error(enum ERROS_TYPE e, const char *where,
                            const matrix *A, const matrix *B)
{
    if (e == OK) return;

    if (e == SIZE_MISMATCH_ERROR)
    {
        fprintf(stderr, "[WARN] %s: размеры не подходят.\n", where);
        if (A && B)
            fprintf(stderr, "       A=%zux%zu, B=%zux%zu\n", A->rows, A->colums, B->rows, B->colums);
        return;
    }

    if (e == SINGULAR_MATRIX_ERROR)
    {
        fprintf(stderr, "[WARN] %s: матрица вырождена (обратной нет).\n", where);
        return;
    }

    matrix_fatal(e, where);
}

static void input_and_fill(matrix *m, const char *name)
{
    printf("\n=== Ввод матрицы %s ===\n", name);
    enum INPUT_CHOISE ch = input_rows_colums(m);

    enum ERROS_TYPE e = type_matrix(m);
    if (e != OK) matrix_fatal(e, "type_matrix");

    if (ch == RANDOM) e = fill_random(m);
    else             e = fill_hands(m);

    if (e != OK) matrix_fatal(e, "fill");
}

int main(void)
{
    srand((unsigned)time(NULL));

    matrix A, B, R;
    constructor(&A);
    constructor(&B);
    constructor(&R);

    input_and_fill(&A, "A");
    input_and_fill(&B, "B");

    // после scanf-ов останется '\n'
    eat_stdin_line();

    for (;;)
    {
        printf("\n================= МЕНЮ =================\n");
        printf("1) Показать A\n");
        printf("2) Показать B\n");
        printf("3) Ввести A заново\n");
        printf("4) Ввести B заново\n");
        printf("5) A + B\n");
        printf("6) A - B\n");
        printf("7) A * B\n");
        printf("8) A * k\n");
        printf("9) det(A) (Гаусс)\n");
        printf("10) rank(A) (Гаусс)\n");
        printf("11) inv(A) (Гаусс–Жордан)\n");
        printf("12) inv(A) (adj(A)/det(A))\n");
        printf("0) Выход\n");

        int cmd = -1;
        if (!read_int_line("Ваш выбор: ", &cmd)) break;

        enum ERROS_TYPE e = OK;

        switch (cmd)
        {
            case 0:
                destructor(&A);
                destructor(&B);
                destructor(&R);
                return 0;

            case 1:
                printf("\nA (%zux%zu):\n", A.rows, A.colums);
                print_matrix(&A);
                break;

            case 2:
                printf("\nB (%zux%zu):\n", B.rows, B.colums);
                print_matrix(&B);
                break;

            case 3:
                free_matrix(&A);
                input_and_fill(&A, "A");
                eat_stdin_line();
                break;

            case 4:
                free_matrix(&B);
                input_and_fill(&B, "B");
                eat_stdin_line();
                break;

            case 5:
                e = matrix_add(&A, &B, &R);
                handle_op_error(e, "matrix_add(A,B)", &A, &B);
                if (e == OK) { printf("\nA + B:\n"); print_matrix(&R); }
                break;

            case 6:
                e = matrix_sub(&A, &B, &R);
                handle_op_error(e, "matrix_sub(A,B)", &A, &B);
                if (e == OK) { printf("\nA - B:\n"); print_matrix(&R); }
                break;

            case 7:
                e = matrix_mul(&A, &B, &R);
                handle_op_error(e, "matrix_mul(A,B)", &A, &B);
                if (e == OK) { printf("\nA * B:\n"); print_matrix(&R); }
                break;

            case 8:
            {
                double k = 0.0;
                if (!read_double_line("Введите k: ", &k)) break;
                e = matrix_scale(&A, k, &R);
                handle_op_error(e, "matrix_scale(A,k)", &A, NULL);
                if (e == OK) { printf("\nA * %.3lf:\n", k); print_matrix(&R); }
                break;
            }

            case 9:
            {
                double det = 0.0;
                e = matrix_det_gauss(&A, &det);
                handle_op_error(e, "matrix_det_gauss(A)", &A, NULL);
                if (e == OK) printf("det(A) = %.10lf\n", det);
                break;
            }

            case 10:
            {
                size_t rank = 0;
                e = matrix_rank_gauss(&A, &rank);
                handle_op_error(e, "matrix_rank_gauss(A)", &A, NULL);
                if (e == OK) printf("rank(A) = %zu\n", rank);
                break;
            }

            case 11:
                e = matrix_inv_gauss_jordan(&A, &R);
                handle_op_error(e, "matrix_inv_gauss_jordan(A)", &A, NULL);
                if (e == OK) { printf("\ninv(A) (Gauss–Jordan):\n"); print_matrix(&R); }
                break;

            case 12:
                e = matrix_inv_adjugate(&A, &R);
                handle_op_error(e, "matrix_inv_adjugate(A)", &A, NULL);
                if (e == OK) { printf("\ninv(A) (adj/det):\n"); print_matrix(&R); }
                break;

            default:
                fprintf(stderr, "Неизвестная команда.\n");
                break;
        }
    }

    destructor(&A);
    destructor(&B);
    destructor(&R);
    return 0;
}
