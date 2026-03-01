// main.c
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <string.h>

#include "matrix_lib.h"
#include "SLAU.h"

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

static void eat_stdin_line(void)
{
    int c;
    while ((c = getchar()) != '\n' && c != EOF) {}
}

// -------------------- Общая обработка ошибок для матричных операций -------------------- //

static void handle_op_error(enum ERROS_TYPE e, const char *where,
                            const matrix *A, const matrix *B)
{
    if (e == OK) return;

    if (e == SIZE_MISMATCH_ERROR)
    {
        fprintf(stderr, "[WARN] %s: размеры не подходят.\n", where);
        if (A && B)
            fprintf(stderr, "       A=%zux%zu, B=%zux%zu\n",
                    A->rows, A->colums, B->rows, B->colums);
        return;
    }

    if (e == SINGULAR_MATRIX_ERROR)
    {
        fprintf(stderr, "[WARN] %s: матрица вырождена (обратной нет).\n", where);
        return;
    }

    matrix_fatal(e, where);
}

// -------------------- Ввод матрицы (использует твою матрикс-либу, которая внутри scanf) -------------------- //

static void input_and_fill(matrix *m, const char *name)
{
    printf("\n=== Ввод матрицы %s ===\n", name);

    enum INPUT_CHOISE ch = input_rows_colums(m);

    enum ERROS_TYPE e = type_matrix(m);
    if (e != OK) matrix_fatal(e, "type_matrix");

    if (ch == RANDOM) e = fill_random(m);
    else             e = fill_hands(m);

    if (e != OK) matrix_fatal(e, "fill");

    eat_stdin_line();
}

// -------------------- Ввод СЛАУ: A (m×n), b (m×1) через fgets -------------------- //
// Ввод построчно: a11 a12 ... a1n b1

static int parse_doubles_from_line(const char *line, double *arr, size_t need)
{
    size_t got = 0;
    const char *p = line;

    while (got < need)
    {
        while (*p == ' ' || *p == '\t') p++;

        if (*p == '\0' || *p == '\n') break;

        errno = 0;
        char *end = NULL;
        double v = strtod(p, &end);
        if (end == p || errno != 0) return 0;

        arr[got++] = v;
        p = end;
    }

    return (got == need);
}

static void slau_input_system(slau *S)
{
    printf("\n=== Ввод СЛАУ: A*x = b ===\n");

    int m = 0, n = 0;
    if (!read_int_line("Введите m (кол-во уравнений): ", &m) || m <= 0)
    {
        fprintf(stderr, "m должно быть >= 1\n");
        return;
    }
    if (!read_int_line("Введите n (кол-во неизвестных): ", &n) || n <= 0)
    {
        fprintf(stderr, "n должно быть >= 1\n");
        return;
    }

    // Пересоздаём A и b
    destructor(&S->A); constructor(&S->A);
    destructor(&S->b); constructor(&S->b);

    enum ERROS_TYPE e = alloc_matrix(&S->A, (size_t)m, (size_t)n);
    if (e != OK) matrix_fatal(e, "alloc_matrix(S.A)");

    e = alloc_matrix(&S->b, (size_t)m, 1);
    if (e != OK) matrix_fatal(e, "alloc_matrix(S.b)");

    printf("\nВвод: в каждой строке %d коэффициентов и затем свободный член.\n", n);
    printf("Пример (n=3):  1 2 3 10   значит: 1*x1 + 2*x2 + 3*x3 = 10\n\n");

    // буфер под числа строки
    double *tmp = (double*)malloc(sizeof(double) * (size_t)(n + 1));
    if (!tmp) matrix_fatal(MEMORY_ALLOCATION_ERROR, "malloc tmp");

    char buf[1024];

    for (int i = 0; i < m; i++)
    {
        for (;;)
        {
            printf("Строка %d: ", i + 1);
            if (!fgets(buf, sizeof(buf), stdin))
            {
                fprintf(stderr, "EOF при вводе.\n");
                free(tmp);
                return;
            }

            if (!parse_doubles_from_line(buf, tmp, (size_t)(n + 1)))
            {
                fprintf(stderr, "Нужно ввести ровно %d чисел: a1..a%d и b.\n", n + 1, n);
                continue;
            }

            for (int j = 0; j < n; j++)
                S->A.matrix[i][j] = tmp[j];

            S->b.matrix[i][0] = tmp[n];
            break;
        }
    }

    free(tmp);

    printf("\nA (%dx%d):\n", m, n);
    print_matrix(&S->A);
    printf("\nb (%dx1):\n", m);
    print_matrix(&S->b);
}

static int slau_choose_method(void)
{
    for (;;)
    {
        printf("\nВыберите метод решения СЛАУ:\n");
        printf("  1) Гаусс (ступенчатый вид + обратный ход) — подходит для любых m×n\n");
        printf("  2) Крамер — только если A квадратная (m==n) и det(A) != 0\n");

        int meth = 0;
        if (!read_int_line("Ваш выбор: ", &meth)) return 0;

        if (meth == 1 || meth == 2) return meth;
        printf("Некорректный выбор.\n");
    }
}

static void slau_print_result(const slau *S, OSHIBKA_NAFIK st)
{
    switch (st)
    {
        case SLAU_OK:
            printf("\nРешение единственное. x:\n");
            print_matrix(&S->x);
            break;

        case SLAU_INFINITE_SOL:
            printf("\nСистема совместна, но решений бесконечно много.\n");
            printf("Показано одно частное решение (свободные переменные = 0):\n");
            print_matrix(&S->x);
            break;

        case SLAU_INCOMPATIBLE:
            printf("\nСистема несовместна (решений нет).\n");
            break;

        case SLAU_SINGULAR:
            printf("\nМатрица вырождена (det(A)=0) — нет единственного решения.\n");
            break;

        case SLAU_METHOD_UNSUPPORTED:
            printf("\nВыбранный метод неприменим к данной системе.\n");
            break;

        case SLAU_DIM_MISMATCH:
            printf("\nОшибка размерностей: A должно быть m×n, b должно быть m×1.\n");
            break;

        default:
            printf("\nОшибка решения: code=%d\n", (int)st);
            break;
    }
}

// -------------------- MAIN -------------------- //

int main(void)
{
    srand((unsigned)time(NULL));

    matrix A, B, R;
    constructor(&A);
    constructor(&B);
    constructor(&R);

    // СЛАУ
    slau S;
    slau_constructor(&S);
    S.eps = 1e-12;

    input_and_fill(&A, "A");
    input_and_fill(&B, "B");

    for (;;)
    {
        printf("\n================= МЕНЮ =================\n");
        printf("---- Матрицы A/B ----\n");
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

        printf("\n---- СЛАУ A*x=b ----\n");
        printf("13) Напоминалка/инструкция по СЛАУ\n");
        printf("14) Ввести систему (A и b)\n");
        printf("15) Показать систему (A и b)\n");
        printf("16) Задать eps для СЛАУ (сейчас eps=%.3e)\n", S.eps);
        printf("17) Решить систему (выбор метода: Гаусс/Крамер)\n");
        printf("18) Показать последнее решение x\n");

        printf("\n0) Выход\n");

        int cmd = -1;
        if (!read_int_line("Ваш выбор: ", &cmd)) break;

        enum ERROS_TYPE e = OK;

        switch (cmd)
        {
            case 0:
                destructor(&A);
                destructor(&B);
                destructor(&R);
                slau_destructor(&S);
                return 0;

            // ---------- Матрицы ----------
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
                break;

            case 4:
                free_matrix(&B);
                input_and_fill(&B, "B");
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
                if (e == OK) { printf("\nA * %.6lf:\n", k); print_matrix(&R); }
                break;
            }

            case 9:
            {
                double det = 0.0;
                e = matrix_det_gauss(&A, &det);
                handle_op_error(e, "matrix_det_gauss(A)", &A, NULL);
                if (e == OK) printf("det(A) = %.12lf\n", det);
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

            // ---------- СЛАУ ----------
            case 13:
                instruction();
                reminder();
                break;

            case 14:
                slau_input_system(&S);
                break;

            case 15:
                if (validate_matrix_ro(&S.A) != OK || validate_matrix_ro(&S.b) != OK)
                {
                    printf("\nСистема ещё не введена или введена некорректно.\n");
                }
                else
                {
                    printf("\nA (%zux%zu):\n", S.A.rows, S.A.colums);
                    print_matrix(&S.A);
                    printf("\nb (%zux1):\n", S.b.rows);
                    print_matrix(&S.b);
                }
                break;

            case 16:
            {
                double eps = 0.0;
                if (!read_double_line("Введите eps (например 1e-12): ", &eps)) break;
                if (eps <= 0.0) { printf("eps должен быть > 0\n"); break; }
                S.eps = eps;
                printf("eps установлен: %.3e\n", S.eps);
                break;
            }

            case 17:
            {
                OSHIBKA_NAFIK vb = slau_validate_basic(&S);
                if (vb != SLAU_OK)
                {
                    printf("\nСначала введите систему (пункт 14).\n");
                    break;
                }

                int meth = slau_choose_method();
                if (meth == 0) break;

                if (meth == 1)
                {
                    S.choice = METHOD_GAUSS_JORDAN;
                }
                else
                {
                    S.choice = METHOD_KRAMER;
                }

                OSHIBKA_NAFIK st = slau_solve(&S);

                if (S.choice == METHOD_KRAMER && S.A.rows != S.A.colums)
                {
                    printf("\nКрамер неприменим: A не квадратная (%zux%zu).\n", S.A.rows, S.A.colums);
                    break;
                }

                slau_print_result(&S, st);
                break;
            }

            case 18:
                if (validate_matrix_ro(&S.x) != OK)
                    printf("\nРешение ещё не получено.\n");
                else
                {
                    printf("\nx (%zux1):\n", S.x.rows);
                    print_matrix(&S.x);
                }
                break;

            default:
                fprintf(stderr, "Неизвестная команда.\n");
                break;
        }
    }

    destructor(&A);
    destructor(&B);
    destructor(&R);
    slau_destructor(&S);
    return 0;
}