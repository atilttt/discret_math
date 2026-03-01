#ifndef SLAU_H
#define SLAU_H

#include "matrix_lib.h"
#include <stddef.h>

typedef enum {
    SLAU_OK = 0,
    SLAU_BAD_INPUT,
    SLAU_DIM_MISMATCH,
    SLAU_INCOMPATIBLE,     // несовместна
    SLAU_INFINITE_SOL,     // бесконечно много решений (совместна, но не единственна)
    SLAU_SINGULAR,         // det(A)=0 для квадратной
    SLAU_METHOD_UNSUPPORTED
} OSHIBKA_NAFIK;

typedef enum {
    METHOD_NONE = 0,
    METHOD_GAUSS_JORDAN = 1,
    METHOD_KRAMER = 2
} SLAU_METHOD;

typedef struct
{
    SLAU_METHOD choice;

    matrix A;   // m x n
    matrix b;   // m x 1
    matrix x;   // n x 1 (результат)

    double eps; // порог
    OSHIBKA_NAFIK err;
} slau;

void reminder(void);
void instruction(void);

void slau_constructor(slau *s);
void slau_destructor(slau *s);

OSHIBKA_NAFIK slau_validate_basic(const slau *s);
OSHIBKA_NAFIK slau_check_kronecker_capelli(const slau *s);

OSHIBKA_NAFIK slau_solve(slau *s);

OSHIBKA_NAFIK slau_solve_gauss_jordan(slau *s);
OSHIBKA_NAFIK slau_solve_kramer(slau *s);

#endif // SLAU_H