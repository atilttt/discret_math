#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// -------------------------------------------------
// === ВЫБОР МЕТОДА (раскомментируйте ОДИН) ===
//#define METHOD_TAYLOR
#define METHOD_CHEBYSHEV
// -------------------------------------------------

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// -------------------------------
// Общая функция: точное значение
// -------------------------------
double exact_arctanh(double x) {
    return 0.5 * log((1.0 + x) / (1.0 - x));
}

// =================================================
// МЕТОД ТЕЙЛОРА
// =================================================
#ifdef METHOD_TAYLOR

int main(void) {
    double x;
    const double eps = 1e-2;

    printf("=== Ряд Тейлора для f(x) = 0.5 * ln((1+x)/(1-x)) ===\n");
    printf("Введите x в интервале (-1, 1): ");
    if (scanf("%lf", &x) != 1) {
        printf("Ошибка ввода.\n");
        return 1;
    }

    if (fabs(x) >= 1.0) {
        printf("Ошибка: x должно быть строго внутри (-1, 1).\n");
        return 1;
    }

    double sum = 0.0;
    double term = x; // первый член: x^(2*0+1)/(2*0+1) = x
    int terms_used = 0;
    int n = 0; // индекс члена: n = 0, 1, 2, ...

    while (fabs(term) > eps) {
        sum += term;
        terms_used++;
        n++;
        // Следующий член: x^(2n+1)/(2n+1)
        // Рекуррентно: term_{n} = term_{n-1} * x^2 * (2n-1)/(2n+1)
        term = term * x * x * (2 * n - 1) / (2 * n + 1);
    }

    double exact = exact_arctanh(x);
    double error = fabs(sum - exact);

    printf("\n--- Результат (ряд Тейлора) ---\n");
    printf("x = %.6f\n", x);
    printf("Приближение: %.6f\n", sum);
    printf("Точное значение: %.6f\n", exact);
    printf("Погрешность: %.6f\n", error);
    printf("Количество членов: %d\n", terms_used);

    return 0;
}

// =================================================
// МЕТОД ЧЕБЫШЕВА
// =================================================
#elif defined(METHOD_CHEBYSHEV)

// Рекуррентное вычисление T_k(x)
double chebyshev(int k, double x) {
    if (k == 0) return 1.0;
    if (k == 1) return x;
    double t0 = 1.0, t1 = x, t2;
    for (int i = 2; i <= k; i++) {
        t2 = 2.0 * x * t1 - t0;
        t0 = t1;
        t1 = t2;
    }
    return t1;
}

// Подынтегральная функция для коэффициента c_k
double integrand(int k, double x) {
    if (fabs(x) >= 1.0) return 0.0;
    double f_x = exact_arctanh(x);
    double weight = 1.0 / sqrt(1.0 - x * x);
    return f_x * chebyshev(k, x) * weight;
}

// Численное интегрирование (трапеции)
double trapezoidal(int k, int n_points) {
    // Избегаем концов отрезка из-за особенности
    const double delta = 1e-12;
    double a = -1.0 + delta;
    double b =  1.0 - delta;
    double h = (b - a) / (n_points - 1);
    double sum = 0.0;

    for (int i = 0; i < n_points; i++) {
        double xi = a + i * h;
        double fx = integrand(k, xi);
        if (i == 0 || i == n_points - 1)
            sum += fx * 0.5;
        else
            sum += fx;
    }
    return sum * h;
}

void compute_coeffs(double* coeffs, int N, int n_points) {
    for (int k = 0; k < N; k++) {
        double I = trapezoidal(k, n_points);
        coeffs[k] = (k == 0) ? I / M_PI : 2.0 * I / M_PI;
    }
}

double evaluate_cheb(double x, double* coeffs, int N) {
    double s = 0.0;
    for (int k = 0; k < N; k++) {
        s += coeffs[k] * chebyshev(k, x);
    }
    return s;
}

int main(void) {
    double x;
    const double eps = 1e-2;
    const int n_points = 20000; // для интегрирования
    const int MAX_N = 15;

    printf("=== Разложение по полиномам Чебышева ===\n");
    printf("Введите x в интервале (-1, 1): ");
    if (scanf("%lf", &x) != 1) {
        printf("Ошибка ввода.\n");
        return 1;
    }

    if (fabs(x) >= 1.0) {
        printf("Ошибка: x должно быть строго внутри (-1, 1).\n");
        return 1;
    }

    double exact = exact_arctanh(x);
    int N_used = 0;
    double approx = 0.0;
    double* coeffs = malloc(MAX_N * sizeof(double));

    if (!coeffs) {
        printf("Ошибка выделения памяти.\n");
        return 1;
    }

    printf("Подбор минимального числа членов для точности %.2e...\n", eps);

    for (int N = 1; N <= MAX_N; N++) {
        compute_coeffs(coeffs, N, n_points);
        approx = evaluate_cheb(x, coeffs, N);
        double err = fabs(approx - exact);

        printf("  N = %2d → приближение = %8.6f, погрешность = %8.6f\n", N, approx, err);

        if (err <= eps) {
            N_used = N;
            break;
        }
    }

    printf("\n--- Результат (Чебышев) ---\n");
    if (N_used > 0) {
        printf("x = %.6f\n", x);
        printf("Приближение: %.6f\n", approx);
        printf("Точное значение: %.6f\n", exact);
        printf("Погрешность: %.6f\n", fabs(approx - exact));
        printf("Количество членов: %d\n", N_used);
    } else {
        printf("Не удалось достичь точности %.2e даже с %d членами.\n", eps, MAX_N);
    }

    free(coeffs);
    return 0;
}

#endif