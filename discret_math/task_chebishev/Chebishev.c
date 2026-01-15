#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MAX_DEGREE 60

// ============================================
// 1. СПОСОБЫ ВЫЧИСЛЕНИЯ ПОЛИНОМОВ ЧЕБЫШЕВА
// ============================================

double chebyshev_horner(int n, double x) {
    if (n == 0) return 1.0;
    if (n == 1) return x;
    
    double Tn_2 = 1.0; 
    double Tn_1 = x;   
    double Tn = 0.0;
    
    for (int i = 2; i <= n; i++) {
        Tn = 2.0 * x * Tn_1 - Tn_2;
        Tn_2 = Tn_1;
        Tn_1 = Tn;
    }
    
    return Tn;
}

double chebyshev_trig(int n, double x) {
    if (x < -1.0 || x > 1.0) {
        return cosh(n * acosh(x));
    }
    return cos(n * acos(x));
}

// ============================================
// 2. КОЭФФИЦИЕНТЫ ПОЛИНОМОВ ЧЕБЫШЕВА
// ============================================

void chebyshev_coefficients(int n, double coeffs[]) {
    for (int i = 0; i <= n; i++) {
        coeffs[i] = 0.0;
    }
    
    if (n == 0) {
        coeffs[0] = 1.0;
        return;
    }
    
    if (n == 1) {
        coeffs[0] = 0.0;
        coeffs[1] = 1.0;
        return;
    }
    
    double *T0 = (double*)malloc(sizeof(double) * 1);
    double *T1 = (double*)malloc(sizeof(double) * 2);
    double *Tn = NULL;
    
    T0[0] = 1.0;  
    T1[0] = 0.0;  
    T1[1] = 1.0;
    
    for (int i = 2; i <= n; i++) {
        Tn = (double*)malloc(sizeof(double) * (i + 1));
        for (int j = 0; j <= i; j++) {
            Tn[j] = 0.0;
        }
    
        for (int j = 0; j <= i-1; j++) {
            if (j + 1 <= i) {
                Tn[j + 1] += 2.0 * T1[j];
            }
        }
        
        for (int j = 0; j <= i-2; j++) {
            Tn[j] -= T0[j];
        }
        
        free(T0);
        T0 = T1;
        T1 = Tn;
    }
    
    for (int i = 0; i <= n; i++) {
        coeffs[i] = T1[i];
    }
    
    free(T0);
    free(T1);
}

double evaluate_polynomial(double coeffs[], int degree, double x) {
    double result = coeffs[degree];
    for (int i = degree - 1; i >= 0; i--) {
        result = result * x + coeffs[i];
    }
    return result;
}

void print_polynomial(double coeffs[], int degree) {
    int first = 1;
    
    for (int i = degree; i >= 0; i--) {
        if (fabs(coeffs[i]) < 1e-10) continue;
        
        if (!first) {
            if (coeffs[i] >= 0) {
                printf(" + ");
            } else {
                printf(" - ");
            }
        } else {
            first = 0;
            if (coeffs[i] < 0) {
                printf("-");
            }
        }
        
        double abs_coeff = fabs(coeffs[i]);
        
        if (i == 0) {
            printf("%.3f", abs_coeff);
        } else if (i == 1) {
            if (fabs(abs_coeff - 1.0) < 1e-10) {
                printf("x");
            } else {
                printf("%.3f*x", abs_coeff);
            }
        } else {
            if (fabs(abs_coeff - 1.0) < 1e-10) {
                printf("x^%d", i);
            } else {
                printf("%.3f*x^%d", abs_coeff, i);
            }
        }
    }
    
    if (first) {
        printf("0");
    }
}

// ============================================
// ПУНКТ 1: ПЕРВЫЕ 8 ПОЛИНОМОВ ЧЕБЫШЕВА
// ============================================
#ifdef PART1

int main() {
    printf("===========================================\n");
    printf("ПУНКТ 1: ПЕРВЫЕ 8 ПОЛИНОМОВ ЧЕБЫШЕВА\n");
    printf("===========================================\n");
    
    double test_points[] = {-0.9, -0.5, 0.0, 0.3, 0.7, 1.0};
    int num_points = sizeof(test_points) / sizeof(test_points[0]);
    
    for (int n = 0; n < 8; n++) {
        printf("\nПолином Чебышева T_%d(x):\n", n);
        
        double coeffs[MAX_DEGREE + 1];
        chebyshev_coefficients(n, coeffs);
        
        printf("Формула: T_%d(x) = ", n);
        print_polynomial(coeffs, n);
        printf("\n");
        
        printf("Значения в точках:\n");
        printf("%-10s %-20s %-20s %-10s\n", 
               "x", "Горнер", "Тригонометрия", "Разница");
        printf("------------------------------------------------------------\n");
        
        for (int i = 0; i < num_points; i++) {
            double x = test_points[i];
            double val_horner = chebyshev_horner(n, x);
            double val_trig = chebyshev_trig(n, x);
            double diff = fabs(val_horner - val_trig);
            
            printf("%-10.2f %-20.15f %-20.15f %-10.2e\n", 
                   x, val_horner, val_trig, diff);
        }
    }
    
    printf("\n\n===========================================\n");
    printf("ВЫЧИСЛЕНИЕ В ПРОИЗВОЛЬНЫХ ТОЧКАХ\n");
    printf("===========================================\n");
    
    double x;
    char continue_calc = 'y';
    
    while (continue_calc == 'y' || continue_calc == 'Y') {
        printf("\nВведите степень полинома (0-7): ");
        int n;
        scanf("%d", &n);
        
        if (n < 0 || n > 7) {
            printf("Ошибка: степень должна быть от 0 до 7\n");
            continue;
        }
        
        printf("Введите значение x: ");
        scanf("%lf", &x);
        
        double val_horner = chebyshev_horner(n, x);
        double val_trig = chebyshev_trig(n, x);
        
        printf("\nT_%d(%.3f) = %.15f (метод Горнера)\n", n, x, val_horner);
        printf("T_%d(%.3f) = %.15f (тригонометрический метод)\n", n, x, val_trig);
        printf("Разница: %e\n", fabs(val_horner - val_trig));
        
        printf("\nХотите вычислить еще одно значение? (y/n): ");
        scanf(" %c", &continue_calc);
    }
    
    return 0;
}

#endif

// ============================================
// ПУНКТ 2: ПОЛИНОМ 60-Й СТЕПЕНИ
// ============================================
#ifdef PART2

int main() {
    printf("===========================================\n");
    printf("ПУНКТ 2: ПОЛИНОМ ЧЕБЫШЕВА 60-Й СТЕПЕНИ\n");
    printf("===========================================\n");
    
    const int n = 60;
    double coeffs[MAX_DEGREE + 1];
    
    printf("\nВычисление коэффициентов полинома T_60(x)...\n");
    chebyshev_coefficients(n, coeffs);
    
    printf("\nКоэффициенты полинома T_60(x):\n");
    printf("(формат: степень: коэффициент)\n");
    printf("------------------------------\n");
    
    for (int i = 0; i <= n; i++) {
        if (fabs(coeffs[i]) > 1e-10) {
            printf("x^%2d: % .10e\n", i, coeffs[i]);
        }
    }
    
    printf("\nПолином в сокращенном виде:\n");
    printf("T_60(x) = ");
    int terms_shown = 0;
    
    for (int i = n; i >= 0; i--) {
        if (fabs(coeffs[i]) > 1e-10) {
            if (terms_shown < 3 || i < 3) {
                if (terms_shown > 0) {
                    if (coeffs[i] >= 0) printf(" + ");
                    else printf(" - ");
                }
                
                double abs_coeff = fabs(coeffs[i]);
                if (i == 0) printf("%.3e", abs_coeff);
                else if (i == 1) printf("%.3e*x", abs_coeff);
                else printf("%.3e*x^%d", abs_coeff, i);
                
                terms_shown++;
            } else if (terms_shown == 3) {
                printf(" + ...");
                terms_shown++;
            }
        }
    }
    printf("\n");
    
    printf("\n===========================================\n");
    printf("ВЫЧИСЛЕНИЕ ЗНАЧЕНИЙ ПОЛИНОМА T_60(x)\n");
    printf("===========================================\n");
    
    double x;
    char continue_calc = 'y';
    
    while (continue_calc == 'y' || continue_calc == 'Y') {
        printf("\nВведите значение x: ");
        scanf("%lf", &x);
        
        double val_horner = chebyshev_horner(n, x);
        double val_trig = chebyshev_trig(n, x);
        double val_coeffs = evaluate_polynomial(coeffs, n, x);
        
        printf("\nT_60(%.3f) = %.15e (метод Горнера)\n", x, val_horner);
        printf("T_60(%.3f) = %.15e (тригонометрический метод)\n", x, val_trig);
        printf("T_60(%.3f) = %.15e (по коэффициентам)\n", x, val_coeffs);
        printf("\nРазницы:\n");
        printf("Горнер vs Тригонометрия: %e\n", fabs(val_horner - val_trig));
        printf("Горнер vs Коэффициенты:  %e\n", fabs(val_horner - val_coeffs));
        printf("Тригонометрия vs Коэффициенты: %e\n", fabs(val_trig - val_coeffs));
        
        printf("\nХотите вычислить еще одно значение? (y/n): ");
        scanf(" %c", &continue_calc);
    }
    
    printf("\n===========================================\n");
    printf("ТЕСТОВЫЕ ЗНАЧЕНИЯ\n");
    printf("===========================================\n");
    
    double test_points[] = {-0.9, -0.5, 0.0, 0.3, 0.7, 1.0, 1.5, 2.0};
    int num_points = sizeof(test_points) / sizeof(test_points[0]);
    
    printf("\n%-10s %-25s %-25s %-25s\n", 
           "x", "Горнер", "Тригонометрия", "Коэффициенты");
    printf("--------------------------------------------------------------------------------\n");
    
    for (int i = 0; i < num_points; i++) {
        x = test_points[i];
        double val_horner = chebyshev_horner(n, x);
        double val_trig = chebyshev_trig(n, x);
        double val_coeffs = evaluate_polynomial(coeffs, n, x);
        
        printf("%-10.2f %-25.15e %-25.15e %-25.15e\n", 
               x, val_horner, val_trig, val_coeffs);
    }
    
    return 0;
}

#endif

// ============================================
// ПУНКТ 4: НУЛИ И ЭКСТРЕМУМЫ T7(x)
// ============================================
#ifdef PART4

int main() {
    printf("===========================================\n");
    printf("ПУНКТ 4: НУЛИ И ЭКСТРЕМУМЫ T7(x)\n");
    printf("===========================================\n\n");
    
    const int n = 7;
    
    printf("1. НУЛИ T7(x) (аналитическое решение):\n");
    printf("   Формула: x_k = cos(pi*(2k-1)/14), k = 1..7\n");
    printf("   k   |   Формула    |     x_k     |  T7(x_k) (вычисл.)\n");
    printf("  -----|--------------|-------------|-------------------\n");
    
    double zeros[7];
    for (int k = 1; k <= n; k++) {
        double theta = M_PI * (2.0 * k - 1.0) / 14.0;
        zeros[k-1] = cos(theta);
        double T_value = chebyshev_horner(n, zeros[k-1]);
        printf("   %d   | cos(%dπ/14) | %11.8f | %+.6e\n", 
               k, (2*k-1), zeros[k-1], T_value);
    }
    
    printf("\n\n2. ЭКСТРЕМУМЫ T7(x) (аналитическое решение):\n");
    printf("   Формула: x_k = cos(pi*k/7), k = 0..7\n");
    printf("   T7(x_k) = (-1)^k\n");
    printf("   k   |   Формула    |     x_k     | T7(x_k) |   Тип\n");
    printf("  -----|--------------|-------------|---------|--------\n");
    
    double extrema[8];
    for (int k = 0; k <= n; k++) {
        double theta = M_PI * k / n;
        extrema[k] = cos(theta);
        double expected = (k % 2 == 0) ? 1.0 : -1.0;
        double T_value = chebyshev_horner(n, extrema[k]);
        const char* type;
        if (k == 0) type = "Максимум";
        else if (k == n) type = "Минимум";
        else if (k % 2 == 0) type = "Максимум";
        else type = "Минимум";
        
        if (k == 0) {
            printf("   %d   |   cos(0)    | %11.8f |   %+4.1f  | %s\n", 
                   k, extrema[k], expected, type);
        } else if (k == n) {
            printf("   %d   |   cos(π)    | %11.8f |   %+4.1f  | %s\n", 
                   k, extrema[k], expected, type);
        } else {
            printf("   %d   | cos(%dπ/7)  | %11.8f |   %+4.1f  | %s\n", 
                   k, k, extrema[k], expected, type);
        }
    }
    
    
    return 0;
}

#endif

// ============================================
// ЕСЛИ НИ ОДИН ИЗ ПУНКТОВ НЕ ВЫБРАН
// ============================================
#if !defined(PART1) && !defined(PART2) && !defined(PART4)

int main() {
    printf("Ошибка: Не выбран режим компиляции.\n");
    printf("Используйте одну из опций:\n");
    printf("  -DPART1  для первого пункта\n");
    printf("  -DPART2  для второго пункта\n");
    printf("  -DPART4  для четвертого пункта\n");
    return 1;
}

#endif