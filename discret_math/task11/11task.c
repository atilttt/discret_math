#include <stdio.h>
#include <math.h>

// Твоя функция: y = x^2 * cos(x)
double f(double x) {
    return x * x * cos(x);
}

// Точная производная: y' = 2x*cos(x) - x^2*sin(x)
double f_exact_derivative(double x) {
    return 2*x*cos(x) - x*x*sin(x);
}

// Три метода численного дифференцирования
double forward_difference(double x, double h) {
    return (f(x + h) - f(x)) / h;
}

double central_difference(double x, double h) {
    return (f(x + h) - f(x - h)) / (2*h);
}

double fourth_order_difference(double x, double h) {
    return (-f(x + 2*h) + 8*f(x + h) - 8*f(x - h) + f(x - 2*h)) / (12*h);
}

// Функция для вычисления производных и ошибок
void compute_derivatives(double x, double h, double results[7]) {
    // results[0] = точное значение
    // results[1] = правая разность
    // results[2] = центральная разность
    // results[3] = 4-го порядка
    // results[4] = ошибка правой
    // results[5] = ошибка центральной
    // results[6] = ошибка 4-го порядка
    
    results[0] = f_exact_derivative(x);
    results[1] = forward_difference(x, h);
    results[2] = central_difference(x, h);
    results[3] = fourth_order_difference(x, h);
    
    results[4] = fabs(results[1] - results[0]);
    results[5] = fabs(results[2] - results[0]);
    results[6] = fabs(results[3] - results[0]);
}

int main() {
    double x_point = 1.5;  // Твоя точка
    double exact = f_exact_derivative(x_point);
    
    printf("============================================================\n");
    printf("ЧИСЛЕННОЕ ДИФФЕРЕНЦИРОВАНИЕ\n");
    printf("Функция: y = x^2 * cos(x)\n");
    printf("Точка x = %.2f\n", x_point);
    printf("Точная производная: y' = 2x*cos(x) - x^2*sin(x) = %.10f\n", exact);
    printf("============================================================\n\n");
    
    // Часть 1: Вычисление при разных фиксированных h
    printf("1. ЗНАЧЕНИЯ ПРОИЗВОДНОЙ ПРИ РАЗНЫХ h:\n");
    printf("================================================================================\n");
    printf("|   h    |   Точное   |   Правая   |  Центральная |  4-й порядок |\n");
    printf("================================================================================\n");
    
    double h_values[] = {0.5, 0.1, 0.05, 0.01, 0.001, 0.0001};
    int n = 6;
    
    for(int i = 0; i < n; i++) {
        double h = h_values[i];
        double results[7];
        compute_derivatives(x_point, h, results);
        
        printf("| %6.4f | %10.6f | %10.6f | %12.6f | %12.6f |\n",
               h, results[0], results[1], results[2], results[3]);
    }
    printf("================================================================================\n\n");
    
    // Часть 2: Ошибки при разных h
    printf("2. ОШИБКИ (разница с точным значением):\n");
    printf("=================================================================\n");
    printf("|   h    |  Ошибка правой  | Ошибка центральной | Ошибка 4-го пор.|\n");
    printf("=================================================================\n");
    
    for(int i = 0; i < n; i++) {
        double h = h_values[i];
        double results[7];
        compute_derivatives(x_point, h, results);
        
        printf("| %6.4f |   %12.6e   |   %12.6e    |   %12.6e   |\n",
               h, results[4], results[5], results[6]);
    }
    printf("=================================================================\n\n");
    
    // Часть 3: Поиск оптимальных h для каждого метода
    printf("3. ПОИСК ОПТИМАЛЬНОГО h ДЛЯ КАЖДОГО МЕТОДА:\n");
    printf("(перебираем h от 1e-1 до 1e-15)\n\n");
    
    // Для правой разности
    printf("а) Правая разность (несимметричная схема):\n");
    double best_h_forward = 0.1;
    double best_err_forward = 1e10;
    
    for(double h = 1e-1; h > 1e-15; h *= 0.1) {
        double approx = forward_difference(x_point, h);
        double err = fabs(approx - exact);
        printf("   h = %7.1e → значение = %12.8f → ошибка = %7.1e\n", 
               h, approx, err);
        
        if(err < best_err_forward) {
            best_err_forward = err;
            best_h_forward = h;
        }
        
        if(h < 1e-10) break;  // Не перебирать слишком долго
    }
    printf("   ЛУЧШИЙ: h = %7.1e, ошибка = %7.1e\n\n", best_h_forward, best_err_forward);
    
    // Для центральной разности
    printf("б) Центральная разность (симметричная схема):\n");
    double best_h_central = 0.1;
    double best_err_central = 1e10;
    
    for(double h = 1e-1; h > 1e-15; h *= 0.1) {
        double approx = central_difference(x_point, h);
        double err = fabs(approx - exact);
        printf("   h = %7.1e → значение = %12.8f → ошибка = %7.1e\n", 
               h, approx, err);
        
        if(err < best_err_central) {
            best_err_central = err;
            best_h_central = h;
        }
        
        if(h < 1e-12) break;
    }
    printf("   ЛУЧШИЙ: h = %7.1e, ошибка = %7.1e\n\n", best_h_central, best_err_central);
    
    // Для схемы 4-го порядка
    printf("в) Схема 4-го порядка:\n");
    double best_h_fourth = 0.1;
    double best_err_fourth = 1e10;
    
    for(double h = 1e-1; h > 1e-10; h *= 0.1) {
        double approx = fourth_order_difference(x_point, h);
        double err = fabs(approx - exact);
        printf("   h = %7.1e → значение = %12.8f → ошибка = %7.1e\n", 
               h, approx, err);
        
        if(err < best_err_fourth) {
            best_err_fourth = err;
            best_h_fourth = h;
        }
    }
    printf("   ЛУЧШИЙ: h = %7.1e, ошибка = %7.1e\n\n", best_h_fourth, best_err_fourth);
    
    // Часть 4: Итоговая таблица
    printf("4. ИТОГОВАЯ ТАБЛИЦА (ОПТИМАЛЬНЫЕ ЗНАЧЕНИЯ):\n");
    printf("===========================================================\n");
    printf("|       Метод        |  Оптимальный h  | Минимальная ошибка |\n");
    printf("===========================================================\n");
    printf("| Правая разность    |    %7.1e     |      %7.1e      |\n", 
           best_h_forward, best_err_forward);
    printf("| Центральная        |    %7.1e     |      %7.1e      |\n", 
           best_h_central, best_err_central);
    printf("| 4-го порядка       |    %7.1e     |      %7.1e      |\n", 
           best_h_fourth, best_err_fourth);
    printf("===========================================================\n\n");
    
    // Часть 5: Выводы
    printf("5. ВЫВОДЫ:\n");
    printf("   1. Точное значение производной в x=1.5: %.10f\n", exact);
    printf("   2. Самый точный метод: схема 4-го порядка\n");
    printf("   3. Оптимальные шаги h разные для каждого метода:\n");
    printf("      - Правая разность: h ≈ 1e-8\n");
    printf("      - Центральная: h ≈ 1e-6\n");
    printf("      - 4-й порядок: h ≈ 1e-4\n");
    printf("   4. При слишком маленьких h ошибка растет из-за округления\n");
    
    return 0;
}