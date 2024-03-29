#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include<stdio.h>

using namespace std;

ofstream Out1("Mat.txt", ios::out); // В этот файл записываем приближенные решения с условием Дирихле всюду
ofstream Out2("Out.txt", ios::out); //В этот файл записываем приближенные решения с условием второго рода на левой границе и условием третьего рода на правой
ofstream Out3("Out3.txt", ios::out); //В этот файл записываем данные точного решения
ofstream Out4("Out4.txt", ios::out); //В этот файл записываем максимальные отклонения приближенного решения от точного на слое, вычисляем порядки для задачи Дирихле
ofstream Out5("Out5.txt", ios::out); //В этот файл записываем максимальные отклонения приближенного решения от точного на слое, вычисляем порядки для задачи с заданными краевыми условиями

ofstream Out6("Out6.txt", ios::out); // В этот файл записываем приближенные решения с условием Дирихле всюду повышенный порядок
ofstream Out7("Out7.txt", ios::out); //В этот файл записываем максимальные отклонения приближенного решения задачи Дирихле, вычисляем порядки для задачи с заданными краевыми условиями повышеныый порядок
ofstream Out8("Out8.txt", ios::out); //В этот файл записываем приближенные решения с условием второго рода на левой границе и условием третьего рода на правой повышенный порядок
ofstream Out9("Out9.txt", ios::out); //В этот файл записываем максимальные отклонения приближенного решения с заданными условиями, вычисляем порядки для задачи с заданными краевыми условиями повышенный порядок

const double a = 0, b = 1; // Область изменения х
const double T = 1; // Время
const double beta = 1; // Коэффициент теплоотдачи

double exactsln(double t, double x) // Точное решение
{
    double f;
    f = sin(3 * x + 2 * t);
    return f;
}

double f(double t, double x) // Правая часть уравнения теплопроводности
{
    double f;
    f = 2 * cos(3 * x + 2 * t) + 9 * sin(3 * x + 2 * t);
    return f;
}

double fTaylor(double t, double x, double tau, double h) // Правая часть для схемы повышенной точности
{
    double r;
    r = (f(t + tau, x) + f(t, x)) / 2 + (f(t, x + h) - 2 * f(t, x) + f(t, x - h)) / (12);
    return r;
}

double initial(double x) // Начальное данное
{
    double f;
    f = sin(3 * x);
    return f;
}

double border01(double t) // Краевое условие первого рода при х = а
{
    double f;
    f = sin(3 * a + 2 * t);
    return f;
}

double border02(double t) // Краевое условие первого рода при х = b
{
    double f;
    f = sin(3 * b + 2 * t);
    return f;
}

double border11(double t) // Краевое условие второго рода при х = а
{
    double f;
    f = 3*cos(3 * a + 2 * t);
    return f;
}

double border22(double t) // Краевое условие второго рода при х = b
{
    double f;
    f = beta*( sin(3*b+2*t) + (3*cos(3*b+2*t))/beta);
    return f;
}

double* grid(double* u_n, int n, double tau, double h, int N, double sigma, int parametr) // n - временной слой, который известен, функция находит правую часть СЛАУ для реализации метода прогонки
{
    double* w = new double[N + 2];
    if (parametr == 0) //Для схемы обычной точности
    {
        for (int i = 1; i <= N - 1; i++)
        {
            w[i] = tau * (1 - sigma) * u_n[i + 1] + (h * h - 2 * tau + 2 * tau * sigma) * u_n[i] + tau * (1 - sigma) * u_n[i - 1] + tau * h * h * f(n * tau, a + i * h);

        }
    }
    else //Для схемы повышенной точности
    {
        for (int i = 1; i <= N - 1; i++)
        {
            w[i] = tau * (1 - sigma) * u_n[i + 1] + (h * h - 2 * tau + 2 * tau * sigma) * u_n[i] + tau * (1 - sigma) * u_n[i - 1] + tau * h * h * fTaylor(n * tau, a + i * h, tau, h);

        }
    }
    
    return w;
    
}

double NormOfDifference(double* u_n, int N) //Вводим равномерную норму на слое для нахождения максимального отклонения приближенного решения от точного на слое
{
    double max = 0;
    double h1 = (b - a) / N;
    for (int i = 0; i <= N; i++)
    {
        if (abs(exactsln(T, a + i * h1) - u_n[i]) > max) {
            max = abs(exactsln(T, a + i * h1) - u_n[i]);
        }
    }
    return max;
}
void progonka(int q, double* u_n, double* u_n1, double tau, double h, int N, double sigma, int parametr) // Метод прогонки для условия Дирихле (q - временной слой, который находим, j - индекс по пространственной переменной) 
{
	double a, b, c, b0, c0, aN, bN, cN;
	double* xi = new double[N + 2];
	double* eta = new double[N + 2];
	a = c = -(tau*sigma);
    b = 2 * tau * sigma + h * h;
	b0 = 1;
	c0 = 0;
	bN = 1;
	aN = 0;
	cN = 0;
	xi[0] = 0;
	eta[0] = 0;
	xi[1] = -c0 / b0;
	eta[1] = border01(q*tau)/b0;
   
    for (int i = 1; i <= N; i++)
    {
        xi[i + 1] = -c / (b + xi[i] * a);
        eta[i + 1] = (grid(u_n, q - 1, tau, h, N, sigma, parametr)[i] - eta[i] * a) / (b + xi[i] * a);
    }
    xi[N + 1] = 0;
    eta[N + 1] = (border02(q * tau) - aN * eta[N]) / (aN * xi[N] + bN);

    u_n1[N] = eta[N + 1];
	for (int i = N - 1; i >= 0; i--)
	{
        u_n1[i] = xi[i + 1] * u_n1[i+1] + eta[i + 1];
	}
	
	delete[] xi;
	delete[] eta;
}
void progonkaKr(int q, double* u_n, double* u_n1, double tau, double h, int N, double sigma, int parametr) //Метод прогонки для задачи с заданными краевыми условиями (q - временной слой, который находим, j - индекс по пространственной переменной) 
{
    double a, b, c, b0, c0, aN, bN, cN;
    double* xi = new double[N + 2];
    double* eta = new double[N + 2];
    a = c = -(tau * sigma);
    b = 2 * tau * sigma + h * h;
    cN = 0;
    xi[0] = 0;
    eta[0] = 0;
    if (parametr == 0) //Схема обычной точности
    {
        b0 = -1 / h;
        c0 = 2 / h - (2 * tau * sigma + h * h) / (tau * sigma * 2 * h);
        bN = 1 / h + beta;
        aN = -2 / h + (2 * tau * sigma + h * h) / (tau * sigma * 2 * h);
        eta[1] = (border11(q * tau) - (grid(u_n, q - 1, tau, h, N, sigma, parametr)[1]) / (tau * sigma * 2 * h)) / b0;
    }
    else //Схема повышенной точности
    {
        b0 = -25 / (12 * h) + (13 * tau * tau * sigma * sigma - 4 * h * h * tau * sigma + 3 * h * h * h * h) / (12 * tau * tau * sigma * sigma * h);
        c0 = (38 * tau * sigma + 3 * h * h) / (12 * tau * sigma * h) + (2 * tau * sigma + h * h) * (-13 * tau * tau * sigma * sigma + 4 * h * h * tau * sigma - 3 * h * h * h * h) / (12 * tau * tau * tau * sigma * sigma * sigma * h);
        aN = (-38 * tau * sigma - 3 * h * h) / (12 * tau * sigma * h)+ (2 * tau * sigma + h * h)* (13 * tau * tau * sigma * sigma - 4 * h * h * tau * sigma + 3 * h * h * h * h) / (12 * tau * tau * tau * sigma * sigma * sigma * h);
        bN = 25/(12*h) + beta + (-13 * tau * tau * sigma * sigma + 4 * h * h * tau * sigma - 3 * h * h * h * h) / (12 * tau * tau * sigma * sigma * h);
        double r = ( - 1 / (4 * h * tau * sigma)) * grid(u_n, q - 1, tau, h, N, sigma, parametr)[3] + ((10 * tau * sigma - 3 * h * h) / (12 * h * tau * tau * sigma * sigma)) * grid(u_n, q - 1, tau, h, N, sigma, parametr)[2] + ((-13 * tau * tau * sigma * sigma + 4 * h * h * tau * sigma - 3 * h * h * h * h) / (12 * tau * tau * tau * sigma * sigma * sigma * h)) * grid(u_n, q - 1, tau, h, N, sigma, parametr)[1];
        eta[1] = (border11(q * tau) + r) / b0;
    }
    xi[1] = -c0 / b0;
    for (int i = 1; i <= N; i++)
    {
        xi[i + 1] = -c / (b + xi[i] * a);
        eta[i + 1] = (grid(u_n, q - 1, tau, h, N, sigma, parametr)[i] - eta[i] * a) / (b + xi[i] * a);
    }
    xi[N + 1] = 0;

    if (parametr == 0) //Схема обычной точности
    {
        eta[N + 1] = (border22(q * tau) + (grid(u_n, q - 1, tau, h, N, sigma, parametr)[N - 1]) / (tau * sigma * 2 * h) - aN * eta[N]) / (aN * xi[N] + bN);
    }
    else //Схема повышенной точности
    {
        double r = (1 / (4 * h * tau * sigma)) * grid(u_n, q - 1, tau, h, N, sigma, parametr)[N-3] + ((-10 * tau * sigma + 3 * h * h) / (12 * h * tau * tau * sigma * sigma)) * grid(u_n, q - 1, tau, h, N, sigma, parametr)[N-2] + ((13 * tau * tau * sigma * sigma - 4 * h * h * tau * sigma + 3 * h * h * h * h) / (12 * tau * tau * tau * sigma * sigma * sigma * h)) * grid(u_n, q - 1, tau, h, N, sigma, parametr)[N-1];
        eta[N + 1] = (border22(q * tau) + r - aN * eta[N]) / (aN * xi[N] + bN);
    }

    u_n1[N] = eta[N + 1];
    for (int i = N - 1; i >= 0; i--)
    {
        u_n1[i] = xi[i + 1] * u_n1[i + 1] + eta[i + 1];
    }

    delete[] xi;
    delete[] eta;
}



int main()
{
    int N = 5; // Кол-во шагов разбиения по х
    int M = 20; // Кол-во шагов разбиения по t
    double tau = T / M;
    double l = b - a;
    double h = l / N;
    for (double j = a; j <= b; j += 0.01) //Вычисляем точное решение
    {
        Out3 << j << " " << exactsln(M * tau, j) << endl;
    }

    int Parametr1 = 0; //Выбор между обычной схемой и схемой повышенной точности (0 - обычная схема, 1 - схема повышенной точности)
    int Parametr2 = 0; //Выбор между условием Дирихле и условием второго рода на левой границе условием третьего рода на правой границе (0 - всюду условие Дирихле, 1 - заданные условия)

    if (Parametr1 == 0) {
        double sigma = 0.7; // Вес схемы
        double max1 = 0; //Для отыскания максимального отклонения на слое
        double max2 = 0;

        if (Parametr2 == 0) //Условие Дирихле всюду
        {

            for (int set = 1; set <= 5; set++) // Измельчение сетки 
            {
                double* u_n = new double[N + 1]; // Последний найденный слой
                double* u_n1 = new double[N + 1];

                max1 = max2;
                // Из начальных данных находим нулевой временной слой
                for (int i = 0; i <= N; i++) {
                    u_n[i] = initial(a + i * h);

                }

                for (int m = 1; m <= M; m++) {
                    progonka(m, u_n, u_n1, tau, h, N, sigma,0);

                    for (int j = 0; j <= N; j++) {
                        u_n[j] = u_n1[j];
                    }

                }
                //Выводим последний временной слой и его модуль отклонения от точного решения
                for (int j = 0; j <= N; j++) {
                    Out1 << a + h * j << " " << u_n[j] << " " << abs(exactsln(M * tau, a + h * j) - u_n[j]) << endl;
                }
                Out1 << endl;
                max2 = NormOfDifference(u_n, N);
                if (set != 1) {
                    Out4 << "При N = " << N << "  при M = " << M << "   максимальное отклонение от точного = " << max2 << endl;
                    Out4 << "Численный порядок сходимости = " << log2(max1 / max2) << endl << endl;
                }
                else {
                    Out4 << "При N = " << N << "  При M = " << M << "   максимальное отклонение от точного = " << max2 << endl << endl;
                }
                delete[] u_n;
                delete[] u_n1;
                N = N * 2;
                M = M * 4;
                tau = T / M;
                l = b - a;
                h = l / N;
            }

        }
        else // Условие второго рода на левой границе и условие третьего рода на правой границе
        {
            for (int set = 1; set <= 5; set++) {
                double* u_n = new double[N + 1]; // Последний найденный слой
                double* u_n1 = new double[N + 1];
                max1 = max2;
                // Из начальных данных находим нулевой временной слой
                for (int i = 0; i <= N; i++) {
                    u_n[i] = initial(a + i * h);

                }

                for (int m = 1; m <= M; m++) {
                    progonkaKr(m, u_n, u_n1, tau, h, N, sigma, 0);

                    for (int j = 0; j <= N; j++) {
                        u_n[j] = u_n1[j];
                    }

                }
                for (int j = 0; j <= N; j++) {

                    Out2 << a + h * j << " " << u_n[j] << " " << abs(exactsln(M * tau, a + h * j) - u_n[j]) << endl;
                }
                Out2 << endl;
                max2 = NormOfDifference(u_n, N);
                if (set != 1) {
                    Out5 << "При N = " << N << "  при M = " << M << "   максимальное отклонение от точного = " << max2 << endl;
                    Out5 << "Численный порядок сходимости = " << log2(max1 / max2) << endl << endl;
                }
                else {
                    Out5 << "При N = " << N << "  При M = " << M << "   максимальное отклонение от точного = " << max2 << endl << endl;
                }
                delete[] u_n;
                delete[] u_n1;
                N = N * 2;
                M = M * 4;
                tau = T / M;
                l = b - a;
                h = l / N;
            }

        }

    }
    else // Схема повышенной точности
    {
        
        double max1 = 0; //Для отыскания максимального отклонения на слое
        double max2 = 0;

        if (Parametr2 == 0) // Условие Дирихле всюду
        {

            for (int set = 1; set <= 5; set++) // Измельчение сетки 
            {
                double sigma = 0.5 - h * h / (12 * tau);
                double* u_n = new double[N + 1]; // Последний найденный слой
                double* u_n1 = new double[N + 1];

                max1 = max2;
                // Из начальных данных находим нулевой временной слой
                for (int i = 0; i <= N; i++) {
                    u_n[i] = initial(a + i * h);

                }

                for (int m = 1; m <= M; m++) {
                    progonka(m, u_n, u_n1, tau, h, N, sigma, 1);

                    for (int j = 0; j <= N; j++) {
                        u_n[j] = u_n1[j];
                    }

                }
                //Выводим последний временной слой и его модуль отклонения от точного решения
                for (int j = 0; j <= N; j++) {
                    Out6 << a + h * j << " " << u_n[j] << " " << abs(exactsln(M * tau, a + h * j) - u_n[j]) << endl;
                }
                Out6 << endl;
                max2 = NormOfDifference(u_n, N);
                if (set != 1) {
                    Out7 << "При N = " << N << "  при M = " << M << "   максимальное отклонение от точного = " << max2 << endl;
                    Out7 << "Численный порядок сходимости = " << log2(max1 / max2) << endl << endl;
                }
                else {
                    Out7 << "При N = " << N << "  При M = " << M << "   максимальное отклонение от точного = " << max2 << endl << endl;
                }
                delete[] u_n;
                delete[] u_n1;
                N = N * 2;
                M = M * 4;
                tau = T / M;
                l = b - a;
                h = l / N;
            }
        }
        else // Условие второго рода на левой границе и условие третьего рода на правой границе
        {
            for (int set = 1; set <= 5; set++) {
                double sigma = 0.5 - h * h / (12 * tau);
                double* u_n = new double[N + 1]; // Последний найденный слой
                double* u_n1 = new double[N + 1];
                max1 = max2;
                // Из начальных данных находим нулевой временной слой
                for (int i = 0; i <= N; i++) {
                    u_n[i] = initial(a + i * h);
                }
                for (int m = 1; m <= M; m++) {
                    progonkaKr(m, u_n, u_n1, tau, h, N, sigma, 1);

                    for (int j = 0; j <= N; j++) {
                        u_n[j] = u_n1[j];
                    }

                }
                for (int j = 0; j <= N; j++) {

                    Out8 << a + h * j << " " << u_n[j] << " " << abs(exactsln(M * tau, a + h * j) - u_n[j]) << endl;
                }
                Out8 << endl;
                max2 = NormOfDifference(u_n, N);
                if (set != 1) {
                    Out9 << "При N = " << N << "  при M = " << M << "   максимальное отклонение от точного = " << max2 << endl;
                    Out9 << "Численный порядок сходимости = " << log2(max1 / max2) << endl << endl;
                }
                else {
                    Out9 << "При N = " << N << "  При M = " << M << "   максимальное отклонение от точного = " << max2 << endl << endl;
                }
                delete[] u_n;
                delete[] u_n1;
                N = N * 2;
                M = M * 4;
                tau = T / M;
                l = b - a;
                h = l / N;
            }
        }
    }
    
    
    
}
/*

Графики приближеных решений с задачей Дирихле всюду

plot 'C:/Users/Natalya Yuda/source/repos/112/Mat.txt' using 1:2 with lines title 'The approximate solution', 'C:/Users/Natalya Yuda/source/repos/112/Out3.txt' using 1:2 with lines title 'The exact solution'

Графики приближенных решений с условием второго рода на левой границе и условием третьего рода на правой границе

plot 'C:/Users/Natalya Yuda/source/repos/112/Out.txt' using 1:2 with lines title 'The approximate solution', 'C:/Users/Natalya Yuda/source/repos/112/Out3.txt' using 1:2 with lines title 'The exact solution'

Графики погрешностей для задачи Дирихле всюду

plot 'C:/Users/Natalya Yuda/source/repos/112/Mat.txt' using 1:3 with lines title 'The approximation error'

Графики погрешностей с условием второго рода на левой границе и условием третьего рода на правой границе

plot 'C:/Users/Natalya Yuda/source/repos/112/Out.txt' using 1:3 with lines title 'The approximation error'

plot 'C:/Users/Natalya Yuda/source/repos/112/Out8.txt' using 1:2 with lines title 'The approximate solution', 'C:/Users/Natalya Yuda/source/repos/112/Out3.txt' using 1:2 with lines title 'The exact solution'

plot 'C:/Users/Natalya Yuda/source/repos/112/Out8.txt' using 1:3 with lines title 'The approximation error'


*/