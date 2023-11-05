#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <omp.h>


double function(double x1, double x2, double x3) {
    // Здесь нужно вставить вашу функцию
    return x1 * x1 + x2 * x2 + x3 * x3;
}

int main(int argc, char** argv) {
    int method = atoi(argv[1]);
    int n = atoi(argv[2]);

    int threads = atoi(argv[3]); // Количество потоков OpenMP
    omp_set_num_threads(threads);

    double a, b, c, d, e, f;
    double total_result = 0.0;

    /* Задаём пределы интегрирования */
    if (argc == 6) {
        a = atof(argv[4]);
        b = atof(argv[5]);
        c = e = a;
        d = f = b;
    }
    else if (argc == 10) {
        a = atof(argv[4]);
        b = atof(argv[5]);
        c = atof(argv[6]);
        d = atof(argv[7]);
        e = atof(argv[8]);
        f = atof(argv[9]);
    }
    else {
        std::cerr << "Usage: " << argv[0] << " a b c d e f " << std::endl;
        return 1;
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    if (method == 1) {
        // Метод Монте-Карло
#pragma omp parallel reduction(+:total_result)
        {
            srand(time(NULL) + omp_get_thread_num());

#pragma omp for
            for (int i = 0; i < n; ++i) {
                double x1 = a + static_cast<double>(rand()) / RAND_MAX * (b - a);
                double x2 = c + static_cast<double>(rand()) / RAND_MAX * (d - c);
                double x3 = e + static_cast<double>(rand()) / RAND_MAX * (f - e);
                total_result += function(x1, x2, x3);
            }
        }

    total_result *= (b - a) * (d - c) * (f - e) / n;
    }

    else if (method == 2) {
        // Метод трапеций
        double dx = (b - a) / n;
        double dy = (d - c) / n;
        double dz = (f - e) / n;

#pragma omp parallel reduction(+:total_result)
        {
#pragma omp for
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    for (int k = 0; k < n; ++k) {
                        double x0 = a + i * dx;
                        double x1 = a + (i + 1) * dx;
                        double y0 = c + j * dy;
                        double y1 = c + (j + 1) * dy;
                        double z0 = e + k * dz;
                        double z1 = e + (k + 1) * dz;

                        total_result += 0.125 * (function(x0, y0, z0) + function(x1, y0, z0) + function(x0, y1, z0) + function(x1, y1, z0)
                            + function(x0, y0, z1) + function(x1, y0, z1) + function(x0, y1, z1) + function(x1, y1, z1));
                    }
                }
            }
        }

        total_result *= dx * dy * dz;
    }

    else if (method == 3) {
        // Метод трапеций
        double dx = (b - a) / n;
        double dy = (d - c) / n;
        double dz = (f - e) / n;

#pragma omp parallel reduction(+:total_result)
        {
#pragma omp for
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    for (int k = 0; k < n; ++k) {
                        double x0 = a + i * dx;
                        double x1 = a + (i + 1) * dx;
                        double y0 = c + j * dy;
                        double y1 = c + (j + 1) * dy;
                        double z0 = e + k * dz;
                        double z1 = e + (k + 1) * dz;

                        total_result += function((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2);
                    }
                }
            }
        }

        total_result *= dx * dy * dz;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    std::cout << total_result << std::endl;
    std::cout << elapsed.count() << std::endl;

    return 0;
}
