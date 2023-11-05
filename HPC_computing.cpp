#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <mpi.h>


double function(double x1, double x2, double x3) {
    // Здесь нужно вставить вашу функцию
    return x1 * x1 + x2 * x2 + x3 * x3;
}


int main(int argc, char** argv) {
    int method = atof(argv[1]);
    int n = atof(argv[2]);

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    auto start_time = std::chrono::high_resolution_clock::now();

    int local_n = n / size;
    double a, b, c, d, e, f;
    double result = 0.0;
    double total_result;

    /* Задаём пределы интегрирования */
    if (argc == 5) {
        a = atof(argv[3]); // Нижние пределы интеграла
        b = atof(argv[4]); // Верхние пределы интеграла
        
        c = e = a;
        d = f = b;
    }
    else if (argc == 9)
    {
        a = atof(argv[3]); // Нижние пределы интеграла
        b = atof(argv[4]); // Верхние пределы интеграла

        c = atof(argv[5]); // Нижние пределы интеграла
        d = atof(argv[6]); // Верхние пределы интеграла

        e = atof(argv[7]); // Нижние пределы интеграла
        f = atof(argv[8]); // Верхние пределы интеграла
    }
    else {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " a b c d e f " << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    if (method == 1) {
        /* Метод Монте-Карло */
        srand(time(NULL) + rank);
        for (int  i = rank * local_n; i < (rank + 1) * local_n; ++i) {
            double x1 = a + static_cast<double>(rand()) / RAND_MAX * (b - a);
            double x2 = c + static_cast<double>(rand()) / RAND_MAX * (d - c);
            double x3 = e + static_cast<double>(rand()) / RAND_MAX * (f - e);
            result += function(x1, x2, x3);
        }

        MPI_Reduce(&result, &total_result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        total_result *= (b - a) * (d - c) * (f - e) / n;
    }

    else if (method == 2) {
        /* Метод трапеции */
        double dx = (b - a) / n;
        double dy = (d - c) / n;
        double dz = (f - e) / n;

        for (int i = rank * local_n; i < (rank + 1) * local_n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    double x0 = a + i * dx;
                    double x1 = a + (i + 1) * dx;
                    double y0 = c + j * dy;
                    double y1 = c + (j + 1) * dy;
                    double z0 = e + k * dz;
                    double z1 = e + (k + 1) * dz;

                    result += 0.125 * (function(x0, y0, z0) + function(x1, y0, z0) + function(x0, y1, z0) + function(x1, y1, z0)
                        + function(x0, y0, z1) + function(x1, y0, z1) + function(x0, y1, z1) + function(x1, y1, z1));
                }
            }
        }

        MPI_Reduce(&result, &total_result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        total_result *= dx * dy * dz;
    }

    else if (method == 3) {
        /* Метод прямоугольников */
        double dx = (b - a) / n;
        double dy = (d - c) / n;
        double dz = (f - e) / n;

        for (int i = rank * local_n; i < (rank + 1) * local_n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    double x0 = a + i * dx;
                    double x1 = a + (i + 1) * dx;
                    double y0 = c + j * dy;
                    double y1 = c + (j + 1) * dy;
                    double z0 = e + k * dz;
                    double z1 = e + (k + 1) * dz;

                    result += function((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2);
                }
            }
        }

        MPI_Reduce(&result, &total_result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        total_result *= dx * dy * dz;

    }
    

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;


    if (rank == 0) {
        std::cout << total_result << std::endl;
        std::cout << elapsed.count() << std::endl;
    }

    MPI_Finalize();
    return 0;
}
