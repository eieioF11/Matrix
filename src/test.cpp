#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <tuple>
#include "../include/Matrix.hpp"
#include <eigen3/Eigen/Dense>
#include "../include/matplotlibcpp.h"
//using namespace Eigen;

void showvec(std::vector<double> x)//ベクトル表示
{
    for(int i=0;i<(int)x.size();i++)
        std::cout << " " << x[i];
    std::cout << std::endl;
}

constexpr int LENGTH = 4;
constexpr int WIDTH = 6;

std::tuple<std::vector<std::vector<double>>, std::vector<double>> randomMatrix(int n)
{
    const int randlim = 100;
    std::vector<std::vector<double>> a(n,std::vector<double>(n,0.0));
    std::srand(std::time(nullptr));
    for (auto &item : a)
        for (auto &i : item)
            i = (rand() % randlim);
    std::vector<double> b(n,0);
    for(int i=0;i<(int)b.size();i++)
        b[i] = i;
    return {a,b};
}

const static int LOOPMAX = 1000;
const static double TIMEOUT = 10000.0;

int main(void)
{
    std::vector<double> dim;
    std::vector<double> det_t;
    std::vector<double> equ_t;
    std::vector<double> eig_t;
    std::chrono::system_clock::time_point start, end;
    for(int i=2;i<LOOPMAX;i++)
    {
        auto [a,b] = randomMatrix(i);
        Matrix A(a);
        //A.show();
        //行列式
        start = std::chrono::system_clock::now();
        long double d=A.det();
        end = std::chrono::system_clock::now();
        double time1 = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        //std::cout << i << "次 det:" << d << " t:" << time << "[s]" << std::endl;
        if(d==0)
        {
            i--;
            continue;
        }
        if (time1 > TIMEOUT)
        {
            std::cerr << "time out :" << time1 << "[s]" << std::endl;
            break;
        }

        start = std::chrono::system_clock::now();
        auto Ainv=A.inv();//逆行列
        auto x=Ainv*b;
        end = std::chrono::system_clock::now();
        double time2 = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        //std::cout << i << "次 t:" << time2 << "[s]" << std::endl;
        //std::cout << "A^-1" << std::endl;
        //Ainv.show();
        //std::cout << "x" << std::endl;
        //showvec(x);
        if (time2 > TIMEOUT)
        {
            std::cerr << "time out :" << time2 << "[s]" << std::endl;
            break;
        }

        start = std::chrono::system_clock::now();
        auto [PAP,P] = A.culc_eigen();//固有値，固有ベクトル計算
        end = std::chrono::system_clock::now();
        double time3 = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        //std::cout << i << "次 t:" << time3 << "[s]" << std::endl;
        //std::cout << "P^-1AP" << std::endl;
        //PAP.show();
        //std::cout << "P" << std::endl;
        //P.show();
        if (time3 > TIMEOUT)
        {
            std::cerr << "time out :" << time3 << "[s]" << std::endl;
            break;
        }

        det_t.push_back(time1);
        equ_t.push_back(time2);
        eig_t.push_back(time3);

        dim.push_back(i);
        std::cout << i << "次 det:" << d << std::endl;
    }
    namespace plt = matplotlibcpp;
    plt::grid(true);
    plt::plot(dim, det_t, "r");
    plt::plot(dim, equ_t, "b");
    plt::plot(dim, eig_t, "g");
    plt::show();
    plt::subplot(1, 3, 1);
    plt::grid(true);
    plt::plot(dim, det_t, "r");
    plt::subplot(1, 3, 2);
    plt::grid(true);
    plt::plot(dim, equ_t, "b");
    plt::subplot(1, 3, 3);
    plt::grid(true);
    plt::plot(dim, eig_t, "g");
    plt::show();
    return 0;
}