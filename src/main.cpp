#include <iostream>
#include <vector>
#include <tuple>
#include "../include/Matrix.hpp"

void showvec(std::vector<double> x)//ベクトル表示
{
    for(int i=0;i<(int)x.size();i++)
        std::cout << " " << x[i];
    std::cout << std::endl;
}

int main(void)
{
    std::vector<std::vector<double>> a=
    {
        { 1,-2,-1},
        {-2, 1,-1},
        {-1,-1, 0},
    };
    ///std::vector<double> b={-7,35,-19};
    Matrix A(a);
    std::cout << "A" << std::endl;
    A.show();
    ////行列式
    //double d=A.det();
    //std::cout << "det:" << d << std::endl;
    //if(d!=0)
    //{
    //    auto Ainv=A.inv();//逆行列
    //    std::cout << "A^-1" << std::endl;
    //    Ainv.show();
    //    auto x=Ainv*b;
    //    std::cout << "x" << std::endl;
    //    showvec(x);     
    //}
    //対角化
    auto [PAP,P] = A.culc_eigen();//固有値，固有ベクトル計算
    std::cout << "P^-1AP" << std::endl;
    PAP.show();
    std::cout << "P" << std::endl;
    P.show();
    //m.show();
}

//int main(void)
//{
//    std::vector<std::vector<double>> a=
//    {
//        {1,2,3},
//        {2,3,4},
//        {3,2,3},
//    };
//    std::vector<double> b={1,3,5};
//    Matrix A(a);
//    std::cout << "A" << std::endl;
//    A.show();
//    //行列式
//    double d=A.det();
//    std::cout << "det:" << d << std::endl;
//    //3次の連立1次方程式
//    /*
//        1*x0+2*x1+3x2=1
//        2*x0+3*x1+4x2=3
//        3*x0+2*x1+3x2=5
//        x={x0,x1,x2}
//    */
//    auto Ainv=A.inv();//逆行列
//    std::cout << "A^-1" << std::endl;
//    Ainv.show();
//    auto x=Ainv*b;
//    std::cout << "x" << std::endl;
//    showvec(x);
//    //対角化
//    auto [PAP,P] = A.culc_eigen();//固有値，固有ベクトル計算
//    std::cout << "P^-1AP" << std::endl;
//    PAP.show();
//    std::cout << "P" << std::endl;
//    P.show();
//    //m.show();
//    return 0;
//}