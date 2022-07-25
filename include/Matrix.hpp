/* 正方行列演算ライブラリ
 * ファイル名:Matrix.hpp
 * 使用Cコンパイラ:gcc version 9.4.0
 */
#ifndef MATRIX_H_
#define MATRIX_H_
#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <float.h>

class Matrix
{
    private:
        const double calclim = 0.000000000001;//計算精度設定
        const int rooplim = 100;//ループ回数制限設定
        std::vector<std::vector<double>> data;//行列
        std::vector<double> get_column(std::vector<std::vector<double>> x,int n);
        double dot(std::vector<double> x, std::vector<double> y);//内積
        std::vector<double> div(std::vector<double> x, double y);//各要素の除算
        //LR
        std::tuple<Matrix,Matrix> crout(std::vector<std::vector<double>> a);
        // ハウスホルダー変換
        std::tuple<std::vector<std::vector<double>>, std::vector<double>,  std::vector<double>> tridiagonalize(std::vector<std::vector<double>> a, std::vector<double> d, std::vector<double> e);
        // QR分解
        std::tuple<std::vector<std::vector<double>>, std::vector<double>,  std::vector<double>> decomp(std::vector<std::vector<double>> a, std::vector<double> d, std::vector<double> e);
    public:
        Matrix(int n);//要素0で初期化
        Matrix(std::vector<std::vector<double>> data);//要素をdataの値で初期化
        void add_value(int i,int j,double value);//i行j列にvalueを代入
        double get_value(int i,int j);//i行j列のデータ抜き出し
        Matrix& operator = (std::vector<std::vector<double>> val){this->data=val;return *this;};
        Matrix& operator = (Matrix val){this->data=val;return *this;};
        //行列演算
        Matrix operator + (Matrix x);//和
        Matrix operator * (Matrix x);//積
        friend std::vector<double> operator * (Matrix x,std::vector<double> vec);//ベクトルとの積
        Matrix inv();//逆行列
        long double det();//行列式
        std::tuple<Matrix,Matrix> crout();//LU分解 <L,R>
        std::tuple<Matrix,Matrix> culc_eigen();//固有値,固有ベクトルの計算 tupleで対角成分が固有値の行列と固有ベクトルを列ベクトルにした行列を返す
        //値取得
        operator std::vector<std::vector<double>>(){return data;};
        std::vector<double> get_row(int row);//行の抜き出し
        std::vector<double> get_line(int line);//列の抜き出し
        void show();//行列の表示
        void show(std::vector<std::vector<double>> a);
};
//初期化部分
Matrix::Matrix(int n)//n行n列の行列作成
{
    std::vector<std::vector<double>> data(n,std::vector<double>(n,0));
    this->data=data;
}
Matrix::Matrix(std::vector<std::vector<double>> data)//n行n列の行列作成
{
    if(data.size()==data[0].size())
    {
        this->data=data;
    }
    else
        std::cerr << "行と列のサイズが異なっています" << std::endl;
}
//行列表示
void Matrix::show(std::vector<std::vector<double>> a)
{
    for (int i = 0; i < (int)a.size(); i++)
    {
        for(int j = 0; j < (int)a.size(); j ++)
            std::cout << " " << a[i][j];
        std::cout << std::endl;
    } 
}
void Matrix::show()
{
    show(this->data);
}
//行列演算
Matrix Matrix::operator + (Matrix x)//行列同士の加算
{
    std::vector<std::vector<double>> A=this->data;
    std::vector<std::vector<double>> B=x.data;
    std::vector<std::vector<double>> C;
    C = A;

    for(int i=0; i < (int)A.size(); i++)
        for(int j=0; j < (int)A[i].size(); j++)
            C[i][j] = A[i][j] + B[i][j];
    Matrix y(C);
    return y;
}
Matrix Matrix::operator * (Matrix x)//行列同士の乗算
{
    std::vector<std::vector<double>> A=this->data;
    std::vector<std::vector<double>> B=x.data;
    std::vector<std::vector<double>> C(A.size(), std::vector<double>(B[0].size(),0));

    for(int i=0; i < (int)A.size(); i++)
        for(int s=0;s<(int)B[0].size();s++)
            C[i][s] = dot(A[i],get_column(B,s));
    Matrix y(C);
    return y;
}

std::vector<double> operator * (Matrix x,std::vector<double> vec)//ベクトルとの乗算
{
    std::vector<std::vector<double>> X = x.data;
    std::vector<double> C(x.data.size(),0);
    for(int i=0; i < (int)x.data.size(); i++)
        C[i] = x.dot(X[i],vec);
    return C;
}

std::vector<double> Matrix::get_column(std::vector<std::vector<double>> x,int n)
{
    std::vector<double> y(x.size());
    for(int i=0;i<(int)x.size();i++)
        y[i] = x[i][n];
    return y;
}
double Matrix::dot(std::vector<double> x, std::vector<double> y)
{
    double z = 0;
    for(int i=0;i<(int)x.size();i++)
        z += x[i] * y[i];
    return z;
}
std::vector<double> Matrix::div(std::vector<double> x, double y)
{
    for(int i=0;i<(int)x.size();i++)
        x[i] = x[i] / y;
    return x;
}

Matrix Matrix::inv()//逆行列
{
    std::vector<std::vector<double>> A(data.size(), std::vector<double>(data[0].size()*2,0));//A[N][N*2]の行列を作る
    int n = data.size();

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            A[i][j] = data[i][j]; //Aの前半列にはオブジェクトの行列を格納
            A[i][j+n] = (i==j)? 1:0; //Aの後半列には単位行列を格納
        }
    }

    for(int i=0;i<n;i++)
    {
        A[i] = div(A[i],A[i][i]); //Aのi行目についてA[i][i]の値で割る
        for(int j=0;j<n;j++)
        {
            if(i==j)
                continue;
            double t = A[j][i];
            for(int k=0;k<n*2;k++)
                A[j][k] = A[j][k] - A[i][k] * t; //Aのj行目についてA[j][i]=0となるようにAのi行目の定数倍を引く
        }
    }
    //この時点でAの前半は単位行列になっている
    std::vector<std::vector<double>> B(n, std::vector<double>(n));
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
        B[i][j] = A[i][j+n]; //Aの後半列を結果行列に格納

    return (Matrix)B;
}
long double Matrix::det()//行列式
{
    std::vector<std::vector<double>> A(this->data.size(), std::vector<double>(this->data[0].size()));
    A=this->data;
    int n=A.size();
    if(n==1)
        return A[0][0];
    else if(n==2)
        return A[0][0] * A[1][1] - A[0][1] * A[1][0];
    else if(n==3)
        return A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[0][2] * A[1][0] * A[2][1]
            - A[0][2] * A[1][1] * A[2][0] - A[0][1] * A[1][0] * A[2][2] - A[0][0] * A[1][2] * A[2][1];

    //三角行列を作成
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(i<j)
            {
                double buf;
                if(A[i][i]!=0)
                    buf=A[j][i]/A[i][i];
                else
                    buf=A[j][i]/DBL_MAX;
                for(int k=0;k<n;k++)
                    A[j][k]-=A[i][k]*buf;
            }
        }
    }
    
    //対角部分の積
    long double det=1.0;
    for(int i=0;i<n;i++)
        det*=A[i][i];
    return det;
}
//分解LR
std::tuple<Matrix,Matrix> Matrix::crout(std::vector<std::vector<double>> a)
{
    
    const int N = a.size();
    Matrix L(N);
    Matrix R(N);
    L.add_value(0,0,1.0);
    R.add_value(0,0,a[0][0]);
    for(int j=1;j<=N-1;j++)
        R.add_value(0,j,a[0][j]);
    for(int i=1;i<=N-1;i++)
        R.add_value(i,0,0.0);
    for(int i=1;i<=N-1;i++)
        L.add_value(i,0,a[i][0]/R.get_value(0,0));
    for(int j=1;j<=N-1;j++)
        L.add_value(0,j,0.0);
    for(int i=1;i<=N-1;i++)
    {
        L.add_value(i,i,1.0);
        R.add_value(i,i,a[i][i]);
        for(int k=0;k<=i-1;k++)
            R.add_value(i,i,R.get_value(i,i)-(L.get_value(i,k)*R.get_value(k,i)));
        for(int j=i+1;j<=N-1;j++)
        {
            L.add_value(j,i,a[j][i]);
            for(int k=0;k<=i-1;k++)
                L.add_value(j,i,L.get_value(j,i)-(L.get_value(j,k)*R.get_value(k,i)));
            L.add_value(j,i,L.get_value(j,i)/R.get_value(i,i));
            R.add_value(j,i,0.0);
            L.add_value(i,j,0.0);
            R.add_value(i,j,a[i][j]);
            for(int k=0;k<=i-1;k++)
                R.add_value(i,j,R.get_value(i,j)-(L.get_value(i,k)*R.get_value(k,j)));
        }
    }
    
    return {L,R};
}
std::tuple<Matrix,Matrix> Matrix::crout()
{
    return crout(this->data);
}
// ハウスホルダー変換
std::tuple<std::vector<std::vector<double>>, std::vector<double>,  std::vector<double>> Matrix::tridiagonalize(std::vector<std::vector<double>> a, std::vector<double> d, std::vector<double> e)
{
    const int N = a.size();
    std::vector<double> v(N,0.0);

    for (int k = 0; k < N - 2; k++)
    {
        std::vector<double> w(N,0.0);
        d[k] = a[k][k];

        double t = 0.0;
        for (int i = k + 1; i < N; i++)
        {
            w[i] =  a[i][k];
            t += w[i] * w[i];
        }
        double s = std::sqrt(t);
        if (w[k + 1] < 0)
            s = -s;

        if (std::fabs(s) < this->calclim)
            e[k + 1] = 0.0;
        else
        {
            e[k + 1]  = -s;
            w[k + 1] +=  s;
            s = 1 / std::sqrt(w[k + 1] * s);
            for (int i = k + 1; i < N; i++)
                w[i] *= s;

            for (int i = k + 1; i < N; i++)
            {
                s = 0.0;
                for (int j = k + 1; j < N; j++)
                {
                    if (j <= i)
                        s += a[i][j] * w[j];
                    else
                        s += a[j][i] * w[j];
                }
                v[i] = s;
            }

            s = 0.0;
            for (int i = k + 1; i < N; i++)
                s += w[i] * v[i];
            s /= 2.0;
            for (int i = k + 1; i < N; i++)
                v[i] -= s * w[i];
            for (int i = k + 1; i < N; i++)
                for (int j = k + 1; j <= i; j++)
                    a[i][j] -= w[i] * v[j] + w[j] * v[i];
            for (int i = k + 1; i < N; i++)
                a[i][k] = w[i];
        }
    }

    d[N - 2] = a[N - 2][N - 2];
    d[N - 1] = a[N - 1][N - 1];

    e[0]     = 0.0;
    e[N - 1] = a[N - 1][N - 2];

    for (int k = N - 1; k >= 0; k--)
    {
        std::vector<double> w(N,0.0);
        if (k < N - 2)
        {
            for (int i = k + 1; i < N; i++)
                w[i] = a[i][k];
            for (int i = k + 1; i < N; i++)
            {
                double s = 0.0;
                for (int j = k + 1; j < N; j++)
                    s += a[i][j] * w[j];
                v[i] = s;
            }
            for (int i = k + 1; i < N; i++)
                for (int j = k + 1; j < N; j++)
                    a[i][j] -= v[i] * w[j];
        }
        for (int i = 0; i < N; i++)
            a[i][k] = 0.0;
        a[k][k] = 1.0;
    }
    return {a,d,e};
}
// QR分解
std::tuple<std::vector<std::vector<double>>, std::vector<double>,  std::vector<double>> Matrix::decomp(std::vector<std::vector<double>> a, std::vector<double> d, std::vector<double> e)
{
    const int N = a.size();
    e[0] = 1.0;
    int h = N - 1;
    while (fabs(e[h]) < this->calclim) h--;

    while (h > 0)
    {
        e[0] = 0.0;
        int l = h - 1;
        while (std::fabs(e[l]) >= this->calclim) l--;

        for (int j = 1; j <= this->rooplim; j++)
        {
            double w = (d[h - 1] - d[h]) / 2.0;
            double s = std::sqrt(w * w + e[h] * e[h]);
            if (w < 0.0)
                s = -s;

            double x = d[l] - d[h] + e[h] * e[h] / (w + s);
            double y = e[l + 1];
            double z = 0.0;
            double t;
            double u;
            for (int k = l; k < h; k++)
            {
                if (std::fabs(x) >= std::fabs(y))
                {
                    t = -y / x;
                    u = 1 / std::sqrt(t * t + 1.0);
                    s = t * u;
                }
                else
                {
                    t = -x / y;
                    s = 1 / std::sqrt(t * t + 1.0);
                    if (t < 0)
                        s = -s;
                    u = t * s;
                }
                w = d[k] - d[k + 1];
                t = (w * s + 2 * u * e[k + 1]) * s;
                d[k    ] = d[k    ] - t;
                d[k + 1] = d[k + 1] + t;
                e[k    ] = u * e[k] - s * z;
                e[k + 1] = e[k + 1] * (u * u - s * s) + w * s * u;

                for (int i = 0; i < N; i++)
                {
                    x = a[k    ][i];
                    y = a[k + 1][i];
                    a[k    ][i] = u * x - s * y;
                    a[k + 1][i] = s * x + u * y;
                }

                if (k < N - 2)
                {
                    x = e[k + 1];
                    y = -s * e[k + 2];
                    z = y;
                    e[k + 2] = u * e[k + 2];
                }
            }

            // 収束判定
            if (std::fabs(e[h]) < this->calclim) break;
        }

        e[0] = 1.0;
        while (std::fabs(e[h]) < this->calclim) h--;
    }

    for (int k = 0; k < N - 1; k++)
    {
        int l = k;
        for (int i = k + 1; i < N; i++)
            if (d[i] > d[l])
                l = i;

        double t = d[k];
        d[k] = d[l];
        d[l] = t;

        for (int i = 0; i < N; i++)
        {
            t       = a[k][i];
            a[k][i] = a[l][i];
            a[l][i] = t;
        }
    }
    return {a,d,e};
}
std::tuple<Matrix,Matrix> Matrix::culc_eigen()//固有値 固有ベクトル
{
    std::vector<std::vector<double>> a(this->data);
    std::vector<double> d(this->data.size(),0.0);
    std::vector<double> e(this->data.size(),0.0);
    // ハウスホルダー変換
    auto [A, D, E] = tridiagonalize(a, d, e);
    // QR分解
    auto [x,eval,t] = decomp(A, D, E);
    Matrix PAP(this->data.size());
    for (int i=0;i<(int)this->data.size();i++)
        PAP.add_value(i,i,eval[i]);
    std::vector<std::vector<double>> p;
    for (int i = 0; i < (int)x.size(); i++)
    {
        std::vector<double> cut;
        for (int j = 0; j < (int)x.size(); j++)
            cut.push_back(x[j][i]);
        p.push_back(cut);
    }
    Matrix P(p);
    return {PAP,P};
}
//要素の設定
void Matrix::add_value(int i,int j,double value)
{
    this->data[i][j]=value;
}
//行列の要素取得
double Matrix::get_value(int i,int j)//i行j列のデータ抜き出し
{
    return this->data[i][j];
}
std::vector<double> Matrix::get_row(int row)//row行目のデータ抜き出し
{
    return this->data[row];
}
std::vector<double> Matrix::get_line(int line)//line列目のデータ抜き出し
{
    std::vector<double> cut;
    for (int i = 0; i < (int)this->data.size(); i++)
        cut.push_back(this->data[i][line]);
    return cut;
}
#endif
