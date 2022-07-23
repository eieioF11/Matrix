#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>

class Matrix
{
    private:
        const double calclim = 0.0000000001;
        const int rooplim = 100;
        std::vector<std::vector<double>> data;
        std::vector<double> get_column(std::vector<std::vector<double>> x,int n);
        double dot(std::vector<double> x, std::vector<double> y);//内積
        std::vector<double> div(std::vector<double> x, double y);//各要素の除算
        std::pair<Matrix,Matrix> crout(std::vector<std::vector<double>> a);
        // 前進代入
        std::pair<std::vector<double>,std::vector<double>> forward_substitution(std::vector<std::vector<double>> a, std::vector<double> y, std::vector<double> b)
        {
            for (int row = 0; row < (int)a.size(); row++)
            {
                for (int col = 0; col < row; col++)
                    b[row] -= a[row][col] * y[col];
                y[row] = b[row];
            }
            return std::make_pair(y,b);
        }
        // 後退代入
        std::pair<std::vector<double>,std::vector<double>> backward_substitution(std::vector<std::vector<double>> a, std::vector<double> x, std::vector<double> y)
        {
            for (int row = (int)a.size() - 1; row >= 0; row--)
            {
                for (int col = (int)a.size() - 1; col > row; col--)
                    y[row] -= a[row][col] * x[col];
                if(a[row][row]!=0)
                    x[row] = y[row] / a[row][row];
                else
                    x[row] = DBL_MAX;
            }
            return std::make_pair(x,y);
        }
        std::vector<double> inverse(std::vector<std::vector<double>> a, std::vector<double> x0)
        {
            const int N = (int)a.size();
            // 正規化 (ベクトル x0 の長さを１にする)
            x0=normarize(x0);
            double e0 = 0.0;
            for (int i = 0; i < N; i++)
                e0 += x0[i];

            for (int k = 1; k <= rooplim; k++)
            {

                //前進代入
                std::vector<double> b(N,0.0);
                std::vector<double> y(N,0.0);
                for (int i = 0; i < N; i++)
                    b[i] = x0[i];
                auto yb = forward_substitution(a,y,b);
                y=yb.first;
                b=yb.second;
                //後退代入
                std::vector<double> x1(N,0.0);
                auto xy=backward_substitution(a,x1,y);
                x1=xy.first;
                y=xy.second;

                // 正規化 (ベクトル x1 の長さを１にする)
                x1=normarize(x1);
                // 収束判定
                double e1 = 0.0;
                for (int i = 0; i < N; i++)
                    e1 += x1[i];
                if (fabs(e0 - e1) < 0.00000000001) break;

                for (int i = 0; i < N; i++)
                    x0[i] = x1[i];
                e0 = e1;
            }
            return x0;
        }

        std::vector<double> normarize(std::vector<double> x)
        {
            double s = 0.0;
            for (int i = 0; i < (int)x.size(); i++)
                s += x[i] * x[i];
            s = std::sqrt(s);
            for (int i = 0; i < (int)x.size(); i++)
                x[i] /= s;
            return x;
        }
        void showvec(std::vector<double> x)
        {
            for(int i=0;i<(int)x.size();i++)
                std::cout << " " << x[i];
            std::cout << std::endl;
        }
    public:
        void test()
        {
            culc_eigen();
            double d=det();
            std::cout << "det:" << d << std::endl;
        }
        Matrix(int n);
        Matrix(std::vector<std::vector<double>> data);
        void add_value(int i,int j,double value);
        double get_value(int i,int j);
        Matrix& operator = (std::vector<std::vector<double>> val){this->data=val;return *this;};
        Matrix& operator = (Matrix val){this->data=val;return *this;};
        //行列演算
        Matrix operator + (Matrix x);//和
        Matrix operator * (Matrix x);//積
        friend std::vector<double> operator * (Matrix x,std::vector<double> vec);//ベクトルとの積
        Matrix inv();//逆行列
        double det();//行列式
        std::vector<double> equations(Matrix x,std::vector<double> vec);
        std::pair<Matrix,Matrix> crout();//LU分解 <L,R>
        std::pair<Matrix,Matrix> culc_eigen();//固有値,固有ベクトルの計算
        //値取得
        operator std::vector<std::vector<double>>(){return data;};
        std::vector<double> get_row(int row);
        std::vector<double> get_line(int line);
        void show();
        void show(std::vector<std::vector<double>> a);
};
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
Matrix Matrix::operator + (Matrix x)
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
Matrix Matrix::operator * (Matrix x)
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

std::vector<double> operator * (Matrix x,std::vector<double> vec)
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

Matrix Matrix::inv()
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
double Matrix::det()
{
  std::vector<std::vector<double>> A(this->data.size(), std::vector<double>(this->data[0].size()));
  A = this->data;
  int n=A.size();
  if(n==1)
    return A[0][0];
  else if(n==2)
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]; //要素数２までは簡単なので直接計算

  //0行目で余因子展開
  int sum = 0;
  for(int i=0;i<n;i++){ //A[0][i]で余因子展開する
    std::vector<std::vector<double>> B(A.size(), std::vector<double>(A[0].size()));
    B = A;
    for(int j=0;j<n;j++){
      B[j].erase(B[j].begin()+i);//B[j][i]を消す
    }
    B.erase(B.begin()); //0行目は最後に消す
    sum += A[0][i] * pow(-1,i+2) * ((Matrix)B).det(); //Bが余因子行列になっている
  }
  return sum;
}
std::pair<Matrix,Matrix> Matrix::crout(std::vector<std::vector<double>> a)
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
    
    return std::make_pair(L,R);
}
std::pair<Matrix,Matrix> Matrix::crout()
{
    return crout(this->data);
}
std::vector<double> Matrix::equations(Matrix x,std::vector<double> vec)
{
    auto LR=x.crout();
    auto L=LR.first;
    auto R=LR.second;
    x=(R*L);
    return inverse((std::vector<std::vector<double>>)x,vec);
}
std::pair<Matrix,Matrix> Matrix::culc_eigen()
{
    Matrix m(this->data);
    for (int k=0;k<this->rooplim;k++)
    {
        auto LR=m.crout();
        auto L=LR.first;
        auto R=LR.second;
        m=(R*L);
        int c=0;
        for (int j=0;j<(int)(this->data.size()-1);j++)
        {
            for(int i=j;i<(int)(this->data.size()-1);i++)
                if(std::fabs(m.get_value(i+1,j))<this->calclim)
                    c++;
        }
        if(c>=(int)this->data.size())
            break;
    }  
    Matrix PAP(this->data.size());
    //対角成分の抜き取り
    for (int i=0;i<(int)this->data.size();i++)
        PAP.add_value(i,i,m.get_value(i,i));
    std::cout << "P^-1AP" << std::endl;
    PAP.show();//対角化行列
    for (int i=0;i<(int)this->data.size();i++)
    {
        std::cout << "A" << std::endl;
        Matrix A(this->data);
        for (int j=0;j<(int)this->data.size();j++)
            A.add_value(j,j,A.get_value(j,j)-PAP.get_value(i,i));
        A.show();
        std::vector<double> v(this->data.size(),1.0);
        auto x=equations(A,v);
        std::cout << "x" << std::endl;
        showvec(x);
    }
    return std::make_pair(PAP,PAP);
}
void Matrix::add_value(int i,int j,double value)
{
    this->data[i][j]=value;
}
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

int main(void)
{
    //std::vector<std::vector<double>> data = 
    //{
    //    {5.0, 4.0, 1.0, 1.0},
    //    {4.0, 5.0, 1.0, 1.0},
    //    {1.0, 1.0, 4.0, 2.0},
    //    {1.0, 1.0, 2.0, 4.0}
    //};
    std::vector<std::vector<double>> data=
    {
        {1,-2,0},
        {-2,2,-2},
        {0,-2,3},
    };
    //std::vector<std::vector<double>> data=
    //{
    //    {1,2,3},
    //    {4,5,6},
    //    {7,8,9},
    //};
    Matrix m(data);
    m.test();
    //m.show();
    return 0;
}