//单纯形法.cpp
//说明：输入文件中问题需化为标准型(求最大值)
#include <iostream>
using namespace std;
#include <cmath>
#include <iomanip>
#include <fstream>
#include <ctime>

//全局变量-----------------------------------------------------------------------------------------------------------------------------------
const float M=pow(10.0,5);
int row=0,col=0,i,j,*R,*S,counter = 0;//counter：人工变量个数
float **A,*C,*b,*tempb,**e,*CB,*test,*rep,*y,**B,**tempB,*beta;

//函数---------------------------------------------------------------------------------------------------------------------------------------
void newspace();   //动态开辟数组
void delspace(float **X,int n);//释放数组空间
void init(float **X,int m,int n);//初始化数组
void output();//在屏幕输出
int maxof(float *X,int m);//求最大数
int minof(float *X,int m);//求最小正数
bool isAllNegative(float *X,int m);//判断数组是否非正数
void match();//与单位阵匹配获得初始基下标或是添加人工变量
bool isNoSovle();//判断是否无解
bool isNoLimitSolve();//判断是否为无界解
bool isMuchSolves();//判断是否有无穷多解
//------------------------------------------------------------------------------------------------------------------------------------------

int main()
{
    time_t begin,end;  //计算程序运行时间
    begin = clock();

    int subr,subc;
    float temp=0.0;

    //定义读出文件
    ifstream inFile("problem.txt",ios::in);
    if (!inFile)
    {
        cout<<"Could not open the file.\n";
        exit(0);
    }

    //定义写入文件
    ofstream outFile("result.txt",ios::out);
    if (!outFile)
    {
        cout<<"Could not open the file.\n";
        exit(0);
    }

    //从文件中读入条件
    inFile>>row>>col;

    newspace();
    C = new float [col+row];
    b = new float [row];

    for(i=0;i<col;i++)
        inFile>>C[i];
    for(i=0;i<row;i++)
        inFile>>b[i];
    for (i=0;i<row;i++)
        for (j=0;j<col;j++)
            inFile>>A [i][j];

    //设置单位矩阵，用以实现match函数
    e = new float *[row];
    for (i=0;i<row;i++)
        e[i]=new float [row];
    init(e,row,row);
    for (i=0;i<row;i++)
        e[i][i] = 1;


    B = new float *[row];//初始可行基的逆矩阵
    for (i=0;i<row;i++)
        B[i]=new float [row];
    init(B,row,row);

    tempB = new float *[row];
    for (i=0;i<row;i++)
        tempB[i]=new float [row];
    init(tempB,row,row);

    tempb = new float [row];
    for (i=0;i<row;i++)
    {
        tempb[i]=0;
    }

    y = new float [row];
    beta = new float [row];
    CB = new float [row];//基变量的价值系数
    R = new int [row];//基变量下标
    S = new int [col+row];//非基变量下标
    test = new float [col+row];//检验数
    rep = new float [row];//theta值

    //确定初始基变量和系数
    match();
    for (i=0;i<col+counter;i++)
        S[i]=i+1;
    for (i=0;i<row;i++)
        CB[i]=C[R[i]-1];
    for (i=0;i<row;i++)
    {
        for (j=0;j<row;j++)
        {
            B[i][j]=e[i][j];
        }
    }

    //计算初始基可行解
    for (i=0;i<row;i++)
    {
        for (j=0;j<row;j++)
        {
            tempb[i]+= B[i][j]*b[j];
        }
    }
    for (i=0;i<row;i++)
    {
        b[i]= tempb[i];
    }

    //求检验数
    for (i=0;i<col+counter;i++)
    {
        temp = 0;

        for (j=0;j<row;j++)
        {
            temp += CB[j]*A[j][i];
        }
        test[i] = C[i]-temp;

    }
    subc = maxof(test,col+counter);   //获得检验数的最大值下标

    //计算beta向量
    for (i=0;i<row;i++)
    {
        for (j=0;j<row;j++)
        {
            beta[i]+= B[i][j]*A[j][subc];
        }
    }

    cout<<"初始单纯形表:\n";
    output();
//--------------算法循环开始--------------//
    while(!isAllNegative(test,col+counter))
    {
        //计算基变量的逆矩阵B乘以资源向量b
        for (i=0;i<row;i++)
        {
            temp=0;

            for (j=0;j<row;j++)
            {
                temp += B[i][j]*beta[j];
            }
            tempb[i] = temp;

        }
        //求θ值
        for(i=0;i<row;i++)
        {
            if (beta[i] > 0)
            {
                rep[i] = tempb[i]/beta[i];
            }
            else
                rep[i] = -1;
        }
        subr=minof(rep,row);//获得θ最小值对应的下标

        //处理下标（换基）
        R[subr] = S[subc];
        for (i=0;i<row;i++)
            CB[i]=C[R[i]-1];

        //计算新基的逆矩阵
        temp=beta[subr];
        beta[subr]=1/temp;
        for (i=0;i<row;i++)
        {
            if (i!=subr)
            {
                y[i]=-beta[i]/temp;
            }
        }

        for (i=0;i<row;i++)
        {
            for (j=0;j<row;j++)
            {
                tempB[i][j]=1*B[i][j]+y[i]*B[subc][j];
            }
        }
        for (i=0;i<row;i++)
        {
            for (j=0;j<row;j++)
            {
                B[i][j]=tempB[i][j];
            }
        }

        //求检验数
        for (i=0;i<col+counter;i++)
        {
            temp = 0;
            for (j=0;j<row;j++)
                temp += CB[j]*A[j][i];
            test[i] = C[i]-temp;
        }
        subc=maxof(test,col+counter);   //获得检验数的最大值下标

        //判断是否为无界解
        if (isNoLimitSolve())
        {
            outFile<<"线性方程具有无界解。\n";
            break;//若为无界解，则跳出while循环
        }
    }

    if (!isNoLimitSolve())
    {
        if (isNoSovle())//判断无解
        {
            outFile<<"线性方程无可行解。\n";
        }

        //存在可行解
        else
        {
            //求值
            temp=0;
            for (i=0;i<row;i++)
            {
                temp += b[i]*CB[i];
            }//for

            outFile<<"问题存在最优解，最优解为: \n";
            for(i=0;i<row;i++)
                outFile<<"x"<<R[i]<<'='<<b[i]<<'\t';
            outFile<<"\n目标函数值为：z="<</*setprecision(8)<<*/temp;

            if (isMuchSolves())//判断是否为无穷多最优解
            {
                outFile<<"\n非基检验数中含0，问题具有无穷多最优解。\n";
            }//if
        }//else
    }//if
//-----------算法循环终止---------------//
    cout<<"最终单纯形表:\n";
    output();

    //释放内存空间
    delspace(A,row);
    delspace(e,row);
    delete[]C;
    delete[]b;
    delete[]CB;
    delete[]test;
    delete[]rep;
    delete[]R;
    delete[]S;

    cout<<"感谢使用本程序，你要求解的方程结果已存入result.txt\n";

    end = clock();
    outFile<<"\nruntime:   "<<float(end-begin)/CLOCKS_PER_SEC<<'s'<<endl;//计算程序运行时间(文件)
    cout<<"\nruntime:   "<<float(end-begin)/CLOCKS_PER_SEC<<'s'<<endl;//计算程序运行时间(屏幕)

    return 0;
}

void newspace()   //动态开辟数组
{
    A = new float *[row];
    if (A==NULL)
    {
        cout<<"开辟数组不成功！程序结束！"<<endl;
        return ;
    }
    for (i=0;i<row;i++)
    {
        A[i]=new float[row+col];
        if(A[i]==NULL)
        {
            cout<<"开辟数组不成功！程序结束！"<<endl;
            return ;
        }
    }
}

void delspace(float **X,int m)  //释放数组空间
{
    int i;
    for(i=0;i<m;i++)
        delete []X[i];
    delete []X;
}

void init(float **X,int m,int n)
{
    for (i=0;i<m;i++)
        for (j=0;j<n;j++)
            X[i][j]=0;
}

int maxof(float *X,int m)
{
    int sub = 0;
    float temp=X[0];
    for (i=1;i<m;i++)
    {
        if (X[i]>temp)
            temp = X[i];
    }
    for(i=0;i<m;i++)
        if(X[i]==temp)
            sub = i;
    return sub;
}

int minof(float *X,int m)//获得非负θ值的最小值
{
    int sub = 0;
    float temp;
    if (X[maxof(X,m)]>0)
    {
        temp = X[maxof(X,m)];
    }
         for (i=0;i<m;i++)
        {
            if (X[i]>0 && X[i]<temp)
                temp = X[i];
        }
        for(i=0; i<m; i++)
            if(X[i] == temp)
                sub = i;
    return sub;
}

bool isAllNegative(float *X,int m)
{
    int flag =0;
    for (i=0;i<m;i++)
        if(X[i] > 0)
            flag += 1;
    if (flag == 0)
        return true;
    else
        return false;
}

//求出基变量的下标（利用与单位矩阵匹配）
void match()
{
    //判断现有矩阵是否包含单位矩阵，从而得到基变量
    for (i=0;i<row;i++)//初始化
        R[i]=0;

    int k=0,flag=0;
    for(k=0;k<row;k++)
    {
        for (i=0;i<col;i++)
        {
            for(j=0;j<row;j++)
            {
                if (A[j][i] != e[j][k])
                    flag +=1;
            }
            if (flag == 0)
            {
                R[k] = i+1;
                break;
            }
            flag = 0;
        }
        flag = 0;
    }

    //若不存在单位向量，则构造缺少向量，并为系数向量添上-M
    i = 0;
    while(i<row)
    {
        if (R[i]==0)
        {
            for (j=0;j<row;j++)
            {
                A[j][col+counter] = e[j][counter];
                C[col+counter] = -M;
                R[i] = col+counter+1;
            }//for
            counter++;
        }//if
        i++;
    }//while
}

//判断是否无解
bool isNoSovle()
{
    int flag = 0;
    for(i=0;i<row;i++)
    {
        if(R[i]>col)
            flag += 1;
    }
    if (flag == 0)
        return false;
    else
        return true;
}

bool isNoLimitSolve()//判断是否为无界解
{
    if (isAllNegative(test,col+counter))
        return 0;
    else
    {
        int flag=0,bol=0;
        for (i=0;i<col+counter;i++)
        {
            flag=0;
            if (test[i]>0)
                for (j=0;j<row;j++)
                    if (A[j][i]>0)//A[j][i]全部<0，则为无界解
                        flag += 1;
            if (flag != 0)
                bol += 1;
        }

        if (bol == 0)
            return 1;
        else
            return 0;
    }
}

bool isMuchSolves()//判断是否有无穷多解
{
    int flag=0;
    for (i=0;i<col+counter;i++)
    {
        for (j=0;j<row;j++)
            if (i+1 == R[j])//R大于1至(col+counter)
                flag = 1;
        if (flag == 0&&test[i] == 0)
            return 1;
        else
            flag = 0;
    }
    return 0;
}

void output()
{
    cout<<"-------------------------------------------------------------------------------\n";
    cout<<"                ";
    for (i=0;i<col+counter;i++)
        cout<<C[i]<<'\t';
    cout<<endl;
    for(i=0;i<row;i++)
    {
        cout<<CB[i]<<"   "<<b[i]<<"   "<<'\t';
        for(j=0;j<col+counter;j++)
            cout<<A[i][j]<<'\t';
        cout<<endl;
    }
    cout<<"                ";
    for(i=0;i<col+counter;i++)
        cout<<test[i]<<'\t';
    cout<<"\n-------------------------------------------------------------------------------\n";

}