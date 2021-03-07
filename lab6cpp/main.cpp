#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

using namespace std;

double potential[101][101], charge[101][101] = {0};
double const h = 0.5; //cell size
double const k = M_PI * (100*h); //let non-zero integer be 1

void dirichletBC(double [][101]);
void chargeInitialization(double [][101]);
void jacobiT1(double[][101], double[][101]);
void jacobiT2(double [][101], double [][101]); //this one accounts for charge
double gausslaw(double [][101]);

int main()
{
    //task 1
    ofstream numplot; //numerical plot
    numplot.open("theoretical task1.csv");
    ofstream analyticalplot; //analytical plot
    analyticalplot.open("analytical task1.csv");

    //numerical plot calculation
    dirichletBC(potential);
    //initialization of potential arrays
    for(int i=1; i<100; i++)
    {
        for(int j=1; j<100; j++)
        {
            potential[i][j] = 1;
        }
    }
    //initialization of charge arrays
    chargeInitialization(charge);
    jacobiT1(potential, charge);
    for(int i=0; i<101; i++)
    {
        for(int j=0; j<101; j++)
        {
            numplot<<potential[i][j]<<",";
        }
        numplot<<endl;
    }

    //analytical plot of task 1
    for(int m=0; m<101; m++)
    {
        for(int n=0; n<101; n++)
        {
            potential[m][n] = exp(-k*n)*sin(k*m);
            analyticalplot << potential[m][n] << ",";
        }
        analyticalplot<<endl;
    }


    //task 2
    ofstream result;
    result.open("gauss law verification.txt");
    //initializing potential arrays
    for(int i=0; i<101; i++)
    {
        for(int j=0; j<101; j++)
        {
            potential[i][j] = 0;
        }
    }
    chargeInitialization(charge);
    //jacobi iteration for task 2
    jacobiT2(potential, charge);
    result<<"The computational result is "<<gausslaw(potential);
    cout<<gausslaw(potential)<<endl;

    return 0;
}

double gausslaw(double potential[][101])
{
    double sum =0; //total sum of e-field
    double sum1=0, sum2=0, sum3=0, sum4=0;
    double fieldx[100][100], fieldy[100][100];
    for(int i=1; i<100; i++)
    {
        for(int j=1; j<100; j++)
        {
            fieldx[i][j] = -(potential[i+1][j] - potential[i-1][j])/(2*h);
            fieldy[i][j] = -(potential[i][j+1] - potential[i][j-1])/(2*h);
        }
    }

    //boundary is the square area from (30,30) to (70,70)
    for(int i=31; i<70; i++)
    {
        sum1 += fieldx[70][i] * h;
    }
    for(int j=30; j<69; j++)
    {
        sum2 -= fieldx[30][j] * h;
    }
    for(int l=31; l<70; l++)
    {
        sum3 += fieldy[l][70] * h;
    }
    for(int m=30; m<69; m++)
    {
        sum4 -= fieldy[m][30] * h;
    }

    return sum = sum1+sum2+sum3+sum4;
}

void jacobiT2(double potential[][101], double charge[][101])
{
    double dv = 1;
    double ptmp[101][101];

    int trackcount = 0; //for counting the iteration loops
    do
    {
        dv = 0;
        for(int i=1; i<100; i++)
        {
            for(int j=1; j<100; j++)
            {
                ptmp[i][j] = (potential[i+1][j] + potential[i-1][j] + potential[i][j+1] + potential[i][j-1] +h*h*charge[i][j])/4;

                if( fabs(ptmp[i][j] - potential[i][j]) > dv)
                {
                    dv = fabs(ptmp[i][j] - potential [i][j]);
                }
            }
        }
        for(int a=1; a<100; a++)
        {
            for(int b=1; b<100; b++)
            {
                potential[a][b] = ptmp[a][b];
            }
        }
        trackcount++;
    }while((trackcount<100000) && dv >1e-20);
    cout<<trackcount<<endl;
}

void jacobiT1(double potential[][101], double charge[][101])
{
    double dv=1; //diff in potential, set to 1 to ensure loop is ran
    double ptmp[101][101]; //temp potential array

    int trackcount = 0; //for counting iterations

    do
    {
        for(int i=1; i<100; i++)
        {
            for(int j=1; j<100; j++)
            {
                ptmp[i][j] = (potential[i+1][j] + potential[i-1][j] + potential[i][j+1] + potential[i][j-1])/4;

                if( fabs(ptmp[i][j] - potential[i][j]) > dv)
                {
                    dv = fabs(ptmp[i][j] - potential[i][j]);
                }
            }
        }

        for(int a=1; a<100; a++)
        {
            for(int b=1; b<100; b++)
            {
                potential[a][b] = ptmp[a][b];
            }
        }
        trackcount++;
    }while((trackcount<10000) && dv > 1e-10);
}

void chargeInitialization(double charge[][101])
{
    charge[50][50] = charge[50][51] = charge[51][50] = charge[51][51] = 0.01;
}

void dirichletBC(double potential[][101])
{
     //using the particular solution
    for(int i=0; i<101; i++)
    {
        potential[i][0] = 0;
        potential[i][100] = exp(-k*i)*sin(k*100);
    }
    for(int j=0; j<101; j++)
    {
        potential[0][j] = sin(k*j);
        potential[100][j] = exp(-k*100)*sin(k*j);
    }

    return;
}
