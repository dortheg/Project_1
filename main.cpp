#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <time.h>
#include <armadillo>

using namespace std;
using namespace arma;

double func_f(double x, double h);
double func_u(double x);
void general_solver(double *u, int n);
void special_solver(double *u, int n);
void LU_arma(int n);


//Hovedprogram
int main(int argc, char *argv[])
{
    int n = 10; //atoi(argv[1]);                                  //Tar inn kommandolinjeargument
    clock_t start_1, start_2, start_3, end_1, end_2, end_3;

    double* u1  = new double[n];
    start_1 = clock();
    general_solver(u1, n);
    end_1 = clock();

    double* u2 = new double[n];
    start_2 = clock();
    special_solver(u2,n);
    end_2 = clock();


//    start_3 = clock();
//    LU_arma(n);
//    end_3 = clock();

    //double t_1, t_2, t_3;
    //t_1 = (end_1 - start_1)/((double)CLOCKS_PER_SEC);
    //t_2 = (end_2 - start_2)/((double)CLOCKS_PER_SEC);
    //t_3 = (end_3 - start_3)/((double)CLOCKS_PER_SEC);


    //Skriver kjoretidene til fil
//    ofstream myfile;
//    myfile.open("../time.txt"); // "../time.txt", hvis legge fil i mappe opp
//    myfile << t_1 << "\n";
//    myfile << t_2 << "\n";
//    myfile << t_3 << "\n";
//    myfile.close();

    return 0;
}


double func_f(double x, double h)
{
    double f;

    f = pow(h,2)*100* exp(-10*x);


    return f;
}

double func_u(double x)
{
    double u;

    u = 1 - (1-exp(-10))*x - exp(-10*x);

    return u;
}

void general_solver(double *u, int n)
{
    // Trenger arrayene som definerer matrisen (problemet), og arrayer aa legge losningene i
    //Gjor array dynamiske, fyll automatisk med n elementer
    double* a = new double[n];
    double* b = new double[n-1];
    double* c = new double[n-1];

    //Fyller arrayene med n elementer
    fill_n(a, n, 2);
    fill_n(b, n-1, -1);
    fill_n(c, n-1, -1);

    double* a_ = new double[n];
    double* f_ = new double[n];

    //Skal evaluere for n skritt i intervallet x (0,1), lager en slik x-array
    double* x = new double [n];

    double h; h = 1.0/(n+1); //Har (n+1) intervaller

    //Loper til n
    for(int i=0; i < n; i++)
    {
        x[i] = (i+1)*h;
    }

    //Lager f-arrayen (funksjonsverdien multiplisert med h²)

    double *f = new double[n];

    for(int i=0; i < n; i++)
    {
        f[i] = func_f(x[i], h);
    }


    // Setter initialbetingelser
    a_[0] = a[0];
    f_[0] = f[0];

    // Gjor Forward Stubstitution
    for(int i=1; i<n; i++)
    {
        double c_i_a_i_;
        c_i_a_i_ = c[i-1]/a_[i-1];

        a_[i] = a[i] - b[i-1]*c_i_a_i_;

        f_[i] = f[i] - f_[i-1]*c_i_a_i_;
    }

    //Setter initialbetingelse for u, og gjor Backward Substitution
    u[n-1] = f_[n-1]/a_[n-1];

    for(int i=2; i<=n; i++)
    {
        u[n-i] = (f_[n-i] - b[n-i]*u[n-i+1])/a_[n-i];
    }

    //Limer paa initialbetingelsene u(0) = u(1) = 0
    double* u_new = new double [n+2];

    //Looper til n+1, fyller opp slik at u_new[1:n+1] = u[0:n], kan jeg i stedet skrive dette?
    for(int i=1; i < n+1; i++)
    {
        u_new[i] = u[i-1];
    }

    u_new[0] =0;
    u_new[n+1]= 0;


    ofstream myfile;
    myfile.open("../u_file.txt");

    for(int i=0; i<n+2; i++)
    {
        myfile << u_new[i] << "\n";
    }
    myfile.close();
}

void special_solver(double *u, int n)
{
    double a = 2.0;

    double* a_ = new double[n];
    double* f_ = new double[n];

    //Skal evaluere for n skritt i intervallet x (0,1), lager en slik x-array
    double* x = new double [n];

    double h; h = 1.0/(n+1); //Har (n+1) intervaller

    //Loper til (n-1)
    for(int i=0; i < n; i++)
    {
        x[i] = (i+1)*h;
    }

    //Lager f-arrayen (funksjonsverdien multiplisert med h²)

    double *f = new double[n];

    for(int i=0; i < n; i++)
    {
        f[i] = func_f(x[i], h);
    }


    // Setter initialbetingelser
    a_[0] = a;
    f_[0] = f[0];

    // Gjor Forward Stubstitution
    for(int i=1; i<n; i++)
    {
        a_[i] = (i+2.0)/(i+1.0); //Mismatch med indeksene

        f_[i] = f[i] + f_[i-1]/a_[i-1];;
    }

    //Setter initialbetingelse for u, og gjor Backward Substitution
    u[n-1] = f_[n-1]/a_[n-1];

    for(int i=2; i<=n; i++)
    {
        u[n-i] = (f_[n-i] + u[n-i+1])/a_[n-i];
    }

    double* u_new = new double [n+2];
    //u_new[1:n+1] = u[0:n]
    for(int i=1; i < n+1; i++)
    {
        u_new[i] = u[i-1];
    }

    //Skriver u-arrayen til fil
    ofstream myfile;
    myfile.open("../u_file.txt");
    for(int i=0; i<n+2; i++)
    {
        myfile << u_new[i] << "\n";
    }
    myfile.close();
}

void LU_arma(int n)
{
    //Skal lose problemet Au = f vha Armadillo sin LU-dekomposisjon
    mat A = zeros<mat>(n, n);
    mat f = zeros<mat>(n, 1);

    //Lager matrisen A
    for(int i=0; i<n; i++)
    {
        A(i,i) = 2;

        if(i < (n-1))
        {
            A(i,i+1) = -1;
        }

        if (i < (n-1))
        {
            A(i+1,i) = -1;
        }
    }

    //Lager arrayen f

    double* x = new double [n];

    double h; h = 1.0/(n+1); //Har (n+1) intervaller

    //Loper til (n-1)
    for(int i=0; i < n; i++)
    {
        x[i] = (i+1)*h;
    }

    //Lager f-arrayen (funksjonsverdien multiplisert med h²)

    for(int i=0; i < n; i++)
    {
        f[i] = func_f(x[i], h);
    }


    mat L = zeros<mat>(n,n);
    mat U = zeros<mat>(n,n);
    mat P = zeros<mat>(n,n);

    lu(L, U, P, A);

    //Loser Ly = f
    mat y = solve(L, P*f);

    //Loser Uu = y
    mat u = solve(U, y);

    //Limer paa 0 i endepunktene, slik at u(0) = u(1) = 0
    double* u_new = new double [n+2];

    //u_new[1:n+1] = u[0:n]
    for(int i=1; i < n+1; i++)
    {
        u_new[i] = u[i-1];
    }

    u_new[0] =0;
    u_new[n+1]= 0;

    //Skriver u-arrayen til fil
    ofstream myfile;
    myfile.open("../u_file.txt");
    for(int i=0; i<n+2; i++)
    {
        myfile << u_new[i] << "\n";
    }
    myfile.close();

}
