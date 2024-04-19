#include <iostream>
#include <cmath>
#include <functional>
#include <random>
#include <string>
#include <iomanip>

const double EPSILON = 1e-12; //For the sin function

//a lower bound, b upper bound, N # slabs
//this is trapezoidal rule
double Do_Integral(const double &a, const double &b, const int &N, std::function<double(double)> integrand)
{
    double sum{0.0}, dx{(double)(b-a)/(double)N};
    for(int i=0;i<N;i++)
    {
        sum += integrand(a+(double)i*dx) + integrand(a+((double)i+1.)*dx);
    }
    return (dx/2.) * sum;
}

//Monte Carlo, drawing random samples in our range
double Do_MC(const double &a, const double &b, const int &N, std::function<double(double)> integrand)
{
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distribution(a,b + EPSILON); //Let's draw random numbers b/n a and b (+ EPSILON to include b)
    double sum{0.0};
    for(int i=0;i<N;i++)
    {
        sum += integrand(distribution(generator));
    }
    return (b-a)/((double)N) * sum;
}

//Gaussian Function
double gaus(double x)
{
    return std::exp(-x*x);
}

//Nasty Sinusoid
double sin_func(double x)
{
    double arg{x*(2.-x)+EPSILON};
    double s{std::sin(std::pow(arg,-1))};
    return s*s;
}

//Printing results nicely in tabular form
void print_table(const double &a, const double &b, std::function<double(double)> integrand, const std::string &name)
{
    std::cout << "Integrating " << name << " for a = " << a << " b = " << b << std::endl;
    std::cout << std::setw(12) << "N , " << std::setw(6) << "Trapezoidal , " << std::setw(3) << "Monte Carlo\n";
    int N{8};
    while(N<=1024)
    {
        std::cout << std::setw(12) << "N = " << N << std::setw(3) << " , " << Do_Integral(a,b,N,integrand) << std::setw(3) << " , " << Do_MC(a,b,N,integrand) << std::endl;
        N *= 2;
    }
}

int main() //MC seems to behave just worse for this range
{
    print_table(-5.,5.,gaus,"Gaussian"); //Integrating e^(-x^2)
    std::cout << std::endl;
    print_table(0.,2.,sin_func,"Sinusoid"); //Integrating gross sin^2 term
    return 0;
}