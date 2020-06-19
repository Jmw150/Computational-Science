//
#include <stdio.h>
#include <stdlib.h>

double function(double x)
{
    return x * x;
}

double steffensenMethod(double(*f)(double), 
    double seed, double precision, int iterations)
{
    if(precision == 0)
    {
        exit(-1);
    }

    double h = precision;
    double x = seed;
    // legal gnu c
    double g(double x)
    {
        return (f(x + h) - f(x)) / h;
    }

    int i = iterations;

    while(i > 0)
    {
        x -= f(x) / ((float)(g(x)));
        i -= 1;
    }
    return x;

}


int main()
{
    printf("%f\n", 
        steffensenMethod(&function, 3, 0.1, 50));
    return 0;
}
