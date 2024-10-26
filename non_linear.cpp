#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double function(vector<double> &coeff, double x)
{
    double result = 0.0;
    int deg = coeff.size();
    for (int i = 0; i < deg; i++)
    {
        result += coeff[i] * pow(x, deg - 1 - i);
    }
    return result;
}

double bisectionMethod(vector<double> &coeff, double a, double b, double error_accept)
{
    if (function(coeff, a) * function(coeff, b) >= 0)
    {
        cout << "Bisection method won't work with the given range." << endl;
        return -1;
    }

    double c;
    double error = abs(b - a);

    while (error > error_accept)
    {
        c = (a + b) / 2;
        if (function(coeff, c) == 0.0)
            break;

        if (function(coeff, c) * function(coeff, a) < 0)
            b = c;
        else
            a = c;

        error = abs(b - a);
    }

    cout << "Root found at: " << c << " with error: " << error << endl;
    return c;
}

double falsePositionMethod(vector<double> &coeff, double a, double b, double error_accept)
{
    if (function(coeff, a) * function(coeff, b) >= 0)
    {
        cout << "False Position method is invalid in this range." << endl;
        return -1;
    }

    double c_before = 0;
    double c = (a * function(coeff, b) - b * function(coeff, a)) /
               (function(coeff, b) - function(coeff, a));
    double error = abs(c - c_before);

    while (error > error_accept)
    {
        c_before = c;
        c = (a * function(coeff, b) - b * function(coeff, a)) /
            (function(coeff, b) - function(coeff, a));

        if (function(coeff, c) * function(coeff, a) < 0)
            b = c;
        else
            a = c;

        error = abs(b - a);
    }

    cout << "Root found at: " << c << " with error: " << error << endl;
    return c;
}

double secantMethod(vector<double> &coeff, double x0, double x1, int n)
{
    double xi = x1;

    for (int i = 1; i <= n; i++)
    {
        double fx0 = function(coeff, x0);
        double fx1 = function(coeff, x1);

        if (fx1 == fx0)
        {
            cout << "Division by zero error in Secant method." << endl;
            return -1;
        }

        xi = x1 - fx1 * (x1 - x0) / (fx1 - fx0);

        x0 = x1;
        x1 = xi;
    }

    cout << "Root found at: " << xi << " after " << n << " iterations." << endl;
    return xi;
}

double newtonRaphsonMethod(vector<double> &coeff, vector<double> &deriv, double x0, int n)
{
    double x1;

    for (int i = 1; i <= n; i++)
    {
        double df_x = function(deriv, x0);
        if (df_x == 0)
        {
            cout << "Derivative is zero; Newton-Raphson method fails." << endl;
            return -1;
        }

        x1 = x0 - function(coeff, x0) / df_x;
        x0 = x1;
    }

    cout << "Root found at: " << x1 << " after " << n << " iterations." << endl;
    return x1;
}

int main()
{
    int deg;
    cout << "Enter the degree of the polynomial: ";
    cin >> deg;

    vector<double> coeff(deg + 1);
    cout << "Enter the coefficients from highest to lowest degree: ";
    for (int i = 0; i <= deg; i++)
    {
        cin >> coeff[i];
    }

    double a, b, error_accept;
    cout << "Enter the interval [a, b] for Bisection and False Position methods: ";
    cin >> a >> b;
    cout << "Enter the acceptable error: ";
    cin >> error_accept;

    bisectionMethod(coeff, a, b, error_accept);
    falsePositionMethod(coeff, a, b, error_accept);

    double x0, x1;
    int n;
    cout << "Enter initial guesses for Secant method (x0 and x1): ";
    cin >> x0 >> x1;
    cout << "Enter the number of iterations for Secant method: ";
    cin >> n;
    secantMethod(coeff, x0, x1, n);

    cout << "Enter initial guess for Newton-Raphson method: ";
    cin >> x0;
    cout << "Enter the number of iterations for Newton-Raphson method: ";
    cin >> n;

    vector<double> deriv(deg);
    for (int i = 0; i < deg; i++)
    {
        deriv[i] = coeff[i] * (deg - i);
    }

    newtonRaphsonMethod(coeff, deriv, x0, n);

    return 0;
}
