#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef PI
#define PI    3.141592653589793238463
#endif

double u1_1(double t, double x1, double x2, double h1, double tv)
{
    if (x1<h1/2) return tv;
    else return 0;
}
double u2_1(double t, double x1, double x2)
{
return 0;
}
double g_1(double t, double x1, double x2, double h1, double tg)
{
    if (x1<h1/2) return tg;
    else return 0.0;
}
double theta_1(double t, double x1, double x2, double h1, double tt)
{
    if (x1<h1/2) return tt;
    else return 293;
}
double f0(double t, double x1, double x2)
{
return 0;
}
// tochnoe reshenie

double u1(double t, double x1, double x2)
{
    double res=sin(2*PI*x1)*sin(2*PI*x2)*exp(t);
    return res;
}
double u2(double t, double x1, double x2)
{
    double res=sin(2*PI*x1)*sin(2*PI*x2)*exp(-t);
    return res;
}
double g(double t, double x1, double x2)
{
    double res=log((cos(2*PI*x1)+1.5)*(sin(2*PI*x2)+1.5)*exp(t));
    return res;
}
double theta(double t, double x1, double x2)
{
    double res=(cos(3*PI*x1)+1.5)*(sin(3*PI*x2)+1.5)*exp(t);
    return res;
}
// pravie chasti
double fu1(double t, double x1, double x2, double R, double mu, double koef)
{
    double res1 = sin(2*PI*x1)*sin(2*PI*x2)*exp(t); // du1/dt = u1
    double res2 = sin(2*PI*x1)*sin(2*PI*x2)*exp(2*t) * 2*PI*cos(2*PI*x1)*sin(2*PI*x2) + sin(2*PI*x1)*sin(2*PI*x2) * 2*PI*sin(2*PI*x1)*cos(2*PI*x2);
    double res3 = R*(cos(3*PI*x1)+1.5)*(sin(3*PI*x2)+1.5)*exp(t)*(-2*PI*sin(2*PI*x1))/(cos(2*PI*x1)+1.5) - 3*PI*R*sin(3*PI*x1)*(sin(3*PI*x2)+1.5)*exp(t);
    double res4 = koef*mu*4*PI*PI*( (-4./3.)*sin(2*PI*x1)*sin(2*PI*x2)*exp(t) - sin(2*PI*x1)*sin(2*PI*x2)*exp(t) + (1./3.)*cos(2*PI*x1)*cos(2*PI*x2)*exp(-t) ) / ( (cos(2*PI*x1)+1.5)*(sin(2*PI*x2)+1.5)*exp(t) );
    double res5 = 0;//u1(t,x1,x2);
return res1+res2+res3-res4 + res5;
}
double fu2(double t, double x1, double x2, double R, double mu, double koef)
{
    double res1 = -sin(2*PI*x1)*sin(2*PI*x2)*exp(-t); // du2/dt = -u2
    double res2 = sin(2*PI*x1)*sin(2*PI*x2)*exp(-2*t) * 2*PI*sin(2*PI*x1)*cos(2*PI*x2) + sin(2*PI*x1)*sin(2*PI*x2) * 2*PI*cos(2*PI*x1)*sin(2*PI*x2);
    double res3 = R*(cos(3*PI*x1)+1.5)*(sin(3*PI*x2)+1.5)*exp(t)*(2*PI*cos(2*PI*x2))/(sin(2*PI*x2)+1.5) + 3*PI*R*cos(3*PI*x2)*(cos(3*PI*x1)+1.5)*exp(t);
    double res4 = koef*mu*4*PI*PI*( (-4./3.)*sin(2*PI*x1)*sin(2*PI*x2)*exp(-t) - sin(2*PI*x1)*sin(2*PI*x2)*exp(-t) + (1./3.)*cos(2*PI*x1)*cos(2*PI*x2)*exp(t) ) / ( (cos(2*PI*x1)+1.5)*(sin(2*PI*x2)+1.5)*exp(t) );
    double res5 = 0;//u2(t,x1,x2);
return res1+res2+res3-res4 + res5;
}
double fg(double t, double x1, double x2)
{
    double res1 = exp(t)*sin(2*PI*x1)*sin(2*PI*x2)*(-2*PI*sin(2*PI*x1))/(cos(2*PI*x1)+1.5) + exp(-t)*sin(2*PI*x1)*sin(2*PI*x2)*(2*PI*cos(2*PI*x2))/(sin(2*PI*x2)+1.5);
    double res2 = 2*PI*exp(t)*cos(2*PI*x1)*sin(2*PI*x2) + 2*PI*exp(-t)*sin(2*PI*x1)*cos(2*PI*x2);
return 1. + res1 + res2;
}
double ft(double t, double x1, double x2, double R, double mu, double cv, double kappa)
{
    double res1 = cv*( (cos(3*PI*x1)+1.5)*(sin(3*PI*x2)+1.5)*exp(t) - 3*PI*sin(3*PI*x1)*(sin(3*PI*x2)+1.5)*exp(2*t)*sin(2*PI*x1)*sin(2*PI*x2) + 3*PI*(cos(3*PI*x1)+1.5)*cos(3*PI*x2)*sin(2*PI*x1)*sin(2*PI*x2) );
    double res2 = (cos(2*PI*x1)+1.5)*(sin(2*PI*x2)+1.5)*exp(t);
    double res3 = -9*PI*PI*kappa*(cos(3*PI*x1)*(sin(3*PI*x2)+1.5) +(cos(3*PI*x1)+1.5)*sin(3*PI*x2))*exp(t);
    double res4 = -(2/3)*mu*4*PI*PI*(cos(2*PI*x1)*sin(2*PI*x2)*cos(2*PI*x1)*sin(2*PI*x2)*exp(2*t) + sin(2*PI*x1)*cos(2*PI*x2)*sin(2*PI*x1)*cos(2*PI*x2)*exp(-2*t) + 2*cos(2*PI*x1)*cos(2*PI*x2)*sin(2*PI*x1)*sin(2*PI*x2));
    double res5 = 2*mu*4*PI*PI*( exp(2*t)*cos(2*PI*x1)*sin(2*PI*x2)*cos(2*PI*x1)*sin(2*PI*x2) + exp(-2*t)*sin(2*PI*x1)*cos(2*PI*x2)*sin(2*PI*x1)*cos(2*PI*x2) + 0.5*(cos(2*PI*x1)*sin(2*PI*x2)*exp(t)+sin(2*PI*x1)*cos(2*PI*x2)*exp(-t))*(cos(2*PI*x1)*sin(2*PI*x2)*exp(t)+sin(2*PI*x1)*cos(2*PI*x2)*exp(-t)) );
    double res6 = -R*theta(t,x1,x2)*2*PI*(cos(2*PI*x1)*sin(2*PI*x2)*exp(t) + sin(2*PI*x1)*cos(2*PI*x2)*exp(-t));
return res1 - res3/res2 - res4/res2 - res5/res2 - res6;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double theta2(double t, double x1, double x2)
{
    double res=(cos(4*PI*x1)+1.5)*(sin(4*PI*x2)+1.5)*exp(t);
    return res;
}
double ft2(double t, double x1, double x2, double R, double mu, double cv, double kappa)
{
    double res1 = cv*( (cos(4*PI*x1)+1.5)*(sin(4*PI*x2)+1.5)*exp(t) - 4*PI*sin(4*PI*x1)*(sin(4*PI*x2)+1.5)*exp(2*t)*sin(2*PI*x1)*sin(2*PI*x2) + 4*PI*(cos(4*PI*x1)+1.5)*cos(4*PI*x2)*sin(2*PI*x1)*sin(2*PI*x2) );
    double res2 = (cos(2*PI*x1)+1.5)*(sin(2*PI*x2)+1.5)*exp(t);
    double res3 = -16*PI*PI*kappa*(cos(4*PI*x1)*(sin(4*PI*x2)+1.5) +(cos(4*PI*x1)+1.5)*sin(4*PI*x2))*exp(t);
    double res4 = -(2/3)*mu*4*PI*PI*(cos(2*PI*x1)*sin(2*PI*x2)*cos(2*PI*x1)*sin(2*PI*x2)*exp(2*t) + sin(2*PI*x1)*cos(2*PI*x2)*sin(2*PI*x1)*cos(2*PI*x2)*exp(-2*t) + 2*cos(2*PI*x1)*cos(2*PI*x2)*sin(2*PI*x1)*sin(2*PI*x2));
    double res5 = 2*mu*4*PI*PI*( exp(2*t)*cos(2*PI*x1)*sin(2*PI*x2)*cos(2*PI*x1)*sin(2*PI*x2) + exp(-2*t)*sin(2*PI*x1)*cos(2*PI*x2)*sin(2*PI*x1)*cos(2*PI*x2) + 0.5*(cos(2*PI*x1)*sin(2*PI*x2)*exp(t)+sin(2*PI*x1)*cos(2*PI*x2)*exp(-t))*(cos(2*PI*x1)*sin(2*PI*x2)*exp(t)+sin(2*PI*x1)*cos(2*PI*x2)*exp(-t)) );
    double res6 = -R*theta2(t,x1,x2)*2*PI*(cos(2*PI*x1)*sin(2*PI*x2)*exp(t) + sin(2*PI*x1)*cos(2*PI*x2)*exp(-t));
return res1 - res3/res2 - res4/res2 - res5/res2 - res6;
}
double fu12(double t, double x1, double x2, double R, double mu, double koef)
{
    double res1 = sin(2*PI*x1)*sin(2*PI*x2)*exp(t); // du1/dt = u1
    double res2 = sin(2*PI*x1)*sin(2*PI*x2)*exp(2*t) * 2*PI*cos(2*PI*x1)*sin(2*PI*x2) + sin(2*PI*x1)*sin(2*PI*x2) * 2*PI*sin(2*PI*x1)*cos(2*PI*x2);
    double res3 = 10*(-2*PI*sin(2*PI*x1))/(cos(2*PI*x1)+1.5);
    double res4 = mu*4*PI*PI*( (-4./3.)*sin(2*PI*x1)*sin(2*PI*x2)*exp(t) - sin(2*PI*x1)*sin(2*PI*x2)*exp(t) + (1./3.)*cos(2*PI*x1)*cos(2*PI*x2)*exp(-t) ) / ( (cos(2*PI*x1)+1.5)*(sin(2*PI*x2)+1.5)*exp(t) );
return res1+res2+res3-res4;
}
double fu22(double t, double x1, double x2, double R, double mu, double koef)
{
    double res1 = -sin(2*PI*x1)*sin(2*PI*x2)*exp(-t); // du2/dt = -u2
    double res2 = sin(2*PI*x1)*sin(2*PI*x2)*exp(-2*t) * 2*PI*sin(2*PI*x1)*cos(2*PI*x2) + sin(2*PI*x1)*sin(2*PI*x2) * 2*PI*cos(2*PI*x1)*sin(2*PI*x2);
    double res3 = 10*(2*PI*cos(2*PI*x2))/(sin(2*PI*x2)+1.5);
    double res4 = mu*4*PI*PI*( (-4./3.)*sin(2*PI*x1)*sin(2*PI*x2)*exp(-t) - sin(2*PI*x1)*sin(2*PI*x2)*exp(-t) + (1./3.)*cos(2*PI*x1)*cos(2*PI*x2)*exp(t) ) / ( (cos(2*PI*x1)+1.5)*(sin(2*PI*x2)+1.5)*exp(t) );
return res1+res2+res3-res4;
}
double fg2(double t, double x1, double x2)
{
    double res0 = exp(g(t,x1,x2));
    double res1 = exp(t)*sin(2*PI*x1)*sin(2*PI*x2)*(-2*PI*sin(2*PI*x1))/(cos(2*PI*x1)+1.5) + exp(-t)*sin(2*PI*x1)*sin(2*PI*x2)*(2*PI*cos(2*PI*x2))/(sin(2*PI*x2)+1.5);
    double res2 = 2*PI*exp(t)*cos(2*PI*x1)*sin(2*PI*x2) + 2*PI*exp(-t)*sin(2*PI*x1)*cos(2*PI*x2);

return res0 + exp(g(t,x1,x2))*res1 + exp(g(t,x1,x2))*res2;
}

