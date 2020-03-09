#include <stdlib.h>
#include <stdio.h>
#include <math.h>

const double dx_taylor = 0.01;


double HE_function(double x)
{
  // calculate H(x) = 1/2 . (x^2  +  arctan(sqrt(1/x^2 - 1)) / sqrt(1/x^2 - 1))

  if(fabs(x - 1.) <= dx_taylor)             // taxlor expansion around x = 1
  {
    return 0.0005624737680297483 + x*(0.7797825941381623 + x*(0.025694875450339216 +
         x*(0.32086022134934217 + x*(-0.19450407011549792 +
               x*(0.09353199402114554 + x*(-0.032113982358561236 + (0.006866621222216209 - 0.0006807274751769723*x)*x))))));
  }
  else if(x > 0. && x < 1. - dx_taylor)     // 0 < x < 1
  {
    double x2 = x * x;
    double arg = sqrt(1./x2  -  1.);
    return (x2  +  atan(arg) / arg) / 2.;
  }
  else if(x > 1. + dx_taylor)               // x > 1
  {
    double x2 = x * x;
    double arg = sqrt(1. - 1./x2);
    return (x2  +  atanh(arg) / arg) / 2.;
  }
  else
  {
    printf("H_function error: x = %lf is out of bounds", x);
    exit(-1);
  }
}


double HT_function(double x)
{
  // calculate H_T(x) = (x^2  +  (1 - 2x^2).arctan(sqrt(1/x^2 - 1)) / sqrt(1/x^2 - 1)) / (1 - x^2)

  if(fabs(x - 1.) <= dx_taylor)             // taxlor expansion around x = 1
  {
    return -0.003202301654317485 + x*(1.6035437558657903 + x*(-0.15518846478629494 +
         x*(-0.32881655667957355 + x*(0.3837277273799249 + x*
                (-0.2432341652146951 + x*(0.09666948118913454 + (-0.0225216819642642 + 0.002355539197625668*x)*x))))));
  }
  else if(x > 0. && x < 1. - dx_taylor)     // 0 < x < 1
  {
    double x2 = x * x;
    double arg = sqrt(1./x2  -  1.);
    return (x2  +  (1. - 2.*x2) * atan(arg) / arg) / (1. - x2);
  }
  else if(x > 1. + dx_taylor)               // x > 1
  {
    double x2 = x * x;
    double arg = sqrt(1. - 1./x2);
    return (x2  +  (1. - 2.*x2) * atanh(arg) / arg) / (1. - x2);
  }
  else
  {
    printf("HT_function error: x = %lf is out of bounds", x);
    exit(-1);
  }
}


double HL_function(double x)
{
  // calculate H_L(x) = x^2 . (-x^2  +  arctan(sqrt(1/x^2 - 1)) / sqrt(1/x^2 - 1)) / (1 - x^2)

  if(fabs(x - 1.) <= dx_taylor)             // taxlor expansion around x = 1
  {
    return 0.004372373412624753 + x*(-0.04433956136745337 + x*(0.2078416939099339 +
         x*(0.9680100429323358 + x*(-0.7695771720535178 + x*
                (0.42777119681106374 + x*(-0.15963396768329577 + (0.03589393063070768 - 0.0036718699257309944*x)*x))))));
  }
  else if(x > 0. && x < 1. - dx_taylor)     // 0 < x < 1
  {
    double x2 = x * x;
    double arg = sqrt(1./x2  -  1.);
    return (-x2  +  atan(arg) / arg) * x2 / (1. - x2);
  }
  else if(x > 1. + dx_taylor)               // x > 1
  {
    double x2 = x * x;
    double arg = sqrt(1. - 1./x2);
    return (-x2  +  atanh(arg) / arg) * x2 / (1. - x2);
  }
  else
  {
    printf("HL_function error: x = %lf is out of bounds", x);
    exit(-1);
  }
}

