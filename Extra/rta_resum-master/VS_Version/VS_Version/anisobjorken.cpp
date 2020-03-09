#include <math.h>
#include <cmath>

const double delta = 0.01;

double compute_zetaLL(double e, double pl, double prefactor)
{
  double x   = pl / e;
  double x2  = x   * x;
  double x3  = x2  * x;
  double x4  = x3  * x;
  double x5  = x4  * x;
  double x6  = x5  * x;
  double x7  = x6  * x;
  double x8  = x7  * x;
  double x9  = x8  * x;
  double x10 = x9  * x;
  double x11 = x10 * x;
  double x12 = x11 * x;
  double x13 = x12 * x;
  double x14 = x13 * x;
  double x15 = x14 * x;
  double x16 = x15 * x;
  double x17 = x16 * x;
  double x18 = x17 * x;
  double x19 = x18 * x;
  double x20 = x19 * x;             // x = pl / e
  double x21 = x20 * x;             // z = xi (from aniso hydro)
  double x22 = x21 * x;             // aL = 1/sqrt(1+xi)

  double aL, aL2, z;

  if(x > 1./3.)
  {
    z = (33920.75424130814 - 72210.38155086128*x - 58150.373221605056*x2 - 123258.68865968155*x3 + 14269.18945991164*x4 + 170364.85584208343*x5 +
      169910.50665817957*x6 + 64799.56668758523*x7 + 262850.50962558796*x8 + 35449.81323782106*x9 - 248808.80620651352*x10 -
      288836.0950617432*x11 - 55525.59817904083*x12 - 249947.48251438234*x13 - 310420.8253593438*x14 - 48317.55700989516*x15 +
      522691.7302032236*x16 + 527504.8150488662*x17 + 219759.88337782127*x18 - 187603.57642353655*x19 - 199506.45061878706*x20 -
      611077.848257917*x21 + 432142.0012199023*x22)/
      (-192.35667843131404 + 43491.082537769165*x + 83354.47448892899*x2 + 45103.07343085356*x3 - 105414.36804542418*x4 + 140186.71296754244*x5 -
      531082.9994828509*x6 - 85658.91194364589*x7 + 377783.60198413196*x8 - 339045.0410056553*x9 - 95837.02795785779*x10 +
      284537.5663725089*x11 + 703062.1998023012*x12 + 223019.9316692852*x13 - 501784.5491947427*x14 - 145230.1534789184*x15 +
      55948.62853147295*x16 + 49679.34805386173*x17 - 641771.3022609851*x18 - 274804.41532698454*x19 + 726388.8998660464*x20 +
      350014.57800287893*x21 - 361748.9148710701*x22);

    if(!(z >= -0.99999999 && z <= 1.e-8))
    {
      z = fmax(-0.99999999, fmin(z, 1.e-8));
    }
    aL  = 1. / sqrt(1. + z);
    aL2 = 1. / (1. + z);
  }
  else
  {
    aL = (2.372796737893896e-62 + 9.355496760751141e-54*x + 3.15985529218801e-46*x2 + 2.29804656071578e-39*x3 + 4.8654069671748624e-33*x4 +
      3.4686835009134695e-27*x5 + 9.052410236842743e-22*x6 + 9.132309729581051e-17*x7 + 3.705485165853083e-12*x8 + 6.240802836268058e-8*x9 +
      0.00044799689605487286*x10 + 1.4025011569370325*x11 + 1953.0793537979494*x12 + 1.229812256787706e6*x13 + 3.543561225712354e8*x14 +
      4.697330865356272e10*x15 + 2.8499566740003765e12*x16 + 7.731610782177606e13*x17 + 8.841791912264315e14*x18 + 3.673281425421166e15*x19 +
      3.042059896930142e15*x20 - 3.4368817938638095e15*x21 - 8.169907788507815e14*x22)/
      (7.281820681114894e-58 + 6.793723008169782e-50*x + 1.0073238134263982e-42*x2 + 3.8567133904345664e-36*x3 + 4.6749055427591935e-30*x4 +
      1.9992235460663164e-24*x5 + 3.2233724058457452e-19*x6 + 2.051207320369606e-14*x7 + 5.334552198382988e-10*x8 + 5.833728132253219e-6*x9 +
      0.02748383972843651*x10 + 56.954350298361284*x11 + 52824.406590310646*x12 + 2.2217655338084057e7*x13 + 4.267549397728813e9*x14 +
      3.733806109621652e11*x15 + 1.459513002063948e13*x16 + 2.4180382199020853e14*x17 + 1.4786509784350255e15*x18 + 1.8740406611426415e15*x19 -
      3.345323820802959e15*x20 - 1.2075997985771218e15*x21 + 1.136213305508547e15*x22);

    if(!(aL >= 0.0001 && aL <= 1.0001))
    {
      aL = fmax(0.0001, fmin(aL, 1.0001));
    }

    aL2 = aL * aL;
    z = 1. / aL2  -  1.;
  }

  double z2 = z  * z;
  double z3 = z2 * z;

  double t_200;
  double t_240;

  if(z > delta)
  {
    double sqrtz = sqrt(z);
    double t = atan(sqrtz) / sqrtz;

    t_200 = 1.  +  (1. + z) * t;
    t_240 = (3.  +  2. * z  -  3. * (1. + z) * t) / z2;
  }
  else if(z < -delta && z > -1.)
  {
    double sqrtmz = sqrt(-z);
    double t = atanh(sqrtmz) / sqrtmz;

    t_200 = 1.  +  (1. + z) * t;
    t_240 = (3.  +  2. * z  -  3. * (1. + z) * t) / z2;
  }
  else if(fabs(z) <= delta)
  {
    double z4 = z3 * z;
    double z5 = z4 * z;
    double z6 = z5 * z;

    t_200 = 2. + 0.6666666666666667*z - 0.1333333333333333*z2 + 0.05714285714285716*z3 - 0.031746031746031744*z4 + 0.020202020202020193*z5 - 0.013986013986013984*z6;

    t_240 = 0.4 - 0.17142857142857149*z + 0.09523809523809523*z2 - 0.06060606060606058*z3 + 0.04195804195804195*z4 - 0.030769230769230785*z5 + 0.023529411764705882*z6;
  }

  double Lambda4 = 2. * e / (aL2 * prefactor * t_200);        // Lambda is effective temperature from ahydro

  double I_240 = prefactor * Lambda4 / 2. * t_240 * aL2;      // see appendix in arxiv.1803.01810

  return I_240  -  3. * pl;
}


double de_dt(double e, double pl, double t, double prefactor, double etas)
{
  return - (e + pl) / t;
}


double dpl_dt(double e, double pl, double t, double prefactor, double etas)
{
  double p = e/3.;
  double T = pow(e / prefactor, 0.25);
  double taupiInv = T / (5. * etas);
  double zetaLL = compute_zetaLL(e, pl, prefactor);

  return - taupiInv * (pl - p)  +  zetaLL / t;
}


void run_aniso_bjorken(double *Temp, double pl0, double t0, double dt, int Nt, double prefactor, double etas)
{
  double t  = t0;
  double T0 = Temp[0];

  double e = prefactor * pow(T0, 4);
  double pl = pl0;

  for(int i = 1; i < Nt; i++)
  {
    double e1  = dt *  de_dt(e, pl, t, prefactor, etas);
    double pl1 = dt * dpl_dt(e, pl, t, prefactor, etas);

    double e2  = dt *  de_dt(e + e1/2., pl + pl1/2., t + dt/2., prefactor, etas);
    double pl2 = dt * dpl_dt(e + e1/2., pl + pl1/2., t + dt/2., prefactor, etas);

    double e3  = dt *  de_dt(e + e2/2., pl + pl2/2., t + dt/2., prefactor, etas);
    double pl3 = dt * dpl_dt(e + e2/2., pl + pl2/2., t + dt/2., prefactor, etas);

    double e4  = dt *  de_dt(e + e3, pl + pl3, t + dt, prefactor, etas);
    double pl4 = dt * dpl_dt(e + e3, pl + pl3, t + dt, prefactor, etas);

    e  += (e1   +  2. * e2   +  2. * e3   +  e4) / 6.;
    pl += (pl1  +  2. * pl2  +  2. * pl3  +  pl4) / 6.;

    t += dt;

    Temp[i] = pow(e / prefactor, 0.25);
  }
}


