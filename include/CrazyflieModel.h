#ifndef CRAZYFLIE_MODEL
#define CRAZYFLIE_MODEL

#include <math.h>


/**
 Input: PWM, range from 0-65535(0%-100%)
 Ouput: thrust of the whole crazyflie drone, in Newton(N)
*/
double PWM2Thrust(double PWM)
{
    return 4*(1.828091e-11*PWM*PWM + 1.187259e-06*PWM -0.000936047);
}

/**
 Input: Thrust needed for the whole drone, in Newton(N)
 Output: PWM, range from 0-65535(0%-100%)
*/
double Thrust2PWM(double thrust)
{
  thrust = thrust/4;
  double a = 1.828091e-11;
  double b = 1.187259e-06;
  double c = -0.000936047;

  c = c - thrust;

  double delta = sqrt(b*b - 4*a*c);
  return (-b+delta)/2/a;
}

#endif /* end of include guard: CRAZYFLIE_MODEL */
