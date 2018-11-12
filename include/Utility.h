/* This header file includes several utilities */
#ifndef UTILITY_H
#define UTILITY_H
#include "TigerTools/TigerEigen.h"
#include <functional>

typedef std::function<void(const double, cRefV y, cRefV p, RefV dy, RefM J)> RKFun;
typedef std::function<int(double t, const double *x, double *dx)> RK4Fun;

//RK4 method, which performs one step of Runge-Kutta 4
//p is additional parameter that you might want to differentiate with
//y and J returns the final state and Jacobian w.r.t time, x0, p
int RK4(RKFun fun, double t0, double tf, cRefV y0, cRefV p, RefV y, RefM J);
int RK4(RKFun fun, double t0, double tf, double *y0, double *p, double *y, double *J, int dimx, int dimp);

int Euler(RK4Fun fun, double t0, double tf, const double *yin, double *yout, int dimx);
int RK4(RK4Fun fun, double t0, double tf, const double *yin, double *yout, int dimx);
#endif
