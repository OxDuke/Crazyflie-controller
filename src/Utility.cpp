/* Implementation of utility functions */
#include "Utility.h"
#include <iostream>

#define NO_DEBUG

int RK4(RKFun fun, double t0, double tf, double *y0, double *p, double *y, double *J, int dimx, int dimp){
    MapV vy(y0, dimx), vp(p, dimp), yout(y, dimx);
    MapM Jout(J, dimx, 1 + dimp + dimx);
    return RK4(fun, t0, tf, vy, vp, yout, Jout);
}

int Euler(RK4Fun fun, double t0, double tf, const double *yin, double *yout, int dimx){
    double h = tf - t0;
    cMapV y0(yin, dimx);
    MapV y(yout, dimx);
    VX dy1(dimx);
    fun(t0, y0.data(), dy1.data());
    y = y0 + h*dy1;
    return 1;
}

int RK4(RK4Fun fun, double t0, double tf, const double *yin, double *yout, int dimx){
    double h = tf - t0;
    cMapV y0(yin, dimx);
    MapV y(yout, dimx);
    VX dy1(dimx), dy2(dimx), dy3(dimx), dy4(dimx);
    //Evaluate at the 4 steps
    fun(t0, y0.data(), dy1.data());
    VX yk1 = y0 + 0.5 * h * dy1;
    double tmid = t0 + 0.5 * h;
    fun(tmid, yk1.data(), dy2.data());
    VX yk2 = y0 + 0.5 * h * dy2;
    fun(tmid, yk2.data(), dy3.data());
    VX yk3 = y0 + h * dy3;
    fun(tf, yk3.data(), dy4.data());
    //Calculate final y
    y = y0 + 1.0/6.0*h*(dy1 + 2*dy2 + 2*dy3 + dy4);
    return 1;
}

int RK4(RKFun fun, double t0, double tf, cRefV y0, cRefV p, RefV y, RefM J){
    //Extract dimensional information
    int dimx = y0.size(), dimp = p.size(), colm = 1 + dimx + dimp;
    double h = tf - t0;
    VX dy1(dimx), dy2(dimx), dy3(dimx), dy4(dimx);
    MX J1(dimx, colm), J2(dimx, colm), J3(dimx, colm), J4(dimx, colm), k1X(dimx, colm), k2X(dimx, colm), k3X(dimx, colm), k4X(dimx, colm);
    //Evaluate at the 4 steps
    fun(t0, y0, p, dy1, J1);
#ifdef DEBUG
    std::cout << "0: " << dy1 << J1 << std::endl;
#endif
    VX yk1 = y0 + 0.5 * h * dy1;
    double tmid = t0 + 0.5 * h;
    fun(tmid, yk1, p, dy2, J2);
#ifdef DEBUG
    std::cout << "1: " << dy2 << " " << J2 << std::endl;
#endif
    VX yk2 = y0 + 0.5 * h * dy2;
    fun(tmid, yk2, p, dy3, J3);
#ifdef DEBUG
    std::cout << "2: " << dy3 <<  " " <<  J3 << std::endl;
#endif
    VX yk3 = y0 + h * dy3;
    fun(tf, yk3, p, dy4, J4);
#ifdef DEBUG
    std::cout << "3: " << dy4 << " " <<  J4 << std::endl;
#endif
    //Calculate final y
    y = y0 + 1.0/6.0*h*(dy1 + 2*dy2 + 2*dy3 + dy4);
    //Calculate J
    J.setZero();
    J.middleCols(1, dimx).setIdentity();
    //Calculate k1X to k4X
    k1X.array() = J1.array();
    k1X.col(0).setZero();
#ifdef DEBUG
    std::cout << "k1X " << k1X << std::endl;
#endif

    k2X.setZero();
    k2X.col(0) = 0.5 * J2.col(0);
    k2X.middleCols(1 + dimx, dimp) = J2.middleCols(1 + dimx, dimp);
    k2X.middleCols(1, dimx) = J2.middleCols(1, dimx);
    k2X += 0.5 * h * J2.middleCols(1, dimx) * k1X;
    k2X.col(0) += 0.5 * J2.middleCols(1, dimx) * dy1;
#ifdef DEBUG
    std::cout << "k2X " << k2X << std::endl;
#endif

    k3X.setZero();
    k3X.col(0) = 0.5 * J3.col(0);
    k3X.middleCols(1 + dimx, dimp) = J3.middleCols(1 + dimx, dimp);
    k3X.middleCols(1, dimx) = J3.middleCols(1, dimx);
    k3X += 0.5 * h * J3.middleCols(1, dimx) * k2X;
    k3X.col(0) += 0.5 * J3.middleCols(1, dimx) * dy2;
#ifdef DEBUG
    std::cout << "k3X " << k3X << std::endl;
#endif

    k4X.setZero();
    k4X.col(0) = J4.col(0);
    k4X.middleCols(1 + dimx, dimp) = J4.middleCols(1 + dimx, dimp);
    k4X.middleCols(1, dimx) = J4.middleCols(1, dimx);
    k4X += h * J4.middleCols(1, dimx) * k3X;
    k4X.col(0) += J4.middleCols(1, dimx) * dy3;
#ifdef DEBUG
    std::cout << "k4X " << k4X << std::endl;
#endif
#ifdef DEBUG
    std::cout << "J " << J << std::endl;
#endif
    //Summation to get final result
    J += 1.0/6.0*h*(k1X + 2*k2X + 2*k3X + k4X);
#ifdef DEBUG
    std::cout << "J " << J << std::endl;
#endif
    //for first col
    J.col(0) += 1.0/6.0*(dy1 + 2*dy2 + 2*dy3 + dy4);
#ifdef DEBUG
    std::cout << "J " << J << std::endl;
#endif
#ifdef DEBUG
    std::cout << "final y=" << y << " J = " << J << std::endl;
#endif

    return 1;
}
