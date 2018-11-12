/*
 * qpMPCUpdate.cpp
 * Copyright (C) 2018 Gao Tang <gt70@duke.edu>
 *
 * Distributed under terms of the MIT license.
 */
/*  Use one-step QP to update the predicted trajectory  */
#include "qpMPCUpdate.h"
#include "qpUpdate.h"
#include "Tools/VecMat.h"
#include "stdio.h"
#include "math.h"
#include "json/json.h"
#include <fstream>
#include "utilityClass.h"
#include "Utility.h"
#include "osqp/osqp.h"
#include "cnpy.h"

TrajOptMPC::dirTranODE *extdT = nullptr;

void funwrapper(const double t, cRefV y, cRefV p, RefV dy, RefM J){
    extdT->RKDynFun(t, y, p, dy, J);
}

int qpMPC::qpMPCUpdater::qpMPCUpdateOSQP(RefV xin, cRefV x0){
    extdT = dT;
    int lenx0 = xin.size();
    VX xout = VX::Zero(lenx0);
    dT->x0 = x0;
    dT->xf.setZero();
    std::cout << "xin " << xin.size() << "xout" << xout.size() << "\n";
    int flag = qpUpdateOSQP(extdT, xin, xout);
    std::cout << "OSQP finish\n";
    xin = xout;  // copy back, so we did an implace update
}

MX qpMPC::acceGen(RefV xin, bool mpc){
    int dimx = 12, dimu = 4;
    int N = (xin.size() + dimu) / (dimx + dimu);
    std::cout << "N =" << N << std::endl;
    Rotor sys;
    MX aout = MX::Zero(3, N);
    VX f = VX::Zero(dimx);
    int tmpint = N * dimx;
    for(int i = 0; i < N; i++){
        if(i < N - 1){
            sys.dyn(0, xin.segment(i*dimx, dimx), xin.segment(tmpint + i*dimu, dimu), f);
            aout.col(i) = f.segment(6, 3);
        }
        else{
            aout.col(i).setZero();
        }
        
    }
    std::cout << "we get " << aout << std::endl;
    return aout;
}