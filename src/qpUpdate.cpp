/*
 * qpUpdate.cpp
 * Copyright (C) 2018 Gao Tang <gt70@duke.edu>
 *
 * Distributed under terms of the MIT license.
 */
/*  Use one-step QP to update the predicted trajectory  */
#include "TrajOpt.h"
#include "Tools/VecMat.h"
#include "stdio.h"
#include "math.h"
#include "json/json.h"
#include <fstream>
#include "utilityClass.h"
#include "Utility.h"
#include "osqp/osqp.h"
#include "cnpy.h"
#include "qpUpdate.h"

dirTranODE *extdT = nullptr;

void funwrapper(const double t, cRefV y, cRefV p, RefV dy, RefM J){
    extdT->RKDynFun(t, y, p, dy, J);
}

int qpUpdater::qpUpdateOSQP(RefV xin, cRefV x0){
    if(extdT == nullptr)
        extdT = dT;
    /* use OSQP to update trajectory, we expand after fixxing xin quickly */
    int N = dT->N, dimx = dT->dimx, dimu = dT->dimu;
    int nSol = dT->numVar;
    int nC = (N - 1) * dimx + (N - 2) * dimx + (N - 1) * dimu;
    int nVar = (N - 1) * dimu + (N - 2) * dimx; // Fix initial and final one right now
    double tf = xin(dT->numVar - 1);
    double h = tf / (N - 1);
    dT->x0 = x0;
    dT->xf.setZero();
    MapM tX(xin.data(), dimx, N);
    MapM tU(xin.data() + N * dimx, dimu, N - 1);

    OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
    settings->verbose = false;
    OSQPWorkspace * work;  // Workspace
    OSQPData * data;  // OSQPData
    // Populate data
    data = (OSQPData *)c_malloc(sizeof(OSQPData));
    data->n = nVar;
    data->m = nC;
    int P_nnz = nVar;
    VX P = VX::Zero(nVar);
    VX q = VX::Zero(nVar);
    typedef Eigen::Matrix<long long, 1, -1> VXl;
    VXl pi = VXl::Zero(nVar);
    VXl pp = VXl::Zero(nVar + 1);
    pp(0) = 0;
    int pind = 0;
    for(int i = 0; i < N - 2; i++){
        for(int j = 0; j < dimx; j++){
            P(pind) = dT->Q(j);
            pi(pind) = pind;
            pp(pind + 1) = pind + 1;
            q(pind) = dT->Q(j) * (tX(j, i + 1) - dT->xbase(j));
            pind += 1;
        }
    }
    for(int i = 0; i < N - 1; i++){
        for(int j = 0; j < dimu; j++){
            P(pind) = dT->R(j);
            pi(pind) = pind;
            pp(pind + 1) = pind + 1;
            q(pind) = dT->R(j) * (tU(j, i) - dT->ubase(j));
            if(dT->enableRDU){
                if(dT->incorrectRDU){
                    P(pind) += dT->RDU(j);
                    if(i == 0){
                        q(pind) += dT->RDU(j) * (tU(j, i) - dT->ubase(j));
                    }
                    else{
                        P(pind - 1) += dT->RDU(j);
                        double tmp = dT->RDU(j) * (tU(j, i) - tU(j, i - 1));
                        q(pind) += tmp;
                        q(pind - 1) -= tmp;
                    }
                }
                else{
                    P(pind) += dT->RDU(j) / h / h;
                    if(i == 0){
                        q(pind) += dT->RDU(j) * (tU(j, i) - dT->ubase(j)) / h / h;
                    }
                    else{
                        P(pind - 1) += dT->RDU(j);
                        double tmp = dT->RDU(j) * (tU(j, i) - tU(j, i - 1)) / h / h;
                        q(pind) += tmp;
                        q(pind - 1) -= tmp;
                    }
                }
            }
            pind += 1;
        }
    }
    // set up A, which is difficult
    int A_nnz = dimx * (dimx + 1) * (N - 2) + (N - 1) * dimx * dimu + (N - 2) *dimx + (N - 1) * dimu;
    VX A = VX::Zero(A_nnz);
    VXl Ai = VXl::Zero(A_nnz);
    VXl Ap = VXl::Zero(nVar + 1);
    VX lbA = VX::Zero(nC);
    VX ubA = VX::Zero(nC);
    int aind = 0;
    int uind = (N - 2) * dimx * (dimx + 1) + (N - 2) * dimx;  //those are already occupied by x
    // assign dynamics
    VX right = VX::Zero(dimx);
    VX y = VX::Zero(dimx);
    MX J = MX::Zero(dimx, 1 + dimx + dimu);
    for(int i = 0; i < N - 1; i++){
        double t0 = i * h, tf = (i + 1) * h;
        if(i == 0) { // special operation for first and last one
            RK4(funwrapper, t0, tf, dT->x0.data(), tU.col(i).data(), y.data(), J.data(), dimx, dimu);
            right = tX.col(i + 1) - y;
        }
        else{
            RK4(funwrapper, t0, tf, tX.col(i).data(), tU.col(i).data(), y.data(), J.data(), dimx, dimu);
            if(i == N - 2)
                right = dT->xf - y;
            else
                right = tX.col(i + 1) - y;
        }
        // if i > 0, we have to assign Jx to A
        if(i > 0){
            int index0 = (N - 1) * dimx + (i - 1) * dimx;
            for(int j = 0; j < dimx; j++){
                A(aind) = -1.0;
                Ai(aind) = (i - 1) * dimx + j;
                aind += 1;
                for(int k = 0; k < dimx; k++){
                    A(aind) = J(k, 1 + j);
                    Ai(aind) = k + i * dimx;
                    aind += 1;
                }
                //Assign bounds for delta x
                A(aind) = 1.0;
                Ai(aind) = index0 + j;
                lbA(index0 + j) = dT->xlb(j) - tX(j, i);  //use i since I am using i
                ubA(index0 + j) = dT->xub(j) - tX(j, i);  //use i since I am using i
                aind += 1;
            }
        }
        // assign to u seems fixed so I will put it here
        for(int j = 0; j < dimu; j++){
            // Assign Ju
            for(int k = 0; k < dimx; k++){
                if(j == 0){
                    lbA(i*dimx + k) = right(k);
                    ubA(i*dimx + k) = right(k);
                }
                A(uind) = J(k, 1 + dimx + j);
                Ai(uind) = i * dimx + k;
                uind += 1;
            }
            // Assign bd
            A(uind) = 1.0;
            int urow = (N - 1) * dimx + (N - 2) * dimx + i*dimu + j;
            Ai(uind) = urow;
            lbA(urow) = dT->ulb(j) - tU(j, i);
            ubA(urow) = dT->uub(j) - tU(j, i);
            uind += 1;
        }
    }
    // Assign Ap
    Ap(0) = 0;
    int toappend = 0;
    for(int i = 0; i < nVar; i++){
        if(i < (N - 2) * dimx)
            toappend = dimx + 2;
        else
            toappend = dimx + 1;
        Ap(i + 1) = Ap(i) + toappend;
    }

    /* start solving the problem */
    data->P = csc_matrix(data->n, data->n, P_nnz, P.data(), pi.data(), pp.data());
    data->q = q.data();
    data->A = csc_matrix(data->m, data->n, A_nnz, A.data(), Ai.data(), Ap.data());
    data->l = lbA.data();
    data->u = ubA.data();
    // Define Solver settings as default
    osqp_set_default_settings(settings);
    settings -> verbose = false;

    // Setup workspace
    work = osqp_setup(data, settings);

    // Solve Problem
    osqp_solve(work);

    OSQPSolution *sol = work->solution;
    OSQPInfo *info = work->info;
    MapV Sol(sol->x, nVar);
    xin.segment(0, dimx) = dT->x0;
    xin.segment((N - 1)*dimx, dimx) = dT->xf;
    xin.segment(dimx, (N - 2) * dimx) += Sol.head((N - 2) * dimx);
    xin.segment(N*dimx, (N - 1) * dimu) += Sol.tail((N - 1) * dimu);
    double dX = Sol.head((N - 2) * dimx).norm();
    double dU = Sol.tail((N - 1) * dimu).norm();
    printf("||dX|| = %lf, ||dU|| = %lf\n", dX, dU);
    // Clean workspace
    osqp_cleanup(work);
    return 1;
}

MX acceGen(RefV xin, bool mpc){
    int dimx = 12, dimu = 4;
    int N = (xin.size() - 1 + dimu) / (dimx + dimu);
    if(mpc)
        N = (xin.size() + dimu) / (dimx + dimu);
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