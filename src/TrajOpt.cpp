/*
 * Implement the class of trajectory optimization
 * The system should be defined and passed in
 */
#include "TrajOpt.h"
#include "json/json.h"

void funwrapper(const double t, cRefV y, cRefV p, RefV dy, RefM J);

namespace TrajOpt{
    dirTranODE::dirTranODE(tigerSys *sys, int N_, double tmin, double tmax){
        N = N_;
        tfmin = tmin;
        tfmax = tmax;
        system = sys;
        conFun = NULL;
        pFun = NULL;
        dimx = sys->getdimx();
        dimu = sys->getdimu();
        numX = N * dimx;
        numU = (N - 1) * dimu;
        numVar = numX + numU + 1;
        //set Q, R, F as eye
        Q = VX::Ones(dimx);
        R = VX::Ones(dimu);
        F = VX::Ones(dimx);
        fixx0 = true;
        fixxf = true;  //Change default from previously false to true
        ubase = VX::Zero(dimu);
        xbase = VX::Zero(dimx);
        tfweight = 0;
    }

    dirTranODE::dirTranODE(Json::Value &jsin) {
        N = jsin["N"].asInt();
        tfmin = jsin["tf"][0].asDouble();
        tfmax = jsin["tf"][1].asDouble();
        tfweight = jsin["tfweight"].asDouble();
        dimx = jsin["dimx"].asInt();
        dimu = jsin["dimu"].asInt();
        numVar = N * dimx + (N - 1) * dimu + 1;
        numX = N * dimx;
        numU = (N - 1) * dimu;
        x0.setZero(dimx);
        xf.setZero(dimx);
        ubase.setZero(dimu);
        xbase.setZero(dimx);
        bool hasxbd = false, hasubd = false, hasubase = false;
        if(jsin.isMember("xlb")){
            hasxbd = true;
            xlb.setZero(dimx);
            xub.setZero(dimx);
        }
        if(jsin.isMember("ulb")){
            hasubd = true;
            ulb.setZero(dimu);
            uub.setZero(dimu);
        }
        if(jsin.isMember("ubase")){
            hasubase = true;
            ubase.setZero(dimu);
        }
        if(jsin.isMember("RDU")){
            enableRDU = true;
            RDU.resize(dimu);
            for(int i = 0; i < dimu; i++){
                RDU(i) = jsin["RDU"][i].asDouble();
            }
        }
        else
            enableRDU = false;
        if(jsin.isMember("incorrectRDU")){
            if(jsin["incorrectRDU"].asInt() > 0)
                incorrectRDU = true;
            else
                incorrectRDU = false;
        }
        F.resize(dimx);
        Q.resize(dimx);
        R.resize(dimu);
        for (int i = 0; i < dimx; i++) {
            F(i) = jsin["F"][i].asDouble();
            Q(i) = jsin["Q"][i].asDouble();
            x0(i) = jsin["x0"][i].asDouble();
            xf(i) = jsin["xf"][i].asDouble();
            if(hasxbd){
                xlb(i) = jsin["xlb"][i].asDouble();
                xub(i) = jsin["xub"][i].asDouble();
            }
        }
        for (int i = 0; i < dimu; i++) {
            R(i) = jsin["R"][i].asDouble();
            if(hasubd){
                ulb(i) = jsin["ulb"][i].asDouble();
                uub(i) = jsin["uub"][i].asDouble();
            }
            if(hasubase)
                ubase(i) = jsin["ubase"][i].asDouble();
        }
        fixx0 = true;
        if(jsin["fixxf"].asInt() > 0)
            fixxf = true;
        else
            fixxf = false;
    }
    /* The function called by SNOPT which return both obj, constraint, and gradients */
    //The correct positions of gradients have already been defined
    //void dirTranODE::userfun(traj &tj, double *F, double *G);

    void dirTranODE::objfun(traj &tj, double &y, RefM g, bool needg){
        VX dxf = tj.X.col(N - 1) - xf;
        y = 0;
        double h = (tj.t(N - 1) - tj.t(0)) / (N - 1);
        if(needg)
            g(0, numVar - 1) = 0;
        for (int i = 0; i < dimx; i++) {
            if(needg)
                g(0, (N - 1)*dimx + i) = 2*F(i)*dxf(i);
            y += F(i)*dxf(i)*dxf(i);
        }
        //Loop over N - 1 states and control
        for (int i = 0; i < N - 1; i++) {
            VX dx = tj.X.col(i) - xf;
            VX du = tj.U.col(i) - ubase;
            for (int j = 0; j < dimx; j++) {
                double tmp = dx(j) * Q(j) * dx(j);
                y += tmp*h;
                if(needg){
                    g(0, i*dimx + j) = 2*Q(j)*dx(j)*h;
                    g(0, numVar - 1) += tmp / (N - 1);
                }
            }
            for (int j = 0; j < dimu; j++) {
                double tmp = du(j) * R(j) * du(j);
                y += tmp * h;
                if(needg){
                    g(0, numX + i*dimu + j) = 2*R(j)*du(j)*h;
                    g(0, numVar - 1) += tmp / (N - 1);
                }
            }
            if(enableRDU){
                if(i == 0){
                    for(int j = 0; j < dimu; j++){
                        double tmp = du(j) * RDU(j) * du(j);
                        if(incorrectRDU)
                            y += tmp * h;
                        else
                            y += tmp / h;
                        if(needg){
                            if(incorrectRDU){
                                g(0, numX + i*dimu + j) += 2*RDU(j)*du(j)*h;
                                g(0, numVar - 1) += tmp / (N - 1);
                            }
                            else{
                                g(0, numX + i*dimu + j) += 2*RDU(j)*du(j)/h;
                                g(0, numVar - 1) -= tmp / h / h / (N - 1);
                            }
                        }
                    }
                }
                else{
                    VX stepu = tj.U.col(i) - tj.U.col(i - 1);
                    for(int j = 0; j < dimu; j++){
                        double tmp = stepu(j) * RDU(j) * stepu(j);
                        if(incorrectRDU)
                            y += tmp * h;
                        else
                            y += tmp / h;
                        if(needg){
                            if(incorrectRDU){
                                g(0, numX + i*dimu + j) += 2*RDU(j)*stepu(j)*h;
                                g(0, numX + (i - 1)*dimu + j) -= 2*RDU(j)*stepu(j)*h;
                                g(0, numVar - 1) += tmp / (N - 1);
                            }
                            else{
                                g(0, numX + i*dimu + j) += 2*RDU(j)*stepu(j)/h;
                                g(0, numX + (i - 1)*dimu + j) -= 2*RDU(j)*stepu(j)/h;
                                g(0, numVar - 1) -= tmp / h / h / (N - 1);
                            }
                        }
                    }
                }
            }
        }
        //Penalty on final time
        y += tfweight * tj.t(N - 1);
        if(needg)
            g(0, numVar - 1) += tfweight; //I'm cheating, x^T Q x + u^T R u is not multiplied by h//Update on Aug 24, problem fixed
    }

    //confun evaluation, assume optionally we can add row
#ifdef __APPLE__
    typedef Eigen::Ref<Eigen::Matrix<long, -1, 1> > RefVi;
#endif
    void dirTranODE::confun(traj &tj, RefV c, RefVi row, RefVi col, RefV value, int rowadd, int nGadd, bool rec, bool needg){
        double h = (tj.t(N - 1) - tj.t(0)) / (N - 1);
        int nc = 0, nceq = dimx * (N - 1);
        if(conFun){
            nc = conFun->dimc * (N - 1);
        }
        else{
            nc = 0;
        }
        /*** Compute state equality constraint ***/
        int nG = nGadd; //already occupied by others
        VX y(dimx);
        MX J(dimx, 1 + dimx + dimu);
        for(int i = 0; i < N - 1; i++){
            int curRow = i*dimx, xind = i*dimx, uind = numX + i*dimu;
            double t0 = tj.t(i), tf = tj.t(i + 1);
            VX x = tj.X.col(i), u = tj.U.col(i);
            //RK4(fun, t0, tf, x, u, y, J);
            y.setZero();
            J.setZero();
            RK4(funwrapper, t0, tf, x.data(), u.data(), y.data(), J.data(), dimx, dimu);
            //Assign values
            c.segment(rowadd + i*dimx, dimx).array() = tj.X.col(i + 1).array() - y.array();
            //w.r.t X(:, i)
            if(needg){
                for (int j = 0; j < dimx; j++) {
                    for (int k = 0; k < dimx; k++) {
                        if(rec){
                            row[nG] = rowadd + curRow + k;
                            col[nG] = xind + j;
                        }
                        value[nG] = -J(k, 1 + j);
                        nG++;
                    }
                }
                //w.r.t X(:, i + 1)
                for (int j = 0; j < dimx; j++) {
                    if(rec){
                        row[nG] = rowadd + curRow + j;
                        col[nG] = xind + dimx + j;//X(i +1)
                    }
                    value[nG] = 1;
                    nG++;
                }
                //w.r.t U(:, i)
                for (int j = 0; j < dimu; j++) {
                    for (int k = 0; k < dimx; k++) {
                        if(rec){
                            row[nG] = rowadd + curRow + k;
                            col[nG] = uind + j;
                        }
                        value[nG] = -J(k, dimx + 1 + j);
                        nG++;
                    }
                }
                //w.r.t tf
                for (int j = 0; j < dimx; j++) {
                    if(rec){
                        row[nG] = rowadd + curRow + j;
                        col[nG] = numVar - 1;
                    }
                    value[nG] = -J(j, 0) / (N - 1);
                    nG++;
                }
            }
        }
        /*** Path constraint***/
        int dimpath = 0;
        if(conFun && conFun->dimc > 0){
            int dimc = conFun->dimc;
            dimpath = dimc;
            VX vc(dimc);
            SpMX cg(dimc, dimx + dimu);
            cg.reserve(conFun->nnz);
            int curRow = (N - 1) * dimx;
            for (int i = 0; i < N - 1; i++) {
                int xind = (i + 1) * dimx, uind = numX + i * dimu;
                conFun->eval(0, tj.X.col(i + 1), tj.U.col(i), vc, cg);
                c.segment(curRow + rowadd, dimc).array() = vc.array();
                if(needg){
                    for (int k=0; k<cg.outerSize(); ++k)
                    {
                        for (Eigen::SparseMatrix<double>::InnerIterator it(cg,k); it; ++it)
                        {
                            if(rec){
                                row[nG] = rowadd + curRow + it.row();
                                int curcol = it.col();
                                if(curcol >= dimx)
                                    col[nG] = uind + curcol - dimx;
                                else
                                    col[nG] = xind + curcol;
                            }
                            value[nG] = it.value();
                            nG++;
                        }
                    }
                }
                curRow += dimc;
            }
        }
        /* point constraints */
        if(pFun && pFun->dimc > 0){
            int dimp = pFun->dimc;
            VX vc(dimp);
            SpMX cg(dimp, dimx + dimu);
            cg.reserve(pFun->nnz);
            int curRow = (N - 1) * dimx + dimpath * (N - 1);
            int xind = pFun->indx * dimx, uind = numX + pFun->indu * dimu;
            pFun->eval(0, tj.X.col(pFun->indx), tj.U.col(pFun->indu), vc, cg);  // I am assuming time is 0
            c.segment(curRow + rowadd, dimp).array() = vc.array();
            if(needg){
                for (int k=0; k<cg.outerSize(); ++k)
                {
                    for (Eigen::SparseMatrix<double>::InnerIterator it(cg,k); it; ++it)
                    {
                        if(rec){
                            row[nG] = rowadd + curRow + it.row();
                            int curcol = it.col();
                            if(curcol >= dimx)
                                col[nG] = uind + curcol - dimx;
                            else
                                col[nG] = xind + curcol;
                        }
                        value[nG] = it.value();
                        nG++;
                    }
                }
            }
            curRow += dimp;
        }
    }
    //Wrapper for RK4, evaluate those
    void dirTranODE::RKDynFun(const double t, cRefV x, cRefV u, RefV f, RefM df){
        system->dyn(t, x, u, f, df);
    }
};
