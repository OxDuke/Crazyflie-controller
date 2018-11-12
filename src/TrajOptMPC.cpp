/*
 * Implement the class of trajectory optimization
 * The system should be defined and passed in
 */
#include "TrajOptMPC.h"
#include "json/json.h"

void funwrapper(const double t, cRefV y, cRefV p, RefV dy, RefM J);

namespace TrajOptMPC{
    dirTranODE::dirTranODE(tigerSys *sys, int N_, double tf_){
        N = N_;
        tf = tf_;
        system = sys;
        conFun = NULL;
        pFun = NULL;
        obsFun = NULL;
        dimx = sys->getdimx();
        dimu = sys->getdimu();
        numX = N * dimx;
        numU = (N - 1) * dimu;
        numVar = numX + numU;
        //set Q, R, F as eye
        Q = VX::Ones(dimx);
        R = VX::Ones(dimu);
        F = VX::Ones(dimx);
        enableRDU = false;
        incorrectRDU = false;
        RDU = VX::Zero(dimu);
        fixx0 = true;
        fixxf = false;  //Change default from previously false to true
        ubase = VX::Zero(dimu);
        xbase = VX::Zero(dimx);
        tfweight = 0;
    }

    dirTranODE::dirTranODE(Json::Value &jsin) {
        N = jsin["N"].asInt();
        tf = jsin["tf"].asDouble();
        dimx = jsin["dimx"].asInt();
        dimu = jsin["dimu"].asInt();
        numVar = N * dimx + (N - 1) * dimu;
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
        incorrectRDU = false;
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
        fixxf = false;
        if(jsin["fixxf"].asInt() > 0)
            fixxf = true;
        else
            fixxf = false;
    }
    /*
    dirTranODE::dirTranODE(YAML::Node &yn) {
        N = yn["N"].as<int>();
        tf = yn["tf"].as<double>();
        dimx = yn["dimx"].as<int>();
        dimu = yn["dimu"].as<int>();
        numVar = N * dimx + (N - 1) * dimu;
        numX = N * dimx;
        numU = (N - 1) * dimu;
        x0.setZero(dimx);
        xf.setZero(dimx);
        ubase.setZero(dimu);
        xbase.setZero(dimx);
        bool hasxbd = false, hasubd = false, hasubase = false;
        if(yn["xlb"]){
            hasxbd = true;
            xlb.setZero(dimx);
            xub.setZero(dimx);
        }
        if(yn["ulb"]){
            hasubd = true;
            ulb.setZero(dimu);
            uub.setZero(dimu);
        }
        if(yn["ubase"]){
            hasubase = true;
            ubase.setZero(dimu);
        }
        if(yn["RDU"]){
            enableRDU = true;
            RDU.resize(dimu);
            for(int i = 0; i < dimu; i++){
                RDU(i) = yn["RDU"][i].as<double>();
            }
        }
        else
            enableRDU = false;
        if(yn["incorrectRDU"]){
            if(yn["incorrectRDU"].as<int>() > 0)
                incorrectRDU = true;
            else
                incorrectRDU = false;
        }
        F.resize(dimx);
        Q.resize(dimx);
        R.resize(dimu);
        for (int i = 0; i < dimx; i++) {
            F(i) = yn["F"][i].as<double>();
            Q(i) = yn["Q"][i].as<double>();
            x0(i) = yn["x0"][i].as<double>();
            xf(i) = yn["xf"][i].as<double>();
            if(hasxbd){
                xlb(i) = yn["xlb"][i].as<double>();
                xub(i) = yn["xub"][i].as<double>();
            }
        }
        for (int i = 0; i < dimu; i++) {
            R(i) = yn["R"][i].as<double>();
            if(hasubd){
                ulb(i) = yn["ulb"][i].as<double>();
                uub(i) = yn["uub"][i].as<double>();
            }
            if(hasubase)
                ubase(i) = yn["ubase"][i].as<double>();
        }
        fixx0 = true;
        if(yn["fixxf"].as<int>() > 0)
            fixxf = true;
        else
            fixxf = false;
    }
    */
    /* The function called by SNOPT which return both obj, constraint, and gradients */
    //The correct positions of gradients have already been defined
    //void dirTranODE::userfun(traj &tj, double *F, double *G);

    void dirTranODE::objfun(traj &tj, double &y, RefM g, bool needg){
        VX dxf = tj.X.col(N - 1) - xf;
        y = 0;
        double h = (tj.t(N - 1) - tj.t(0)) / (N - 1);
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
                }
            }
            for (int j = 0; j < dimu; j++) {
                double tmp = du(j) * R(j) * du(j);
                y += tmp * h;
                if(needg){
                    g(0, numX + i*dimu + j) = 2*R(j)*du(j)*h;
                }
            }
            if(enableRDU){
                if(i == 0){
                    for(int j = 0; j < dimu; j++){
                        double tmp = du(j) * RDU(j) * du(j);
                        y += tmp / h;
                        if(needg){
                            g(0, numX + i*dimu + j) += 2*RDU(j)*du(j)/h;
                        }
                    }
                }
                else{
                    VX stepu = tj.U.col(i) - tj.U.col(i - 1);
                    for(int j = 0; j < dimu; j++){
                        double tmp = stepu(j) * RDU(j) * stepu(j);
                        y += tmp / h;
                        if(needg){
                            g(0, numX + i*dimu + j) += 2*RDU(j)*stepu(j)/h;
                            g(0, numX + (i - 1)*dimu + j) -= 2*RDU(j)*stepu(j)/h;
                        }
                    }
                }
            }
        }
        // penalize the last control, too if enableRDU is on
        if(enableRDU){
            VX stepu = tj.U.col(N - 2) - ubase;
            int i = N - 2;
            for(int j = 0; j < dimu; j++){
                double tmp = stepu(j) * RDU(j) * stepu(j);
                y += tmp / h;
                if(needg){
                    g(0, numX + i*dimu + j) += 2*RDU(j)*stepu(j)/h;
                }
            }
        }
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
            if(obsFun)
                nc += obsFun->dimc * (N - 1);
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
            }
        }
        /*** Path constraint***/
        int dimpath = 0;
        for(int pathi = 0; pathi < 2; pathi++){
            if(pathi == 0){
                if((!conFun) || (conFun->dimc == 0))
                    continue;
            }
            else{
                if((!obsFun) || (obsFun->dimc == 0))
                    continue;
            }
            int dimc = 0;
            if(pathi == 0)
                dimc = conFun->dimc;
            else
                dimc = obsFun->dimc;
            dimpath = dimc;
            VX vc(dimc);
            SpMX cg(dimc, dimx + dimu + 1);
            if(pathi == 0)
                cg.reserve(conFun->nnz);
            else
                cg.reserve(obsFun->nnz);
            int curRow = (N - 1) * dimx;
            for (int i = 0; i < N - 1; i++) {
                int xind = (i + 1) * dimx, uind = numX + i * dimu;
                double curt = i * h;
                if(pathi == 0)
                    conFun->eval(curt, tj.X.col(i + 1), tj.U.col(i), vc, cg);
                else
                    obsFun->eval(curt, tj.X.col(i + 1), tj.U.col(i), vc, cg);
                c.segment(curRow + rowadd, dimc).array() = vc.array();
                if(needg){
                    for (int k=0; k<cg.outerSize(); ++k)
                    {
                        for (Eigen::SparseMatrix<double>::InnerIterator it(cg,k); it; ++it)
                        {
                            int curcol = it.col();
                            if(rec){
                                row[nG] = rowadd + curRow + it.row();
                                int curcol = it.col();
                                if(curcol >= dimx){
                                    col[nG] = uind + curcol - dimx;
                                    if(curcol == dimx + dimu)
                                        col[nG] = numVar - 1;  //time
                                }
                                else
                                    col[nG] = xind + curcol;
                            }
                            if(curcol == dimx + dimu)
                                value[nG] = pFun->indx * it.value() / (N - 1);
                            else
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
            SpMX cg(dimp, dimx + dimu + 1);
            cg.reserve(pFun->nnz);
            int curRow = (N - 1) * dimx + dimpath * (N - 1);
            int xind = pFun->indx * dimx, uind = numX + pFun->indu * dimu;
            double evalt = pFun->indx * h;
            pFun->eval(evalt, tj.X.col(pFun->indx), tj.U.col(pFun->indu), vc, cg);  // I am assuming time is 0
            c.segment(curRow + rowadd, dimp).array() = vc.array();
            if(needg){
                for (int k=0; k<cg.outerSize(); ++k)
                {
                    for (Eigen::SparseMatrix<double>::InnerIterator it(cg,k); it; ++it)
                    {
                        int curcol = it.col();
                        if(rec){
                            row[nG] = rowadd + curRow + it.row();
                            if(curcol >= dimx){
                                col[nG] = uind + curcol - dimx;
                                if(curcol == dimx + dimu)
                                    col[nG] = numVar - 1;  //time
                            }
                            else
                                col[nG] = xind + curcol;
                        }
                        if(curcol == dimx + dimu)
                            value[nG] = pFun->indx * it.value() / (N - 1);
                        else
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
