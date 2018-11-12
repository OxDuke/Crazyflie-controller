/*
 * modelLoader.h
 * Copyright (C) 2017 Gao Tang <gt70@duke.edu>
 *
 * Distributed under terms of the MIT license.
 */

/* Previously I export a model into a json file, now I directly parse that file and get a powerful function for inference */

#ifndef MODELLOADER_H
#define MODELLOADER_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "stdlib.h"
#include "TigerEigen.h"
#include "json/json.h"
#include "qpUpdate.h"
#include "qpMPCUpdate.h"


class mlpModel{
private:
    std::vector<MX*> w;
    std::vector<VX*> b;
    double leaky;
    VX xmean, xstd, ymean, ystd;
    int dimIn, dimOut, nConnect;
public:
    mlpModel(std::string fnm){
        std::ifstream fin(fnm, std::ios::in);
        if(!fin.is_open()){
            std::cout << fnm << " open fail\n";
            exit(0);
        }
        Json::Value jsin;
        fin >> jsin;
        fin.close();
        // parse the json file to construct the model
        dimIn = jsin["dimIn"].asInt();
        dimOut = jsin["dimOut"].asInt();
        nConnect = jsin["nConnect"].asInt();
        leaky = jsin["leaky"].asDouble();
        // read in mean and std for x and y
        xmean.resize(dimIn);
        xstd.resize(dimIn);
        for(int i = 0; i < dimIn; i++){
            xmean(i) = jsin["xmean"][i].asDouble();
            xstd(i) = jsin["xstd"][i].asDouble();
        }
        ymean.resize(dimOut);
        ystd.resize(dimOut);
        for(int i = 0; i < dimOut; i++){
            ymean(i) = jsin["ymean"][i].asDouble();
            ystd(i) = jsin["ystd"][i].asDouble();
        }
        // read in weight and bias
        int prevRow = dimIn;
        for(int i = 0; i < nConnect; i++){
            // read b0, first get its size
            std::string bkey = std::string("b") + std::to_string(i);
            std::string wkey = std::string("w") + std::to_string(i);
            int bsz = jsin[bkey].size();
            VX *lyrb = new VX(bsz);
            MX *lyrw = new MX(bsz, prevRow);
            for(int bj = 0; bj < bsz; bj++){
                (*lyrb)(bj) = jsin[bkey][bj].asDouble();
                for(int wi = 0; wi < prevRow; wi++){
                    (*lyrw)(bj, wi) = jsin[wkey][bj][wi].asDouble();
                }
            }
            b.push_back(lyrb);
            w.push_back(lyrw);
            prevRow = bsz;
        }
    }

    VX eval(double *xin){
        MapV x(xin, dimIn);
        //shift
        VX tmp = x - xmean;
        tmp.array() /= xstd.array();
        VX out(dimOut);
        for(int i = 0; i < nConnect; i++){
            tmp = (*(w[i])) * tmp + (*(b[i]));
            if(i < nConnect - 1)
                doleaky(tmp);
        }
        //shift
        out = tmp.array() * ystd.array() + ymean.array();
        return out;
    }

    void eval(double *xin, double *xout){
        MapV x(xin, dimIn);
        //shift
        VX tmp = x - xmean;
        tmp.array() /= xstd.array();
        MapV out(xout, dimOut);
        for(int i = 0; i < nConnect; i++){
            tmp = (*(w[i])) * tmp + (*(b[i]));
            if(i < nConnect - 1)
                doleaky(tmp);
        }
        //shift
        out = tmp.array() * ystd.array() + ymean.array();
    }

    void doleaky(RefV vx){
        int sz = vx.size();
        for(int i = 0; i < sz; i++){
            if(vx(i) < 0)
                vx(i) *= leaky;
        }
    }

    ~mlpModel(){
        for(auto w_ : w)
            delete w_;
        for(auto b_ : b)
            delete b_;
    }

};

class mlpModelQP : public mlpModel{
public:
    qpUpdater qpUp;

    mlpModelQP(std::string fnm) : mlpModel(fnm){
    }

    VX evalQP(double *xin){
        VX out = eval(xin);
        VX x0 = VX::Zero(12);
        for(int i = 0; i < 3; i++){
            x0(i) = xin[i];
            x0(6 + i) = xin[3 + i];
        }
        qpUp.qpUpdateOSQP(out, x0);
        return out;
    }

};

class mlpModelMPCQP : public mlpModel{
public:
    qpMPC::qpMPCUpdater qpUp;

    mlpModelMPCQP(std::string fnm) : mlpModel(fnm){
        // Test it here
        double xin[6] = {1, 1, 1, 0, 0, 0};
        VX out = evalQP(xin);
        MX acce = qpMPC::acceGen(out, true);
    }

    VX evalQP(double *xin){
        VX out = eval(xin);
        VX x0 = VX::Zero(12);
        for(int i = 0; i < 3; i++){
            x0(i) = xin[i];
            x0(6 + i) = xin[3 + i];
        }
        qpUp.qpMPCUpdateOSQP(out, x0);
        std::cout << "returning traj\n";
        return out;
    }

};

#endif /* !MODELLOADER_H */
