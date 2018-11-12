/*
 * testQPupdate.cpp
 * Copyright (C) 2018 Gao Tang <gt70@duke.edu>
 *
 * Distributed under terms of the MIT license.
 */

/*  testQPupdate.cpp: test if we can use QP to update one solution  */
#include "cnpy.h"
#include "qpUpdate.h"
#include "utilityClass.h"
#include <iostream>
#include "stdlib.h"
#include "time.h"
#include <fstream>
#include "modelLoader.h"


dirTranODE *dT = nullptr;
//define the function to be integrated here


int main(){
    srand(time(NULL));
    mlpModelQP GaoModel("/home/sun/Downloads/rmuAnyVelWeightRDU_7_500_317.json");
    double weight = 1.0;
    double xin[7] = { -3, -3, -1.5, 0, 0, 0, 0}; // initial point
    VX vec = GaoModel.evalQP(xin);
    for(int i = 0; i < 20; i++)
        printf("%lf\n", vec(i*12 + 2));
    return 1;
}
