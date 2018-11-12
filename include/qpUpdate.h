/*
 * qpUpdate.h
 * Copyright (C) 2018 Gao Tang <gt70@duke.edu>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef QPUPDATE_H
#define QPUPDATE_H

#include "TigerEigen.h"
#include "utilityClass.h"
#include <string>
#include <fstream>
#include "yaml-cpp/yaml.h"


class qpUpdater{
public:
	TrajOpt::dirTranODE *dT;
	Rotor sys;
	qpUpdater(){
		//Load the configuration file
		std::ifstream fin("/home/sun/crazyflie_ws/src/crazyflie_ros/crazyflie_controller/config/droneConfig.json", std::ios::in);
	    Json::Value jsin;
	    fin >> jsin;
	    fin.close();
	    dT = new dirTranODE(jsin);
	    dT->system = &sys;
	}
	void updateWeight(double weight){
		dT->tfweight = weight; // I do not need it since tf is not optimized here
	}
	int qpUpdateOSQP(RefV xin, cRefV x0);
	int acceGen(cRefV xin, RefV aout);
	~qpUpdater(){
		delete dT;
	}
};

MX acceGen(RefV xin, bool mpc=false);

#endif /* !QPUPDATE_H */
