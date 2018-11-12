#ifndef QPMPCUPDATE_H
#define QPMPCUPDATE_H

#include "TrajOptMPC.h"
#include "TrajOpt.h"
#include "TigerEigen.h"
#include "utilityClass.h"
#include <string>
#include <fstream>

namespace qpMPC{
	class qpMPCUpdater{
	public:
		TrajOptMPC::dirTranODE *dT;
		Rotor sys;
		qpMPCUpdater(){
			//Load the configuration file
			std::ifstream fin("/home/sun/crazyflie_ws/src/crazyflie_ros/crazyflie_controller/config/droneConfigMPC.json", std::ios::in);
		    Json::Value jsin;
		    fin >> jsin;
		    fin.close();
		    dT = new TrajOptMPC::dirTranODE(jsin);
		    dT->system = &sys;
		}
		int qpMPCUpdateOSQP(RefV xin, cRefV x0);
		int acceGen(cRefV xin, RefV aout);
		~qpMPCUpdater(){
			delete dT;
		}
		
	};

	MX acceGen(RefV xin, bool mpc=true);
};

int qpUpdateOSQP(TrajOptMPC::dirTranODE *dT, RefV xin, RefV xout);
#endif /* !QPUPDATE_H */