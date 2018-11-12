main:
	g++ testModel.cpp qpUpdate.cpp Utility.cpp TrajOpt.cpp -I../include -ljsoncpp -losqp -std=c++11 -g -o main
qp:
	g++ testQPupdate.cpp qpUpdate.cpp Utility.cpp TrajOpt.cpp -I../include -ljsoncpp -losqp -std=c++11 -g -o main

vel:
	g++ calcVel.cpp qpUpdate.cpp Utility.cpp TrajOpt.cpp -I../include -ljsoncpp -losqp -std=c++11 -o main
