#include "modelLoader.h"
#include "TigerEigen.h"
#include "stdio.h"
#include "qpUpdate.h"

int N = 20, dimx = 12, dimu = 4;

int main(){
	mlpModelQP GaoModel("/home/sun/Downloads/anyVel_6_200_317.json");
	double xin[6] = {0, 0, -0.5, 0, 0, 0};
	VX out = GaoModel.evalQP(xin);
	MapM tX(out.data(), dimx, N);
	MapM tU(out.data() + N*dimx, dimu, N - 1);
	double dt = out.tail(1)(0) / (N - 1);
	for(int i = 0; i < N; i++){
		printf("%lf %lf %lf %lf\n", i*dt, tX(0, i), tX(1, i), tX(2, i));
	}
	MX acce = acceGen(out);
	std::cout << acce << std::endl;
	return 1;
}