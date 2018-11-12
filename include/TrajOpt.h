/* The header file for trajectory optimization by direct methods
 * This is designed as the basis of LQRTree algorithm
 */
#ifndef TRAJOPT_H
#define TRAJOPT_H
#include "TigerTools/TigerEigen.h"
#include "Tools/VecMat.h"
#include "Utility.h"
#include "json/json.h"
namespace TrajOpt{
/* A basis class which is inherited by others
 * This is equivalent to the TigerSys class in Matlab
 */
    class tigerSys{
        private:
            int dimx, dimu, dimg;
        public:
            tigerSys(int dimx_, int dimu_, int dimg_ = 0):dimx(dimx_), dimu(dimu_), dimg(dimg_){}
            int getdimx(){return dimx;}
            int getdimu(){return dimu;}
            virtual void dyn(const double t, cRefV x, cRefV u, RefV f, RefM df){};
    };
/* A functor which evaluate constraints */
    class constraintFunctor{
        public:
            int dimstate, dimctrl, dimc;
            int nnz;
            VX clb;
            VX cub;
            virtual void eval(const double t, cRefV x, cRefV u, RefV f, SpMX &fx) = 0;
    };
/* A functor for complex state and control constraints, include indexes */
    class pointConstrFunctor{
        public:
            int dimstate, dimctrl, dimc;  // dimensions
            int nnz;  // sparse jacobian
            int indx, indu; //indexes of x and u in the traj
            VX clb, cub;  // lower and upper bound
            virtual void eval(const double t, cRefV x, cRefV u, RefV f, SpMX &fx) = 0;
    };
/* A class which is a trajectory */
    class traj{
        public:
            MX X;
            MX U;
            rVX t;
            traj(){}
            /* Construct directly by solution of local optimizer */
            traj(double *x, int N, int dimx, int dimu){
                MapM X_(x, dimx, N);
                X = X_;
                MapM U_(x + dimx*N, dimu, N - 1);
                U = U_;
                t = rVX::LinSpaced(N, 0, x[dimx*N + dimu*(N - 1)]);
            }
            /* Construct empty ones */
            traj(int N, int dimx, int dimu){
                X = MX::Zero(dimx, N);
                U = MX::Zero(dimu, N - 1);
                t = rVX::Zero(N);
            }
            /* Copy the instance back to x*/
            void copyto(double *x){
                int N = X.cols();
                int dimx = X.rows();
                int dimu = U.rows();
                V_Copy(x, X.data(), N*dimx);
                V_Copy(x+N*dimx, U.data(), dimu*(N-1));
                x[N*dimx+dimu*(N-1)] = t(N-1) - t(0);
                //MapV mV(x, N*dimx + (N - 1)*dimu + 1);
                //mV.segment(0, N*dimx).array() = X.array();
                //mV.segment(N*dimx, (N-1)*dimu) = U.array();
                //mV(N*dimx + (N-1)*dimu) = t(N - 1) - t(0);
            }
            /* Construct traj according to guess x */
            void copyfrom(const double *x){
                int N = X.cols();
                int dimx = X.rows();
                int dimu = U.rows();
                V_Copy(X.data(), x, N*dimx);
                V_Copy(U.data(), x + N*dimx, (N-1)*dimu);
                //U.col(N - 1) = U.col(N - 2);
                t.setLinSpaced(N, 0, x[N*dimx + (N-1)*dimu]);
            }
    };
/* A class to optimize a trajectory */
    class dirTranODE{
        public:
            tigerSys *system;
            int N; //Number of discretization
            double tfmin, tfmax;
            VX Q, R, F;
            VX RDU; // RDU means penalty on delta U, this is introduced since the drone motor is slow to ramp up
            bool enableRDU;
            bool incorrectRDU;
            VX x0, xf;
            VX xlb, xub, xbase;
            bool fixx0, fixxf;
            VX ulb, uub, ubase;
            int dimx, dimu, numVar, numX, numU;
            constraintFunctor *conFun = nullptr;
            pointConstrFunctor *pFun = nullptr;
            constraintFunctor *obsFun = nullptr;
            double tfweight; //This is used for time-optimal problem
            RKFun fun;

            dirTranODE(tigerSys *sys, int N_, double tmin, double tmax);

            dirTranODE(Json::Value &js);
            /* The function called by SNOPT which return both obj, constraint, and gradients */
            //The correct positions of gradients have already been defined
            void userfun(traj &tj, double *F, double *G);

            //Compute objective function
            void objfun(traj &tj, double &y, RefM g, bool needg=true);

            //Compute constraint function
#ifdef __APPLE__
            typedef Eigen::Ref<Eigen::Matrix<long, -1, 1> > RefVi;
#endif
            void confun(traj &tj, RefV c, RefVi row, RefVi col, RefV value, int rowadd, int nGadd, bool record, bool needg=true);

            //Wrapper for dynamic equation
            void RKDynFun(const double, cRefV y, cRefV p, RefV dy, RefM J);

            //Init a traj
            void ranGenGuess(double *x){
                x[numVar - 1] = (tfmin + tfmax) / 2.0 - (double)rand()/RAND_MAX * (tfmax - tfmin) / 2.0;
                MapM vX(x, dimx, N);
                for (int i = 0; i < dimx; i++) {
                    vX.row(i).setLinSpaced(N, x0(i), xf(i));
                }
                MapM vU(x + numX, dimu, N - 1);
                VX udiff = uub - ulb;
                for (int i = 0; i < dimu; i++) {
                    vU.row(i).setRandom().array() *= 0.2*udiff(i);
                    vU.row(i).array() += ubase(i);
                }
            }

            //Init a traj randomly, output is not x but traj
            traj ranGentj(){
                traj tj(N, dimx, dimu);
                tj.t.setLinSpaced(N, 0, (tfmin + tfmax)/2.);
                for (int i = 0; i < dimx; i++) {
                    tj.X.row(i).setLinSpaced(N, x0(i), xf(i));
                }
                VX udiff = uub - ulb;
                for (int i = 0; i < dimu; i++) {
                    tj.U.row(i).setRandom().array() *= udiff(i);
                    tj.U.row(i).array() += ubase(i);
                }
                return tj;
            }

            void setX0(cRefV x0_){
                x0 = x0_;
            }

            void setXf(cRefV xf_){
                xf = xf_;
            }

            void setulb(cRefV ulb_){
                ulb = ulb_;
            }

            void setuub(cRefV uub_){
                uub = uub_;
            }

            void setxlb(cRefV xlb_){
                xlb = xlb_;
            }

            void setxub(cRefV xub_){
                xub = xub_;
            }

    };

};
typedef TrajOpt::dirTranODE dirTranODE;
typedef TrajOpt::traj traj;
typedef TrajOpt::tigerSys tigerSys;
#endif
