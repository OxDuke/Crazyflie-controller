/* The header file for trajectory optimization by direct methods
 * This is designed as the basis of LQRTree algorithm
 */
#ifndef TRAJOPTMPC_H
#define TRAJOPTMPC_H
#include "TigerTools/TigerEigen.h"
#include "Tools/VecMat.h"
#include "Utility.h"
#include "json/json.h"
#include "yaml-cpp/yaml.h"
#include "TrajOpt.h"
namespace TrajOptMPC{
/* A basis class which is inherited by others
 * This is equivalent to the TigerSys class in Matlab
 */
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
/* A class to optimize a trajectory */
    class dirTranODE{
        public:
            tigerSys *system;
            int N; //Number of discretization
            double tf;
            double tfmin, tfmax;  //deprecated
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
            constraintFunctor *obsFun = nullptr;
            pointConstrFunctor *pFun = nullptr;
            double tfweight; //This is used for time-optimal problem
            RKFun fun;

            dirTranODE(tigerSys *sys, int N_, double t);

            dirTranODE(Json::Value &js);
            dirTranODE(YAML::Node &yn);
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
                tj.t.setLinSpaced(N, 0, tf);
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
#endif
