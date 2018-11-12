//
// Created by Gao Tang on 8/22/17.
//

#ifndef ROTORDB_UTILITYCLASS_H
#define ROTORDB_UTILITYCLASS_H

#include "TrajOpt.h"
#include "TigerTools/TigerEigen.h"
#include <vector>
#include "json/json.h"

/* Define a class for system dynamics */
class Rotor : public tigerSys {
private:
    int dimx = 12, dimu = 4;
    double m = 0.0335, g = 9.81, L = 0.046;
    double In[3][3] = {{16.571710e-6, 0.830806e-6, 0.718277e-6}, {0.830806e-6, 16.655602e-6, 1.800197e-6}, {0.718277e-6, 1.800197e-6, 29.261652e-6}};
    double invIn[3][3] = {{60544.8576958428,-2878.57629514431,-1309.08447480706},{-2878.57629514431,60578.6197781223,-3656.17889390418},{-1309.08447480706,-3656.17889390418,34431.4848507385}};
    double a = 0.3140636424 / 4.0, b = 0.3112328233 / 4.0, c = -0.3744188e-2 / 4.0, d = 0.005964552, e = 1.563383e-5;
    double cg0[204] = {0.0};
public:
    Rotor():tigerSys(12, 4){};
    virtual void dyn(const double t, cRefV x, cRefV u, RefV f, RefM df){
        //First we call dronedyn to flush cg0
        dronedyn(t, x.data(), u.data());
        df.setZero();
        MapV dx(cg0, dimx);
        f = dx;
        MapM J(cg0 + dimx, dimx, dimx + dimu);
        df.middleCols(1, dimx + dimu) = J;
    }
    void dyn(const double tt, cRefV x, cRefV u, RefV f){
        double phi = x[3], theta = x[4], psi = x[5], xd = x[6], yd = x[7], zd = x[8], pp = x[9], q = x[10], r = x[11];
        const double *t = u.data();
        double *cg1 = f.data();
        double t1 = cos(theta);
        double t2 = sin(theta);
        double t3 = sin(phi);
        double t4 = cos(phi);
        double t5 = 0.1e1 / t4;
        t5 = t5 * (pp * t2 - r * t1);
        double t6 = sin(psi);
        double t7 = cos(psi);
        double t8 = t1 * t3;
        double t9 = pow(t[1], 0.2e1);
        double t10 = pow(t[2], 0.2e1);
        double t11 = pow(t[0], 0.2e1);
        double t12 = pow(t[3], 0.2e1);
        double t13 = (double) (4 * c) + (t11 + t12 + t9 + t10) * a + (t[0] + t[1] + t[2] + t[3]) * b;
        double t14 = 0.1e1 / m;
        double t15 = In[2][0] * pp + In[2][1] * q + r * In[2][2];
        double t16 = In[1][0] * pp + q * In[1][1] + In[1][2] * r;
        t9 = L * (-a * t12 + a * t9 + b * t[1] - b * t[3]) - q * t15 + r * t16;
        t12 = pp * In[0][0] + In[0][1] * q + In[0][2] * r;
        t10 = L * (-a * t10 + a * t11 + b * t[0] - b * t[2]) - pp * t15 + r * t12;
        t11 = ((a * t[0] + b) * t[0] - (a * t[1] + b) * t[1] + (a * t[2] + b) * t[2] - (a * t[3] + b) * t[3]) * d - pp * t16 + q * t12;
        cg1[0] = xd;
        cg1[1] = yd;
        cg1[2] = zd;
        cg1[3] = pp * t1 + r * t2;
        cg1[4] = t5 * t3 + q;
        cg1[5] = -t5;
        cg1[6] = t14 * (t2 * t7 + t8 * t6) * t13;
        cg1[7] = t14 * (t2 * t6 - t8 * t7) * t13;
        cg1[8] = t14 * t4 * t1 * t13 - g;
        cg1[9] = -t10 * invIn[0][1] + t11 * invIn[0][2] + t9 * invIn[0][0];
        cg1[10] = -t10 * invIn[1][1] + t11 * invIn[1][2] + t9 * invIn[1][0];
        cg1[11] = -t10 * invIn[2][1] + t11 * invIn[2][2] + t9 * invIn[2][0];
    }
    void dronedyn(const double tt, const double *x, const double *u){
        double phi = x[3], theta = x[4], psi = x[5], xd = x[6], yd = x[7], zd = x[8], pp = x[9], q = x[10], r = x[11];
        const double *t = u;
        double t1 = cos(theta);
        double t2 = sin(theta);
        double t3 = pp * t1 + r * t2;
        double t4 = sin(phi);
        double t5 = cos(phi);
        double t6 = 0.1e1 / t5;
        double t7 = t1 * r;
        double t8 = t2 * pp;
        double t9 = t8 - t7;
        double t10 = t6 * t9;
        double t11 = cos(psi);
        double t12 = sin(psi);
        double t13 = t1 * t12;
        double t14 = t11 * t2;
        double t15 = pow(t[2], 0.2e1);
        double t16 = pow(t[0], 0.2e1);
        double t17 = pow(t[1], 0.2e1);
        double t18 = pow(t[3], 0.2e1);
        double t19 = (4.0 * c) + (t18 + t17 + t16 + t15) * a + (t[0] + t[1] + t[2] + t[3]) * b;
        t11 = t11 * t1;
        t12 = t12 * t2;
        double t20 = -t11 * t4 + t12;
        double t21 = 0.1e1 / m;
        t5 = t21 * t5;
        double t22 = t5 * t1;
        double t23 = In[2][0] * pp;
        double t24 = In[2][1] * q;
        double t25 = r * In[2][2] + t23 + t24;
        double t26 = In[1][0] * pp;
        double t27 = In[1][2] * r;
        double t28 = q * In[1][1] + t26 + t27;
        t17 = L * (a * t17 - a * t18 + b * t[1] - b * t[3]) - q * t25 + r * t28;
        t18 = In[0][1] * q;
        double t29 = In[0][2] * r;
        double t30 = pp * In[0][0] + t18 + t29;
        t15 = L * (-a * t15 + a * t16 + b * t[0] - b * t[2]) - pp * t25 + r * t30;
        t16 = a * t[0];
        t25 = a * t[1];
        double t31 = a * t[2];
        double t32 = a * t[3];
        t28 = ((b + t16) * t[0] - (b + t25) * t[1] + (b + t31) * t[2] - (b + t32) * t[3]) * d - pp * t28 + q * t30;
        t30 = -t15 * invIn[1][1] + t17 * invIn[1][0] + t28 * invIn[1][2];
        double t33 = t15 * invIn[2][1] - t17 * invIn[2][0] - t28 * invIn[2][2];
        double t34 = pow(t6, 0.2e1);
        double t35 = t34 * pow(t4, 0.2e1) + 0.1e1;
        double t36 = t6 * t3;
        double t37 = q * In[2][0] - r * In[1][0];
        double t38 = In[0][0] - In[2][2];
        double t39 = -r * t38 + 0.2e1 * t23 + t24;
        double t40 = -In[0][0] + In[1][1];
        double t41 = -q * t40 - 0.2e1 * t26 - t27;
        double t42 = t37 * invIn[0][0] - t39 * invIn[0][1] - t41 * invIn[0][2];
        double t43 = t37 * invIn[1][0] - t39 * invIn[1][1] - t41 * invIn[1][2];
        t37 = t37 * invIn[2][0] - t39 * invIn[2][1] - t41 * invIn[2][2];
        t39 = In[1][1] - In[2][2];
        t23 = r * t39 - t23 - 0.2e1 * t24;
        t24 = pp * In[2][1] - r * In[0][1];
        t40 = -pp * t40 + 0.2e1 * t18 + t29;
        t41 = t23 * invIn[0][0] + t24 * invIn[0][1] + t40 * invIn[0][2];
        double t44 = t23 * invIn[1][0] + t24 * invIn[1][1] + t40 * invIn[1][2];
        t26 = q * t39 + t26 + 0.2e1 * t27;
        t18 = -pp * t38 - t18 - 0.2e1 * t29;
        t27 = pp * In[1][2] - q * In[0][2];
        t29 = t18 * invIn[0][1] + t26 * invIn[0][0] - t27 * invIn[0][2];
        t38 = t18 * invIn[1][1] + t26 * invIn[1][0] - t27 * invIn[1][2];
        t18 = t18 * invIn[2][1] + t26 * invIn[2][0] - t27 * invIn[2][2];
        t16 = b + 0.2e1 * t16;
        t26 = -t16;
        t27 = invIn[0][1] * L;
        t39 = invIn[0][2] * d;
        double t45 = t39 * t16 + t27 * t26;
        double t46 = invIn[1][2] * d;
        double t47 = invIn[1][1] * L;
        double t48 = invIn[2][2] * d;
        double t49 = invIn[2][1] * L;
        double t50 = t48 * t16 + t49 * t26;
        t25 = b + 0.2e1 * t25;
        double t51 = invIn[0][0] * L;
        double t52 = t25 * (t51 - t39);
        double t53 = invIn[1][0] * L;
        double t54 = invIn[2][0] * L;
        double t55 = t25 * (t54 - t48);
        t31 = b + 0.2e1 * t31;
        t27 = t31 * (t27 + t39);
        t49 = t31 * (t49 + t48);
        t32 = b + 0.2e1 * t32;
        double t56 = -t32;
        t39 = -t39 * t32 + t51 * t56;
        t51 = -t46 * t32 + t53 * t56;
        t48 = -t48 * t32 + t54 * t56;
        t54 = t21 * (t13 * t4 + t14);
        t56 = t54 * t19;
        double t57 = t21 * t20;
        double t58 = t54 * t32;
        double t59 = t54 * t25;
        t20 = t21 * t20 * t19;
        t14 = t21 * (t14 * t4 + t13) * t19;
        t12 = t21 * (-t12 * t4 + t11) * t19;
        double t60 = t1 * t6;
        t6 = t2 * t6;
        cg0[0] = xd;
        cg0[1] = yd;
        cg0[2] = zd;
        cg0[3] = t3;
        cg0[4] = t10 * t4 + q;
        cg0[5] = -t10;
        cg0[6] = t56;
        cg0[7] = t57 * t19;
        cg0[8] = t22 * t19 - g;
        cg0[9] = -t15 * invIn[0][1] + t17 * invIn[0][0] + t28 * invIn[0][2];
        cg0[10] = t30;
        cg0[11] = -t33;
        cg0[12] = 0;
        cg0[13] = 0;
        cg0[14] = 0;
        cg0[15] = 0;
        cg0[16] = 0;
        cg0[17] = 0;
        cg0[18] = 0;
        cg0[19] = 0;
        cg0[20] = 0;
        cg0[21] = 0;
        cg0[22] = 0;
        cg0[23] = 0;
        cg0[24] = 0;
        cg0[25] = 0;
        cg0[26] = 0;
        cg0[27] = 0;
        cg0[28] = 0;
        cg0[29] = 0;
        cg0[30] = 0;
        cg0[31] = 0;
        cg0[32] = 0;
        cg0[33] = 0;
        cg0[34] = 0;
        cg0[35] = 0;
        cg0[36] = 0;
        cg0[37] = 0;
        cg0[38] = 0;
        cg0[39] = 0;
        cg0[40] = 0;
        cg0[41] = 0;
        cg0[42] = 0;
        cg0[43] = 0;
        cg0[44] = 0;
        cg0[45] = 0;
        cg0[46] = 0;
        cg0[47] = 0;
        cg0[48] = 0;
        cg0[49] = 0;
        cg0[50] = 0;
        cg0[51] = 0;
        cg0[52] = -t7 * t35 + t8 * t35;
        cg0[53] = -t34 * t4 * t9;
        cg0[54] = t5 * t13 * t19;
        cg0[55] = -t5 * t11 * t19;
        cg0[56] = -t21 * t4 * t1 * t19;
        cg0[57] = 0;
        cg0[58] = 0;
        cg0[59] = 0;
        cg0[60] = 0;
        cg0[61] = 0;
        cg0[62] = 0;
        cg0[63] = -t9;
        cg0[64] = t36 * t4;
        cg0[65] = -t36;
        cg0[66] = t12;
        cg0[67] = t14;
        cg0[68] = -t5 * t2 * t19;
        cg0[69] = 0;
        cg0[70] = 0;
        cg0[71] = 0;
        cg0[72] = 0;
        cg0[73] = 0;
        cg0[74] = 0;
        cg0[75] = 0;
        cg0[76] = 0;
        cg0[77] = 0;
        cg0[78] = -t20;
        cg0[79] = t56;
        cg0[80] = 0;
        cg0[81] = 0;
        cg0[82] = 0;
        cg0[83] = 0;
        cg0[84] = 1;
        cg0[85] = 0;
        cg0[86] = 0;
        cg0[87] = 0;
        cg0[88] = 0;
        cg0[89] = 0;
        cg0[90] = 0;
        cg0[91] = 0;
        cg0[92] = 0;
        cg0[93] = 0;
        cg0[94] = 0;
        cg0[95] = 0;
        cg0[96] = 0;
        cg0[97] = 1;
        cg0[98] = 0;
        cg0[99] = 0;
        cg0[100] = 0;
        cg0[101] = 0;
        cg0[102] = 0;
        cg0[103] = 0;
        cg0[104] = 0;
        cg0[105] = 0;
        cg0[106] = 0;
        cg0[107] = 0;
        cg0[108] = 0;
        cg0[109] = 0;
        cg0[110] = 1;
        cg0[111] = 0;
        cg0[112] = 0;
        cg0[113] = 0;
        cg0[114] = 0;
        cg0[115] = 0;
        cg0[116] = 0;
        cg0[117] = 0;
        cg0[118] = 0;
        cg0[119] = 0;
        cg0[120] = 0;
        cg0[121] = 0;
        cg0[122] = 0;
        cg0[123] = t1;
        cg0[124] = t6 * t4;
        cg0[125] = -t6;
        cg0[126] = 0;
        cg0[127] = 0;
        cg0[128] = 0;
        cg0[129] = -t42;
        cg0[130] = -t43;
        cg0[131] = -t37;
        cg0[132] = 0;
        cg0[133] = 0;
        cg0[134] = 0;
        cg0[135] = 0;
        cg0[136] = 1;
        cg0[137] = 0;
        cg0[138] = 0;
        cg0[139] = 0;
        cg0[140] = 0;
        cg0[141] = t41;
        cg0[142] = t44;
        cg0[143] = t23 * invIn[2][0] + t24 * invIn[2][1] + t40 * invIn[2][2];
        cg0[144] = 0;
        cg0[145] = 0;
        cg0[146] = 0;
        cg0[147] = t2;
        cg0[148] = -t60 * t4;
        cg0[149] = t60;
        cg0[150] = 0;
        cg0[151] = 0;
        cg0[152] = 0;
        cg0[153] = t29;
        cg0[154] = t38;
        cg0[155] = t18;
        cg0[156] = 0;
        cg0[157] = 0;
        cg0[158] = 0;
        cg0[159] = 0;
        cg0[160] = 0;
        cg0[161] = 0;
        cg0[162] = t54 * t16;
        cg0[163] = t57 * t16;
        cg0[164] = t22 * t16;
        cg0[165] = t45;
        cg0[166] = t46 * t16 + t47 * t26;
        cg0[167] = t50;
        cg0[168] = 0;
        cg0[169] = 0;
        cg0[170] = 0;
        cg0[171] = 0;
        cg0[172] = 0;
        cg0[173] = 0;
        cg0[174] = t59;
        cg0[175] = t57 * t25;
        cg0[176] = t22 * t25;
        cg0[177] = t52;
        cg0[178] = t25 * (t53 - t46);
        cg0[179] = t55;
        cg0[180] = 0;
        cg0[181] = 0;
        cg0[182] = 0;
        cg0[183] = 0;
        cg0[184] = 0;
        cg0[185] = 0;
        cg0[186] = t54 * t31;
        cg0[187] = t57 * t31;
        cg0[188] = t22 * t31;
        cg0[189] = t27;
        cg0[190] = t31 * (t47 + t46);
        cg0[191] = t49;
        cg0[192] = 0;
        cg0[193] = 0;
        cg0[194] = 0;
        cg0[195] = 0;
        cg0[196] = 0;
        cg0[197] = 0;
        cg0[198] = t58;
        cg0[199] = t57 * t32;
        cg0[200] = t22 * t32;
        cg0[201] = t39;
        cg0[202] = t51;
        cg0[203] = t48;
    }
};

class geometry{
public:
    bool out = true;  // means we want to be outside
    int ngrad = 0;
    geometry(bool out_=true):out(out_){};
    virtual void evalConstr(const double *s, double *val, double *grad)=0;
};
//class circles
class circle : public geometry{
public:
    double x, y, r;
    circle(){};
    circle(double _x, double _y, double _r): x(_x), y(_y), r(_r){ngrad=2;}
    void evalConstr(const double *s, double *val, double *grad){
        double sx = s[0], sy = s[1];
        *val = pow(sx - x, 2) + pow(sy - y, 2) - r * r;
        grad[0] = 2 * (sx - x);
        grad[1] = 2 * (sy - y);
    }
};
//class sphere
class sphere : public geometry{
public:
    double x, y, z, r;
    sphere(){};
    sphere(double _x, double _y, double _z, double _r): x(_x), y(_y), z(_z), r(_r){ngrad=3;}
    void update(double _x, double _y, double _z, double _r){
        x = _x;
        y = _y;
        z = _z;
        r = _r;
    }
    void evalConstr(const double *s, double *val, double *grad){
        double sx = s[0], sy = s[1], sz = s[2];
        *val = -(pow(sx - x, 2) + pow(sy - y, 2) + pow(sz - z, 2)) + r * r;  // assume we want to be outside
        grad[0] = -2 * (sx - x);
        grad[1] = -2 * (sy - y);
        grad[2] = -2 * (sz - z);
        if(!out){
            *val *= -1.0;
            grad[0] *= -1.0;
            grad[1] *= -1.0;
            grad[2] *= -1.0;
        }
    }
    double evalConstr(const double *s){
        double sx = s[0], sy = s[1], sz = s[2];
        double val = -(pow(sx - x, 2) + pow(sy - y, 2) + pow(sz - z, 2)) + r * r;  // assume we want to be outside
        return val;
    }
    bool violate(const double *s, double tol=0){
        double val = -(pow(s[0] - x, 2) + pow(s[1] - y, 2) + pow(s[2] - z, 2)) + r * r;  // assume we want to be outside
        if((out && val < -tol) || (!out && val > tol)){
            return false;
        }
        else{
            return true;
        }
    }
    void lazyRepair(double *xin){  // repair a point by moving to the closest point on surface
        if(!violate(xin))
            return;
        VX dir(3);
        dir(0) = xin[0] - x;
        dir(1) = xin[1] - y;
        dir(2) = xin[2] - z;
        double dirnorm = dir.norm();
        if(dirnorm < 1e-6){
            dir.setRandom();
        }
        dir.array() /= dir.norm();
        xin[0] = x + r * dir(0);
        xin[1] = y + r * dir(1);
        xin[2] = z + r * dir(2);
    }
};
//class cylinder
class cylinder : public geometry{
public:
    double x, y, z, r;  // assume it is placed on the ground, so x, y, z and r should be sufficient to describe it, z is height
    cylinder(){};
    cylinder(double _x, double _y, double _z, double _r): x(_x), y(_y), z(_z), r(_r){ngrad=3;}
    void update(double _x, double _y, double _z, double _r){
        x = _x;
        y = _y;
        z = _z;
        r = _r;
    }
    void evalConstr(const double *s, double *val, double *grad){
        double sx = s[0], sy = s[1], sz = s[2];
        double val1 = -(pow(sx - x, 2) + pow(sy - y, 2)) + r * r;
        double val2 = z - sz;
        if(out){
            if(val1 > val2){
                *val = val2;
                grad[0] = 0;
                grad[1] = 0;
                grad[2] = -1.0;
            }
            else{
                *val = val1;
                grad[0] = -2 * (sx - x);
                grad[1] = -2 * (sy - y);
                grad[2] = 0;
            }
        }
        else {
            if(val1 > val2){
                *val = -val2;
                grad[0] = 0;
                grad[1] = 0;
                grad[2] = 1.0;
            }
            else{
                *val *= -val1;
                grad[0] = 2 * (sx - x);
                grad[1] = 2 * (sy - y);
                grad[2] = 0;
            }
        }
    }
};

//Declare the class for collision avoidance
class avoidConstraint : public TrajOpt::constraintFunctor {
public:
    std::vector<geometry*> obs;
    avoidConstraint(int N=1){
        dimstate = 12;
        dimctrl = 4;
        dimc = 1;
        nnz = 3 * N;
        cub = VX::Zero(N);
        clb = -1e10 * VX::Ones(N);
        obs.push_back(new sphere(0, 0, 0, 0));
    }
    void init(Json::Value &js){
        int N = js["radius"].size();
        for (int i = 0; i < N; i++) {
            double x = js["center"][i][0].asDouble(), y = js["center"][i][1].asDouble(), z = js["center"][i][2].asDouble();
            double r = js["radius"][i].asDouble();
            if(i == 0){
                ((sphere*)(obs[0]))->update(x, y, z, r);
            }
            else{
                obs.push_back(new sphere(x, y, z, r));
            }
        }
    }
    void update(const double *ob){
        ((sphere*)(obs[0]))->update(ob[0], ob[1], ob[2], ob[3]);
    }
    void eval(const double t, cRefV x, cRefV u, RefV f, SpMX &fx) {
        const double *s = x.data();
        double val;
        for (int i = 0; i < obs.size(); i++) {
            int ngrad = obs[i]->ngrad;
            VX Egrad(obs[i]->ngrad);
            double *grad = Egrad.data();
            obs[i]->evalConstr(s, &val, grad);
            f(i) = val;
            for(int j = 0; j < ngrad; j++){
                fx.coeffRef(i, j) = grad[j];
            }
        }
    }
    bool violate(const double *x, double tol=0){
        // determine if a given state in x violates constraints
        return ((sphere*)(obs[0]))->violate(x, tol);
    }
    void lazyRepair(double *x){
        return ((sphere*)(obs[0]))->lazyRepair(x);
    }
    double evalConstr(const double *x){
        return ((sphere*)(obs[0]))->evalConstr(x);
    }
    ~avoidConstraint(){
        for(auto i : obs)
            delete i;
    }
};

//Declare the class for collision avoidance
class cylinderConstraint : public TrajOpt::constraintFunctor {
public:
    std::vector<geometry*> obs;
    cylinderConstraint(int N=1){
        dimstate = 12;
        dimctrl = 4;
        dimc = 1;
        nnz = 3 * N;
        cub = VX::Zero(N);
        clb = -1e10 * VX::Ones(N);
        obs.push_back(new cylinder(0, 0, 0, 0));
    }
    void init(Json::Value &js){
        int N = js["radius"].size();
        for (int i = 0; i < N; i++) {
            double x = js["center"][i][0].asDouble(), y = js["center"][i][1].asDouble(), z = js["center"][i][2].asDouble();
            double r = js["radius"][i].asDouble();
            if(i == 0){
                ((cylinder*)(obs[0]))->update(x, y, z, r);
            }
            else{
                obs.push_back(new cylinder(x, y, z, r));
            }
        }
    }
    void update(const double *ob){
        ((cylinder*)(obs[0]))->update(ob[0], ob[1], ob[2], ob[3]);
    }
    void eval(const double t, cRefV x, cRefV u, RefV f, SpMX &fx) {
        const double *s = x.data();
        double val;
        for (int i = 0; i < obs.size(); i++) {
            int ngrad = obs[i]->ngrad;
            VX Egrad(obs[i]->ngrad);
            double *grad = Egrad.data();
            obs[i]->evalConstr(s, &val, grad);
            f(i) = val;
            for(int j = 0; j < ngrad; j++){
                fx.coeffRef(i, j) = grad[j];
            }
        }
    }
    ~cylinderConstraint(){
        for(auto i : obs)
            delete i;
    }
};
#endif //ROTORDB_UTILITYCLASS_H
