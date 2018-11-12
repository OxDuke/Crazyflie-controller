// reading a text file
//#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

#include "TigerEigen.h"

class Gao
{
private:
  std::vector<double> _cumT; // cumulative T

  // positions
  std::vector<double> _Px;
  std::vector<double> _Py;
  std::vector<double> _Pz;

  // velocities
  std::vector<double> _Vx;
  std::vector<double> _Vy;
  std::vector<double> _Vz;

  // accelerations
  std::vector<double> _Ax;
  std::vector<double> _Ay;
  std::vector<double> _Az;

  // euler angles
  std::vector<double> _phi;
  std::vector<double> _theta;
  std::vector<double> _psi;

  // offsets
  double _offsetX;
  double _offsetY;
  double _offsetZ;

  double _timeOffset;


  /**
   * [whichSegment description]
   * @param t          time
   * @param segment    return the segment, from 1 to ...
   * @param timeScaled return scaled(0.0-1.0) time 
   */
  void whichSegment(double t, int* segment, double* timeScaled)
  {
    t = t - _timeOffset;

    if(t < 0.0)
      t = 0.0;
    else if (t > _cumT.back())
      t = _cumT.back();

    std::vector<double>::iterator low;
    low = std::lower_bound(_cumT.begin(), _cumT.end(), t);

    if (low - _cumT.begin() == 0)
    {
      *segment = 1;
      *timeScaled = 0.0;
    }
    else
    {
      *segment = low - _cumT.begin();
      *timeScaled = (t - *(low - 1))/(*low - *(low - 1));
    }

  }

public:
  Gao():
  _offsetX(0.0),
  _offsetY(0.0),
  _offsetZ(0.0),
  _timeOffset(0.0)
  {
    _cumT.clear();

  }

  void clearTrajectory()
  {
    _cumT.clear();

    _Px.clear();
    _Py.clear();
    _Pz.clear();

    _Vx.clear();
    _Vy.clear();
    _Vz.clear();

    _Ax.clear();
    _Ay.clear();
    _Az.clear();

    _phi.clear();
    _theta.clear();
    _psi.clear();

    _offsetX = 0.0;
    _offsetY = 0.0;
    _offsetZ = 0.0;

    _timeOffset = 0.0;

  }

  void generateTrajectory(std::string fileName)
  {
    std::ifstream myFile(fileName.c_str());
    std::string line;

    clearTrajectory();

    if (myFile.is_open())
    {
      while (getline(myFile, line))
      {
        std::stringstream lineStream(line);
        double t;
        double x, y, z;
        double phi, theta, psi;
        double vx, vy, vz;
        double p, q, r;
        double ax, ay, az;
        double u1, u2 ,u3, u4;

        lineStream >> t >> x >> y >> z >> phi >> theta >> psi >> vx >> vy >> vz >> p >> q >> r >> u1 >> u2 >> u3 >> u4 >> ax >> ay >> az;

        _cumT.push_back(t);

        _Px.push_back(x);
        _Py.push_back(y);
        _Pz.push_back(z);

        _Vx.push_back(vx);
        _Vy.push_back(vy);
        _Vz.push_back(vz);

        _Ax.push_back(ax);
        _Ay.push_back(ay);
        _Az.push_back(az);

        _phi.push_back(phi);
        _theta.push_back(theta);
        _psi.push_back(psi);

      }
    }
    else
    {
      std::cout << "Unable to open file, trajectory generation failed!" << std::endl;
    }
  }

  void generateTrajectory(RefV vec, RefM acce, double tfin=0)
  {
    std::cout << "acce = " << acce << std::endl;
    int size = vec.size();
    int N = 20, dimx = 12, dimu = 4;
    double tf = 0;
    if(tfin == 0){
      N = (size - 1 + dimu) / (dimx + dimu);
      tf = vec(N*dimx + (N - 1)*dimu);
    }
    else{
      N = (size + dimu) / (dimx + dimu);
      tf = tfin;
    }
    double dt = tf / (N - 1);

    clearTrajectory();

    for(int i = 0; i < N; i++){
      int addind = i * dimx;
      _cumT.push_back(i * dt);
      _Px.push_back(vec(addind));
      _Py.push_back(vec(addind + 1));
      _Pz.push_back(vec(addind + 2));
      _Vx.push_back(vec(addind + 6));
      _Vy.push_back(vec(addind + 7));
      _Vz.push_back(vec(addind + 8));
      _phi.push_back(vec(addind + 3));
      _theta.push_back(vec(addind + 4));
      _psi.push_back(vec(addind + 5));
      _Ax.push_back(acce(0, i));
      _Ay.push_back(acce(1, i));
      _Az.push_back(acce(2, i));
    }

  }

  void setPositionOffset(double offsetX, double offsetY, double offsetZ)
  {
    _offsetX = offsetX;
    _offsetY = offsetY;
    _offsetZ = offsetZ;
  }

  void setTimeOffset(double timeOffset)
  {
    _timeOffset = timeOffset;
  }

  double samplePositionX(double t)
  {
    int segment;
    double timeScaled;
    whichSegment(t, &segment, &timeScaled);

    return _Px[segment - 1] + timeScaled * (_Px[segment] - _Px[segment - 1]);
  }

  double samplePositionY(double t)
  {
    int segment;
    double timeScaled;
    whichSegment(t, &segment, &timeScaled);

    return _Py[segment - 1] + timeScaled * (_Py[segment] - _Py[segment - 1]);
  }

  double samplePositionZ(double t)
  {
    int segment;
    double timeScaled;
    whichSegment(t, &segment, &timeScaled);

    return _Pz[segment - 1] + timeScaled * (_Pz[segment] - _Pz[segment - 1]);
  }

  double sampleVelocityX(double t)
  {
    int segment;
    double timeScaled;
    whichSegment(t, &segment, &timeScaled);

    return _Vx[segment - 1] + timeScaled * (_Vx[segment] - _Vx[segment - 1]);
  }

  double sampleVelocityY(double t)
  {
    int segment;
    double timeScaled;
    whichSegment(t, &segment, &timeScaled);

    return _Vy[segment - 1] + timeScaled * (_Vy[segment] - _Vy[segment - 1]);
  }

  double sampleVelocityZ(double t)
  {
    int segment;
    double timeScaled;
    whichSegment(t, &segment, &timeScaled);

    return _Vz[segment - 1] + timeScaled * (_Vz[segment] - _Vz[segment - 1]);
  }

  double sampleAccelerationX(double t)
  {
    int segment;
    double timeScaled;
    whichSegment(t, &segment, &timeScaled);

    return _Ax[segment - 1] + timeScaled * (_Ax[segment] - _Ax[segment - 1]);
  }

  double sampleAccelerationY(double t)
  {
    int segment;
    double timeScaled;
    whichSegment(t, &segment, &timeScaled);

    return _Ay[segment - 1] + timeScaled * (_Ay[segment] - _Ay[segment - 1]);
  }

  double sampleAccelerationZ(double t)
  {
    int segment;
    double timeScaled;
    whichSegment(t, &segment, &timeScaled);

    return _Az[segment - 1] + timeScaled * (_Az[segment] - _Az[segment - 1]);
  }

  void sampleTrajectory(double t, double& Px, double& Py, double& Pz, double& Vx, double& Vy, double& Vz, double& Ax, double& Ay, double& Az, double& phi, double& theta, double& psi)
  {
    int segment;
    double timeScaled;
    whichSegment(t, &segment, &timeScaled);

    Px = _Px[segment - 1] + timeScaled * (_Px[segment] - _Px[segment - 1]) + _offsetX;
    Py = _Py[segment - 1] + timeScaled * (_Py[segment] - _Py[segment - 1]) + _offsetY;
    Pz = _Pz[segment - 1] + timeScaled * (_Pz[segment] - _Pz[segment - 1]) + _offsetZ;

    Vx = _Vx[segment - 1] + timeScaled * (_Vx[segment] - _Vx[segment - 1]);
    Vy = _Vy[segment - 1] + timeScaled * (_Vy[segment] - _Vy[segment - 1]);
    Vz = _Vz[segment - 1] + timeScaled * (_Vz[segment] - _Vz[segment - 1]);

    Ax = _Ax[segment - 1] + timeScaled * (_Ax[segment] - _Ax[segment - 1]);
    Ay = _Ay[segment - 1] + timeScaled * (_Ay[segment] - _Ay[segment - 1]);
    Az = _Az[segment - 1] + timeScaled * (_Az[segment] - _Az[segment - 1]);

    phi = _phi[segment - 1] + timeScaled * (_phi[segment] - _phi[segment - 1]);
    theta = _theta[segment - 1] + timeScaled * (_theta[segment] - _theta[segment - 1]);
    psi = _psi[segment - 1] + timeScaled * (_psi[segment] - _psi[segment - 1]);
    
  }

  double getDuration()
  {
    return _cumT.back();
  }

  double getEndTime()
  {
    return _cumT.back() + _timeOffset;
  }

  double getEndPosition(double& Px, double& Py, double &Pz)
  {
    Px = _Px.back() + _offsetX;
    Py = _Py.back() + _offsetY;
    Pz = _Pz.back() + _offsetZ;
  }

  void printVector(std::vector<double>& v)
  {
    for (auto vi : v)
      std::cout << vi << ' ';
    std::cout << '\n' << std::endl;
  }

};

