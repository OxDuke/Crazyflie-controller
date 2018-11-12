#ifndef MINIMUM_SNAP_TRAJECTORY
#define MINIMUM_SNAP_TRAJECTORY

#include <Eigen/Dense>
//#include <math.h>
#include <cmath>

//#include <iostream>

class MinimumSnapTrajectory
{
private:
  double _startPosition;
  double _endPosition;
  int _numberOfSegments;

  Eigen::VectorXd trajectory; //store the coefficients of the trajectory

  double* Times; //time needed for every segment of trajectory
  double* cumulativeTimes; //cumulative time taken for every segment of trajectory

  bool verifyFormat(int numberOfWayPoints)
  {
    if (numberOfWayPoints < 2)
      return false;
    return true;
  }

  void copyTs(double Ts[], int numberOfWayPoints)
  {
    Times = new double [numberOfWayPoints - 1];
    for (int i = 0; i < numberOfWayPoints - 1; ++i)
    {
      Times[i] = Ts[i];
    }
  }

  void calculateCumulativeTimes(double Ts[], int numberOfWayPoints)
  {
    cumulativeTimes = new double [numberOfWayPoints];
    cumulativeTimes[0] = 0;
    for (int i = 1; i < numberOfWayPoints; ++i)
    {
      cumulativeTimes[i] = cumulativeTimes[i - 1] + Ts[i - 1];
    }
  }

  void getSegmentAndTime(double t, int* tIndex, double* newT)
  {
    // t should be greater than 0 and less than the total time of trajectory
    if (t < 0.0)
      t = 0.0;
    else if (t > cumulativeTimes[_numberOfSegments])
      t = cumulativeTimes[_numberOfSegments];


    for (int i = 0; i < _numberOfSegments + 1; ++i)
    {
      if (t <= cumulativeTimes[i]) // figure out which segment does this t belong to
      {
        *tIndex = i;
        break;
      }
    }

    //assign the start point(t=0.0) to the first segment
    if (*tIndex == 0)
    {
      *tIndex = 1;
    }

    // scale t to 0-1
    *newT = (t - cumulativeTimes[*tIndex - 1]) / (Times[*tIndex - 1]);

  }

public:
  MinimumSnapTrajectory()
  {
    _startPosition = 0.0;
    _endPosition = 0.0;
    _numberOfSegments = 0;

  }

  /**
   * [generateTrajectory description]
   * @param wayPoints         an array of way points
   * @param numberOfWayPoints numberOfWayPoints, at least 2
   * @param Ts                time needed for every segment(between two way points)
   * @note                    the length of the wayPoints should be larger than Ts by exactly 1
   */
  void generateTrajectory(double wayPoints[], int numberOfWayPoints, double Ts[])
  {
    // verify the format of input
    if (verifyFormat(numberOfWayPoints) == false)
      return;
    _numberOfSegments = numberOfWayPoints - 1;
    int n = _numberOfSegments;

    // step1. miscellenous works
    copyTs(Ts, numberOfWayPoints);
    calculateCumulativeTimes(Ts, numberOfWayPoints);
    _startPosition = wayPoints[0];
    _endPosition = wayPoints[numberOfWayPoints - 1];

    // step2. solve trajectory coefficients
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(4 + 8 * (n - 1) + 4, 8 * n);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(4 + (2 + 6) * (n - 1) + 4);


    //construct matrix A
    A.topLeftCorner(4, 8) <<
                          0, 0, 0, 0, 0, 0, 0, 1,
                          0, 0, 0, 0, 0, 0, 1, 0,
                          0, 0, 0, 0, 0, 2, 0, 0,
                          0, 0, 0, 0, 6, 0, 0, 0;

    A.bottomRightCorner(4, 8) <<
                              1, 1, 1, 1, 1, 1, 1, 1,
                              7, 6, 5, 4, 3, 2, 1, 0,
                              42, 30, 20, 12, 6, 2, 0, 0,
                              210, 120, 60, 24, 6, 0, 0, 0;

    for (int i = 0; i < n - 1; ++i)
    {
      A.block(4 + 8 * i, 8 * i, 8, 8) << 1, 1, 1, 1, 1, 1, 1, 1,
              0, 0, 0, 0, 0, 0, 0, 0,
              7, 6, 5, 4, 3, 2, 1, 0,
              42, 30, 20, 12, 6, 2, 0, 0,
              210, 120, 60, 24, 6, 0, 0, 0,
              840, 360, 120, 24, 0, 0, 0, 0,
              2520, 720, 120, 0, 0, 0, 0, 0,
              5040, 720, 0, 0, 0, 0, 0, 0;

      A.block(4 + 8 * i, 8 + 8 * i, 8, 8) << 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 1,
              0, 0, 0, 0, 0, 0, -1, 0,
              0, 0, 0, 0, 0, -2, 0, 0,
              0, 0, 0, 0, -6, 0, 0, 0,
              0, 0, 0, -24, 0, 0, 0, 0,
              0, 0, -120, 0, 0, 0, 0, 0,
              0, -720, 0, 0, 0, 0, 0, 0;
    }

    // construct vector b
    b.head(1) << wayPoints[0];
    b.tail(4) << wayPoints[n], 0.0, 0.0, 0.0;
    for (int i = 0; i < n - 1; ++i)
    {
      b.segment(4 + 8 * i, 2) << wayPoints[i + 1], wayPoints[i + 1];
    }

    // time to solve coefficients!
    trajectory = A.fullPivLu().solve(b);

  }

  /**
   * [generateTrajectory description]
   * @param wayPoints         an array of way points
   * @param numberOfWayPoints numberOfWayPoints, at least 2
   * @param Ts                time needed for every segment(between two way points)
   * @note                    the length of the wayPoints should be larger than Ts by exactly 1
   */
  void generateTrajectory(double wayPoints[], double startVelocity, double endVelocity, int numberOfWayPoints, double Ts[])
  {
    // verify the format of input
    if (verifyFormat(numberOfWayPoints) == false)
      return;
    _numberOfSegments = numberOfWayPoints - 1;
    int n = _numberOfSegments;

    // step1. miscellenous works
    copyTs(Ts, numberOfWayPoints);
    calculateCumulativeTimes(Ts, numberOfWayPoints);
    _startPosition = wayPoints[0];
    _endPosition = wayPoints[numberOfWayPoints - 1];

    // step2. solve trajectory coefficients
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(4 + 8 * (n - 1) + 4, 8 * n);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(4 + (2 + 6) * (n - 1) + 4);


    //construct matrix A
    A.topLeftCorner(4, 8) <<
                          0, 0, 0, 0, 0, 0, 0, 1,
                          0, 0, 0, 0, 0, 0, 1, 0,
                          0, 0, 0, 0, 0, 2, 0, 0,
                          0, 0, 0, 0, 6, 0, 0, 0;

    A.bottomRightCorner(4, 8) <<
                              1, 1, 1, 1, 1, 1, 1, 1,
                              7, 6, 5, 4, 3, 2, 1, 0,
                              42, 30, 20, 12, 6, 2, 0, 0,
                              210, 120, 60, 24, 6, 0, 0, 0;

    for (int i = 0; i < n - 1; ++i)
    {
      A.block(4 + 8 * i, 8 * i, 8, 8) << 1, 1, 1, 1, 1, 1, 1, 1,
              0, 0, 0, 0, 0, 0, 0, 0,
              7, 6, 5, 4, 3, 2, 1, 0,
              42, 30, 20, 12, 6, 2, 0, 0,
              210, 120, 60, 24, 6, 0, 0, 0,
              840, 360, 120, 24, 0, 0, 0, 0,
              2520, 720, 120, 0, 0, 0, 0, 0,
              5040, 720, 0, 0, 0, 0, 0, 0;

      A.block(4 + 8 * i, 8 + 8 * i, 8, 8) << 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 1,
              0, 0, 0, 0, 0, 0, -1, 0,
              0, 0, 0, 0, 0, -2, 0, 0,
              0, 0, 0, 0, -6, 0, 0, 0,
              0, 0, 0, -24, 0, 0, 0, 0,
              0, 0, -120, 0, 0, 0, 0, 0,
              0, -720, 0, 0, 0, 0, 0, 0;
    }

    // construct vector b
    b.head(2) << wayPoints[0], startVelocity;
    b.tail(4) << wayPoints[n], endVelocity, 0.0, 0.0;
    for (int i = 0; i < n - 1; ++i)
    {
      b.segment(4 + 8 * i, 2) << wayPoints[i + 1], wayPoints[i + 1];
    }

    // time to solve coefficients!
    trajectory = A.fullPivLu().solve(b);

  }

  double samplePosition(double t)
  {
    int tIndex;
    double newT;

    getSegmentAndTime(t, &tIndex, &newT);

    Eigen::VectorXd ts(8);
    ts << pow(newT, 7), pow(newT, 6), pow(newT, 5), pow(newT, 4), pow(newT, 3), pow(newT, 2), newT, 1.0;
    return ts.dot(trajectory.segment(8 * (tIndex - 1), 8));
  }

  double sampleVelocity(double t)
  {
    int tIndex;
    double newT;

    getSegmentAndTime(t, &tIndex, &newT);

    Eigen::VectorXd ts(8);
    ts << 7 * pow(newT, 6), 6 * pow(newT, 5), 5 * pow(newT, 4), 4 * pow(newT, 3), 3 * pow(newT, 2), 2 * newT, 1.0, 0.0;
    return ts.dot(trajectory.segment(8 * (tIndex - 1), 8)) / (Times[tIndex - 1]);
  }

  double sampleAcceleration(double t)
  {
    int tIndex;
    double newT;

    getSegmentAndTime(t, &tIndex, &newT);

    Eigen::VectorXd ts(8);
    ts << 42 * pow(newT, 5), 30 * pow(newT, 4), 20 * pow(newT, 3), 12 * pow(newT, 2), 6 * newT, 2.0, 0.0, 0.0;
    return ts.dot(trajectory.segment(8 * (tIndex - 1), 8)) / (Times[tIndex - 1]) / (Times[tIndex - 1]);
  }

  double getStartPosition()
  {
    return _startPosition;
  }

  double getEndPosition()
  {
    return _endPosition;
  }

  double getEndTime()
  {
    return cumulativeTimes[_numberOfSegments];
  }

};

class TDMST
{
private:
  MinimumSnapTrajectory trjX, trjY, trjZ;

public:
  TDMST()
  {
  }

  ~TDMST()
  {}

  void generateTrajectory(double wayPointsX[], double wayPointsY[], double wayPointsZ[], int numberOfWayPoints, double averageVelocity)
  {
    double Ts[numberOfWayPoints - 1];
    for (int i = 0; i < numberOfWayPoints - 1; ++i)
    {
      double deltaX = wayPointsX[i + 1] - wayPointsX[i];
      double deltaY = wayPointsY[i + 1] - wayPointsY[i];
      double deltaZ = wayPointsZ[i + 1] - wayPointsZ[i];
      Ts[i] = sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ) / averageVelocity;
    }

    trjX.generateTrajectory(wayPointsX, numberOfWayPoints, Ts);
    trjY.generateTrajectory(wayPointsY, numberOfWayPoints, Ts);
    trjZ.generateTrajectory(wayPointsZ, numberOfWayPoints, Ts);
  }

  void generateTrajectory(double wayPointsX[], double wayPointsY[], double wayPointsZ[],
                          double startVelocityX, double startVelocityY, double startVelocityZ,
                          double endVelocityX, double endVelocityY, double endVelocityZ,
                          int numberOfWayPoints, double averageVelocity)
  {
    double Ts[numberOfWayPoints - 1];
    for (int i = 0; i < numberOfWayPoints - 1; ++i)
    {
      double deltaX = wayPointsX[i + 1] - wayPointsX[i];
      double deltaY = wayPointsY[i + 1] - wayPointsY[i];
      double deltaZ = wayPointsZ[i + 1] - wayPointsZ[i];
      Ts[i] = sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ) / averageVelocity;
    }

    trjX.generateTrajectory(wayPointsX, startVelocityX, endVelocityX, numberOfWayPoints, Ts);
    trjY.generateTrajectory(wayPointsY, startVelocityY, endVelocityY, numberOfWayPoints, Ts);
    trjZ.generateTrajectory(wayPointsZ, startVelocityZ, endVelocityZ, numberOfWayPoints, Ts);
  }

  void sampleTrajectory(double t, double& Px, double& Py, double& Pz, double& Vx, double& Vy, double& Vz, double& Ax, double& Ay, double& Az, double& phi, double& theta, double& psi)
  {
    Px = trjX.samplePosition(t);
    Py = trjY.samplePosition(t);
    Pz = trjZ.samplePosition(t);

    Vx = trjX.sampleVelocity(t);
    Vy = trjY.sampleVelocity(t);
    Vz = trjZ.sampleVelocity(t);

    Ax = trjX.sampleAcceleration(t);
    Ay = trjY.sampleAcceleration(t);
    Az = trjZ.sampleAcceleration(t);
  }

  void getEndPosition(double& Px, double& Py, double& Pz)
  {
    double endTime = getEndTime();

    Px = trjX.samplePosition(endTime);
    Py = trjY.samplePosition(endTime);
    Pz = trjZ.samplePosition(endTime);
  }

  double getEndTime()
  {
    return trjX.getEndTime();
  }

  double getDuration()
  {
    return trjX.getEndTime();
  }


};


#endif /* MINIMUM_SNAP_TRAJECTORY */