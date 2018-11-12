#ifndef MINIMUM_JERK_TRAJECTORY
#define MINIMUM_JERK_TRAJECTORY

#include <Eigen/Dense>
#include <math.h>

class MinimumJerkTrajectory{
private:
  // three constrains that define a min jerk trajectory
  double m_startPosition;
  double m_endPosition;
  double m_duration;

  // conveninet vector
  double C[6];

public:
  MinimumJerkTrajectory()
  {
    m_startPosition = 0;
    m_endPosition = 0;
    m_duration = 0;
  }

  MinimumJerkTrajectory(double startPosition, double endPosition, double duration)
  {
    generateTrajectory(startPosition, endPosition, duration);
  }

  void generateTrajectory(double startPosition, double endPosition, double duration)
  {

    m_startPosition = startPosition;
    m_endPosition = endPosition;
    m_duration = duration;

    double T[6];
    T[5] = pow(duration,5);
    T[4] = pow(duration,4);
    T[3] = pow(duration,3);
    T[2] = pow(duration,2);
    T[1] = duration;
    T[0] = 1;

    // solve Ax=b
    Eigen::Matrix<double, 6, 6> A;
    A << 0, 0, 0, 0, 0, 1,
         T[5], T[4], T[3], T[2], T[1], 1,
         0, 0, 0, 0, 1, 0,
         5*T[4], 4*T[3], 3*T[2], 2*T[1], 1, 0,
         0, 0, 0, 2, 0, 0,
         20*T[3], 12*T[2], 6*T[1], 2, 0, 0;

    Eigen::VectorXd b(6);
    b << startPosition, endPosition, 0, 0, 0, 0;

    Eigen::VectorXd x(6);

    x = A.fullPivLu().solve(b);
    C[5] = x[0];
    C[4] = x[1];
    C[3] = x[2];
    C[2] = x[3];
    C[1] = x[4];
    C[0] = x[5];
  }

  double samplePosition(double time)
  {
    return C[5]*pow(time,5) + C[4]*pow(time,4) + C[3]*pow(time,3) + C[2]*time*time + C[1]*time + C[0]*1;
  }

  double sampleVelocity(double time)
  {
    return 5.0*C[5]*pow(time,4) + 4.0*C[4]*pow(time,3) + 3.0*C[3]*time*time + 2.0*C[2]*time + C[1];
  }

  double sampleAcceleration(double time)
  {
    return 20.0*C[5]*pow(time,3) + 12.0*C[4]*time*time + 6.0*C[3]*time + 2.0*C[2];
  }

  double getStartPosiiton()
  {
    return m_startPosition;
  }

  double getEndPosition()
  {
    return m_endPosition;
  }

  double getDuration()
  {
    return m_duration;
  }

};


#endif /* MINIMUM_JERK_TRAJECTORY */
