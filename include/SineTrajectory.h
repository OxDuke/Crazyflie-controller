#ifndef _SINETRAJECTORY_H_
#define _SINETRAJECTORY_H_ 

class SineTrajectory
{
private:
  char _axis;
  double _omega;
  double _duration;

  double _startX;
  double _startY;
  double _startZ;
  
public:
  SineTrajectory()
  {   
  }

  ~SineTrajectory()
  {}

  void generateTrajectory(double startX, double startY, double startZ, char axis, double omega, double duration)
  {
  	_startX = startX;
  	_startY = startY;
  	_startZ = startZ;

  	_axis = axis;
  	_omega = omega;
  	_duration = duration;

  }

  void sampleTrajectory(double t, double& Px, double& Py, double& Pz, double& Vx, double& Vy, double& Vz, double& Ax, double& Ay, double& Az, double& phi, double& theta, double& psi)
  {

  	double movingAxisPosition = cos(_omega * t) - 1;
  	double movingAxisVelocity = _omega * -sin(_omega * t);
  	double movingAxisAcceleration = _omega * _omega * -cos(_omega *t);

  	if(_axis == 'x')
  	{
  		Px = movingAxisPosition + _startX;
  		Vx = movingAxisVelocity;
  		Ax = movingAxisAcceleration;

  		Py = _startY;
  		Vy = Ay = 0.0;

  		Pz = _startZ;
  		Vz = Az = 0.0;

  	}
  	else if (_axis == 'y')
  	{
  		Py = movingAxisPosition + _startY;
  		Vy = movingAxisVelocity;
  		Az = movingAxisAcceleration;

  		Px = _startX;
  		Vx = Ax = 0.0;
  		
  		Pz = _startZ;
  		Vz = Az = 0.0;
  	}
  	else // _axis == 'z'
  	{
  		Px = _startX;
  		Vx = Ax = 0.0;
  		
  		Py = _startY;
  		Vy = Ay = 0.0;

  		Pz = movingAxisPosition + _startZ;
  		Vz = movingAxisVelocity;
  		Az = movingAxisAcceleration;
  	}
    
  }

  void getEndPosition(double& Px, double& Py, double& Pz)
  {
  	Px = _startX;
  	Py = _startY;
  	Pz = _startZ; 
  }

  double getEndTime()
  {
  	return _duration;
  }

  double omega()
  {
    return _omega;
  }

  char axis()
  {
    return _axis;
  }
    
};


#endif /* _SINETRAJECTORY_H_ */