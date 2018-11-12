#ifndef CIRCLE_TRAJECTORY
#define CIRCLE_TRAJECTORY

#include <math.h>

//# define M_PI           3.14159265358979323846

class CircleTrajectory
{
public:
	double _centerX;   // x-coordinate of the center of the circle
	double _centerY;   // y-coordinate of the center of the circle
	double _centerZ;   // z-coordinate of the center of the circle
	double _radius;    // radius of the circle
	double _startAngle;
	int    _numberOfCircles; // number of circles
	double _angularVelocity;     // angularVelocity
	double _totalTime;
	double _maxLinearAcceleration;  // max linear acceleration allowd, for slow start mechanism

    // sort of a slow start mechanism, velocity will gradually increase, not abruptly.
	double rampUp(double t)
	{
		double maxLinearVelocity = _maxLinearAcceleration * (t - 0.0);
        
        // slow start
		if (maxLinearVelocity >= _angularVelocity * _radius)
			return _angularVelocity;
		// normal velocity
		else
			return maxLinearVelocity / _radius;
	}

public:
	CircleTrajectory()
	{
		_centerX = 0.0;
		_centerY = 0.0;
		_centerZ = 0.0;
		_radius = 0.0;
		_startAngle = 0.0;
		_numberOfCircles = 1;
		_angularVelocity = 0.0;
		_totalTime = 0.0;
		_maxLinearAcceleration = 0.0;

	}

	/**
	 * generate a circle trajectory in the X-Y plane
	 *
	 * @param startX          X-coordinate of start point
	 * @param startY          Y-coordinate of start point
	 * @param startZ          Z-coordinate of start point
	 * @param centerX         X-coordinate of center of the circle
	 * @param centerY         Y-coordinate of center of the circle
	 * @param numberOfCircles number of circles, at least 1
	 * @param time            time for every circle
	 */
	void generateTrajectory(double startX, double startY, double startZ, double centerX, double centerY, int numberOfCircles, double time, double maxLinearAcceleration)
	{
		_centerX = centerX;
		_centerY = centerY;
		_centerZ = startZ; // the circle is in the X-Y plane

		_radius = sqrt(pow((startX - centerX), 2) + pow((startY - centerY), 2));
		_startAngle = atan2((startY - centerY), (startX - centerX));

		if (numberOfCircles < 1)
		{
			numberOfCircles = 1;
		}
		_numberOfCircles = numberOfCircles;

		_angularVelocity = 2 * M_PI / time;
		_totalTime = time * numberOfCircles;
		_maxLinearAcceleration = maxLinearAcceleration;

	}
    
    /**
     * [samplePositionX description]
     * @param  t should start from 0.
     * @return   X-coordinate position at time t
     */
	double samplePositionX(double t)
	{
		double omega = rampUp(t);
		return _centerX + _radius * cos(t * omega + _startAngle);
	}

	double samplePositionY(double t)
	{
		double omega = rampUp(t);
		return _centerY + _radius * sin(t * omega + _startAngle);
	}

	double samplePositionZ(double t)
	{
		return _centerZ;
	}

	double sampleVelocityX(double t)
	{
		double omega = rampUp(t);
		return _radius * omega * -sin(t * omega + _startAngle);
	}

	double sampleVelocityY(double t)
	{
		double omega = rampUp(t);
		return _radius * omega * cos(t * omega + _startAngle);
	}

	double sampleVelocityZ(double t)
	{
		return 0.0;
	}

	double sampleAccelerationX(double t)
	{
		double omega = rampUp(t);
		return _radius * omega * omega * -cos(t * omega + _startAngle);
	}

	double sampleAccelerationY(double t)
	{
		double omega = rampUp(t);
		return _radius * omega * omega * -sin(t * omega + _startAngle);
	}

	double sampleAccelerationZ(double t)
	{
		return 0.0;
	}

	double getDuration()
	{
		return _totalTime;
	}

};

#endif /* CIRCLE_TRAJECTORY */