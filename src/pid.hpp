#pragma once

#include <ros/ros.h>
#include <ros/console.h>

class PID
{
public:
    PID(
        double kp,
        double kd,
        double ki,
        double minOutput,
        double maxOutput,
        double integratorMin,
        double integratorMax,
        const std::string& name)
        : m_kp(kp)
        , m_kd(kd)
        , m_ki(ki)
        , m_minOutput(minOutput)
        , m_maxOutput(maxOutput)
        , m_integratorMin(integratorMin)
        , m_integratorMax(integratorMax)
        , m_integral(0)
        , m_previousError(0)
        , m_previousTime(ros::Time::now())
    {

    }

    void reset()
    {
        m_integral = 0;
        m_previousError = 0;
        m_previousTime = ros::Time::now();
    }

    void setIntegral(double integral)
    {
        m_integral = integral;
    }

    double ki() const
    {
        return m_ki;
    }

    double kp() const
    {
        return m_kp;
    }

    double kd() const
    {
        return m_kd;
    }

    double setKp(const double newK)
    {
        m_kp = newK;
    }

    double setKi(const double newK)
    {
        m_ki = newK;
    }

    double setKd(const double newK)
    {
        m_kd = newK;
    }

    double getIntegral()
    {
        return m_integral;
    }

    double update(double error)
    {
        ros::Time time = ros::Time::now();
        double dt = time.toSec() - m_previousTime.toSec();
        m_integral += error * dt;
        m_integral = std::max(std::min(m_integral, m_integratorMax), m_integratorMin);

        double p = m_kp * error;
        double d = 0;
        if (dt > 0)
        {
            d = m_kd * (error - m_previousError) / dt;
        }
        double i = m_ki * m_integral;

        double output = p + d + i;

        m_previousError = error;
        m_previousTime = time;
        return std::max(std::min(output, m_maxOutput), m_minOutput);
    }

    double update(double error, double dError) // dError: velocity term
    {
        ros::Time time = ros::Time::now();
        double dt = time.toSec() - m_previousTime.toSec();
        m_integral += error * dt;
        m_integral = std::max(std::min(m_integral, m_integratorMax), m_integratorMin);

        double p = m_kp * error;
        double d = m_kd * dError;
        double i = m_ki * m_integral;

        double output = p + d + i;

        m_previousError = error;
        m_previousTime = time;
        return std::max(std::min(output, m_maxOutput), m_minOutput);
    }

        double update(double error, double dError, double fuck) // dError: velocity term
    {
        ros::Time time = ros::Time::now();
        double dt = time.toSec() - m_previousTime.toSec();
        m_integral += error * dt;
        m_integral = std::max(std::min(m_integral, m_integratorMax), m_integratorMin);

        double p = m_kp * error;
        double d = m_kd * dError;
        double i = m_ki * m_integral;

        double output = p + d + i;

        m_previousError = error;
        m_previousTime = time;
        return std::max(std::min(output, m_maxOutput), m_minOutput);
    }

private:
    double m_kp;
    double m_kd;
    double m_ki;
    double m_minOutput;
    double m_maxOutput;
    double m_integratorMin;
    double m_integratorMax;
    double m_integral;
    double m_previousError;
    ros::Time m_previousTime;
};
